/*
 Author: Shohei Kojima @ RIKEN
 Description:
    This reads repeat fasta and saves 11-mers (as uint32_t)
    that appears at least 1 time in a file
    as a 2-bit binary file.
    If a non-ATGC character appears, it converts
    such character to all possible nucleotides.
 Compile:
    g++ -o convert_rep_to_2bit_k11 convert_rep_to_2bit_k11.cpp -O2
    g++ -shared -fPIC -o convert_rep_to_2bit_k11.so convert_rep_to_2bit_k11.cpp -O2
 Usage:
    ./prog input.fa output_prefix
 Output:
    1) output_prefix.mk : 2-bit compressed k-mers
    2) output_prefix.mi : number of k-mers saved in output_prefix.mk
 Misc info:
    This only uses one thread.
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include "dna_to_2bit.hpp"
#include "convert_rep_to_2bit_k11.hpp"
using namespace dna_to_2bit_hpp;

#define REP_KMER_SIZE 11
#define SHIFT_16_TO_11 10

typedef unsigned long long ull;
const char FASTA_HEADER_START = '>';


/*
 Check whether the sequence contains non-ATGC character(s)
 */
inline bool contains_nonATGC(std::string& seq) {
    for (char c : seq) {
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
            return true;
        }
    }
    return false;
}


/*
 Convert non-ATGC characters to any one char.
 */
inline const char* convert_nonATGC_to_any(std::string& seq, char any) {
    std::string converted;
    for (char c : seq) {
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
            converted += any;
        } else {
            converted += c;
        }
    }
    return converted.c_str();
}


/*
 This converts 11-nt DNA to 2bit and stores as uint32_t in the vector.
 Args:
    1) pointer to seq
    2) seq len (ull)
    3) window_size (must be 16)
    4) vector to store the 2bit seq
 */
inline void dna_to_2bit_bidirectional_22(char* seq, ull& seqlen, const int& window_size, std::vector<uint32_t>& v) {
    if (seqlen < REP_KMER_SIZE) {
        return;
    }
    
    // first window_size
    uint32_t bit2f=0;
    uint32_t bit2r=0;
    int nn=0;
    for (int i=0; i < SHIFT_16_TO_11; i++) {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_32[seq[i]];
        bit2r >>= 2;
        bit2r |= dna_to_2bitr_32[seq[i]];
        if (seq[i] == 'N' || seq[i] == 'n') {  // ignore when N or n appears
            nn= REP_KMER_SIZE - 1;
        } else if (nn > 0) {  // within window_size-nt from N or n
            nn--;
        }
    }
    bit2r >>= SHIFT_16_TO_11;
    if (nn == 0) {
        if (bit2f == bit2r) {
            v.push_back(bit2f);
        } else {
            v.push_back(bit2f);
            v.push_back(bit2r);
        }
    }
    
    // rolling calc.
    for (ull i=REP_KMER_SIZE; i < seqlen; i++) {
        // shift back
        bit2r <<= SHIFT_16_TO_11;
        // usual processing
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_32[seq[i]];
        bit2r >>= 2;
        bit2r |= dna_to_2bitr_32[seq[i]];
        // shift forward
        bit2f <<= SHIFT_16_TO_11;
        bit2f >>= SHIFT_16_TO_11;
        bit2r >>= SHIFT_16_TO_11;
        // usual processing
        if (seq[i] == 'N' || seq[i] == 'n') {  // ignore when N or n appears
            nn= REP_KMER_SIZE - 1;
        } else if (nn > 0) {  // within window_size-nt from N or n
            nn--;
        } else if (bit2f == bit2r) {
            v.push_back(bit2f);
        } else {
            v.push_back(bit2f);
            v.push_back(bit2r);
        }
    }
}


/*
 This is the de facto main func.
 Args:
    1) fasta file (e.g., human repeat consensus fasta file)
    2) output file (output_prefix.mk). If file is already present, this overwrites it.
    3) output file (output_prefix.mi). If file is already present, this overwrites it.
 */
int find_and_save_all_kmers(char* fa, char* out_mk, char* out_mi) {
    // for 2bit conversion
    const int window_size=init_dna_to_2bit_32();
    
    // open fasta
    std::ifstream infile;
    infile.open(fa);
    if (! infile.is_open()) { return 1; }
    
    // read fasta
    std::vector<const char*> seqs;
    std::vector<ull> seqlens;
    std::string line;
    std::string seq;
    while (std::getline(infile, line)) {
        if (line[0] == FASTA_HEADER_START) {
            if (! seq.empty()) {
                ull seq_size=seq.size();
                if (seq_size >= REP_KMER_SIZE) {
                    if (contains_nonATGC(seq) == true) {  // convert non-ATGC letter to all nucleoside
                        seqs.push_back(convert_nonATGC_to_any(seq, 'A'));
                        seqs.push_back(convert_nonATGC_to_any(seq, 'T'));
                        seqs.push_back(convert_nonATGC_to_any(seq, 'G'));
                        seqs.push_back(convert_nonATGC_to_any(seq, 'C'));
                        for (int i=0; i < 4; i++) {
                            seqlens.push_back(seq_size);
                        }
                    } else {
                        seqs.push_back(seq.c_str());  // save in the vector
                        seqlens.push_back(seq_size);
                    }
                }
                seq.clear();  // clear for next fasta entry
            }
        } else {
            transform(line.begin(), line.end(), line.begin(),
                [](unsigned char c) {return toupper(c);}); // convert to upper letters
            seq += line;
        }
    }
    if (seq.size() >= REP_KMER_SIZE) {
        if (contains_nonATGC(seq) == true) {  // convert non-ATGC letter to all nucleoside
            seqs.push_back(convert_nonATGC_to_any(seq, 'A'));
            seqs.push_back(convert_nonATGC_to_any(seq, 'T'));
            seqs.push_back(convert_nonATGC_to_any(seq, 'G'));
            seqs.push_back(convert_nonATGC_to_any(seq, 'C'));
            for (int i=0; i < 4; i++) {
                seqlens.push_back(seq.size());
            }
        } else {
            seqs.push_back(seq.c_str());  // save in the vector
            seqlens.push_back(seq.size());
        }
    }
    
    // 2bit conversion
    std::vector<uint32_t> v;
    ull seqs_size=seqs.size();
    for (ull i=0; i < seqs_size; i++) {
        dna_to_2bit_bidirectional_22((char*)seqs[i], seqlens[i], window_size, v);
    }
    
    // sort
    std::cout << "Sorting vector... (k-mers loaded = " << v.size() << ")" << std::endl;
    std::sort(v.begin(), v.end());
    std::cout << "Sorting finished" << std::endl;
    
    // save k-mers in output_prefix.mk
    std::ofstream outfile;
    outfile.open(out_mk, std::ios::binary);
    if (! outfile.is_open()) {
        return 1;
    }
    
    ull num_saved=0;
    ull veclen=v.size();
    size_t type_size=sizeof(uint32_t);
    for (ull i=1; i < veclen; i++) {
        if (! (v[i - 1] == v[i])) {
            outfile.write((char*)&v[i - 1], type_size);
            num_saved++;
        }
    }
    if (! (v[veclen - 2] == v[veclen - 1])) {
        outfile.write((char*)&v[veclen - 1], type_size);
        num_saved++;
    }
    std::cout << num_saved << " k-mers found" << std::endl;
    
    // save redundant k-mers in output_prefix.mi
    std::ofstream outfile2;
    outfile2.open(out_mi);
    if (! outfile.is_open()) {
        return 1;
    }
    outfile2 << num_saved << std::endl;
    
    return 0;
}


/*
 Just a file existence checker
 */
bool file_checker(char* f_path) {
    std::FILE* infile=fopen(f_path, "r");
    if (! infile) { return 1; }
    fclose(infile);
    return 0;
}


/*
 This is a wrapper of find_and_save_all_kmers().
 This checks the existence of files.
 Usage:
    ./prog input.fa output_prefix
 Args:
    1) input.fa
    2) output_prefix
 */
int main(int argc, char* argv[]) {
    // argv check
    if (argc <= 2) {
        std::cerr << "Please specify a fasta file and output prefix." << std::endl;
        return 1;
    } else if (argc >= 4) {
        std::cerr << "Too many arguments. Please specify a fasta file and output prefix." << std::endl;
        return 1;
    }
    
    // summarize paths
    char* fa=argv[1];
        
    std::string out_prefix=std::string(argv[2]);
    std::string str_mk= out_prefix + ".mk";
    std::string str_mi= out_prefix + ".mi";
    char* out_mk=(char *)str_mk.c_str();
    char* out_mi=(char *)str_mi.c_str();
    
    // file check
    assert(file_checker(fa) == 0);
    std::ofstream outfile;
    outfile.open(out_mk, std::ios::binary);
    if (! outfile.is_open()) { return 1; }
    outfile.close();
    
    std::cout << "in.fa: " << fa << std::endl;
    std::cout << "out.mk: " << out_mk << std::endl;
    std::cout << "out.mi: " << out_mi << std::endl;
    
    // process find_and_save_all_kmers()
    int ret=find_and_save_all_kmers(fa, out_mk, out_mi);
    
    return ret;
}

