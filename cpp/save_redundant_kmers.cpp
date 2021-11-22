/*
 Author: Shohei Kojima @ RIKEN
 Description:
    This reads genome fasta and saves 32-mers
    that appears more than 1 time in a genome
    as a 2-bit binary file.
 Compile:
    g++ -o save_redundant_kmers save_redundant_kmers.cpp -O2
    g++ -shared -fPIC -o save_redundant_kmers.so save_redundant_kmers.cpp -O2
 Usage:
    ./prog input.fa output_prefix
 Output:
    1) output_prefix.mk : 2-bit compressed k-mers
    2) output_prefix.mi : number of k-mers saved in output_prefix.mk
 Prerequisites:
    fasta index for the input fasta file.
 Misc info:
    This only uses one thread.
    In the case of GRCh38,
        this requires 48 GB RAM;
        sorts 6,086,521,103 k-mers;
        takes ~10 min to sort;
        there are 241,273,035 redundant k-mers;
        output .mk file is 1.8GB.
 */

#include <fstream>
#include <iostream>
#include <vector>
#include "dna_to_2bit.hpp"
#include "parse_fai.hpp"
#include "save_redundant_kmers.hpp"
using namespace dna_to_2bit_hpp;
using namespace parse_fai_hpp;


/*
 This is the de facto main func.
 Args:
    1) fasta file (e.g., human genome)
    2) fasta index (must be indexed by `samtools faidx` command)
    3) output file (output_prefix.mk). If file is already present, this overwrites it.
    4) output file (output_prefix.mi). If file is already present, this overwrites it.
 */
int find_and_save_red_kmers(char* fa, char* f_fai, char* out_mk, char* out_mi) {
    // read fai
    fai_parser fai;
    int chr_num=fai.load_fai(f_fai);
    ull genome_total_len=0;
    for (int chr=0; chr < chr_num; chr++) {
        fai_info* chr_info=fai.get_chr_info(chr);
        genome_total_len += chr_info->chr_len;
    }
    std::cout << "Genome length: " << genome_total_len << std::endl;
    
    // open fa
    std::FILE* infile=fopen(fa, "r");
    if (! infile) {
        return 1;
    }
    
    // for 2bit conversion
    const int window_size=init_dna_to_2bit_64();
    std::vector<uint64_t> v;
    v.reserve(genome_total_len);
    
    // read chrs
    for (int chr=0; chr < chr_num; chr++) {
        fai_info* chr_info=fai.get_chr_info(chr);
        ull chr_len =chr_info->chr_len;
        ull buf_size=chr_info->byte_per_line;
        char buf[buf_size];
        
        // read chr
        char* seq= new char[chr_info->chr_len + 1];
        fseeko64(infile, chr_info->seq_start, SEEK_SET);
        ull n_char_read=0;
        while (fread(buf, sizeof(char), buf_size, infile)) {  // read one line
            for (int i=0; i + 1 < buf_size; i++) {
                if (n_char_read >= chr_len) {
                    break;
                }
                seq[n_char_read]=buf[i];
                n_char_read++;
            }
        }
        
        // 2bit conversion
        dna_to_2bit_bidirectional_64(seq, chr_info->chr_len, window_size, v);
        delete seq;
    }
    std::cout << chr_num << " chrs loaded" << std::endl;
    
    // sort
    std::cout << "Sorting vector... (k-mers loaded = " << v.size() << ")" << std::endl;
    std::sort(v.begin(), v.end());
    std::cout << "Sorting finished" << std::endl;
    
    // save redundant k-mers in output_prefix.mk
    std::ofstream outfile;
    outfile.open(out_mk, std::ios::binary);
    if (! outfile.is_open()) {
        return 1;
    }
    
    ull num_saved=0;
    ull veclen=v.size();
    ull n=1;
    size_t type_size=sizeof(uint64_t);
    for (ull i=1; i < veclen; i++) {
        if (! (v[i - 1] == v[i])) {
            if (n >= 2) {
                outfile.write((char*)&v[i - 1], type_size);
                num_saved++;
            }
            n=0;
        }
        n++;
    }
    if (n >= 2) {
        outfile.write((char*)&v[veclen - 1], type_size);
        num_saved++;
    }
    std::cout << num_saved << " redundant k-mers found" << std::endl;
    
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
 This is a wrapper of find_and_save_red_kmers().
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
    
    std::string str_fai=std::string(fa);
    str_fai += ".fai";
    char* f_fai=(char *)str_fai.c_str();
    
    std::string out_prefix=std::string(argv[2]);
    std::string str_mk= out_prefix + ".mk";
    std::string str_mi= out_prefix + ".mi";
    char* out_mk=(char *)str_mk.c_str();
    char* out_mi=(char *)str_mi.c_str();
    
    // file check
    assert(file_checker(fa) == 0);
    assert(file_checker(f_fai) == 0);
    std::ofstream outfile;
    outfile.open(out_mk, std::ios::binary);
    if (! outfile.is_open()) { return 1; }
    outfile.close();
    
    std::cout << "in.fa: " << fa << std::endl;
    std::cout << "in.fa.fai: " << f_fai << std::endl;
    std::cout << "out.mk: " << out_mk << std::endl;
    std::cout << "out.mi: " << out_mi << std::endl;
    
    // process find_and_save_red_kmers()
    int ret=find_and_save_red_kmers(fa, f_fai, out_mk, out_mi);
    
    return ret;
}

