/*
 Author: Shohei Kojima @ RIKEN
 This reads a relatively small fasta file (e.g., total 1MB seq)
 and hash the sequences for both forward and reverse strands.
 */


#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>
#include "complementary_seq.hpp"
#include "hash_rep.hpp"
using namespace complementary_seq_hpp;

typedef unsigned long long ull;
const ull HASH_BASE = 100000007;
const char FASTA_HEADER_START = '>';


/*
 Hash a single seq
 Args:
    1) pointer to sequence
    2) length of sequence
    3) t (HASH_BASE of seed_len)
    4) vector that will hold hash values
 */
inline ull hash_seq(const char* seq, int& seed_len, ull& t, std::vector<int>& hsv) {
    // calc. hash value for first SEED_LEN
    ull hv=0;  // hash value
    char prev[seed_len];
    int prev_pos=0;
    for (int i=0; i < seed_len; i++) {
        hv = hv * HASH_BASE + seq[i];
        prev[prev_pos]=seq[i];
        prev_pos++;
    }
    hsv.push_back(hv);
    
    // calc. hash val throughout seq
    ull n_char_read=seed_len;
    prev_pos=0;
    while (seq[n_char_read]) {  // calc. hash value for every SEED_LEN window
        hv = (hv * HASH_BASE) + seq[n_char_read] - (prev[prev_pos] * t);
        hsv.push_back(hv);
        prev[prev_pos]=seq[n_char_read];
        prev_pos = ++prev_pos % seed_len;
        n_char_read++;
    }
    
    return n_char_read;  // return length of input seq
}


/*
 Hash a single seq for both forward and reverse strand.
 Wrapper of hash_seq().
 Args:
    1) pointer to sequence
    2) length of sequence
    3) t (HASH_BASE of seed_len)
    4) vector that will hold hash values
 */
void rolling_hash_bidirectional(const char* seq, int& seed_len, ull& t, std::vector<int>& hsv) {
    ull seq_len;
    
    // hash forward strand
    seq_len=hash_seq(seq, seed_len, t, hsv);
    
    // hash reverse strand
    const char* cseq=(const char*)complementary_seq(seq, seq_len);
    seq_len=hash_seq(seq, seed_len, t, hsv);
    
    delete[] cseq;
}


/*
 Reads a fasta file with multiple seqs and calc. hash values
 Args:
    1) name of fasta file
    2) length to calc. hash value (seed_len)
    3) vector that will hold hash values
 */
int hash_a_fasta_file(char* fa, int seed_len, std::vector<int>& hvs) {
    // to hash
    comp_init();
    ull t=1;  // HASH_BASE of seed_len
    for (int i=0; i < seed_len; i++) {  // calc. hash value for first SEED_LEN
        t  = t  * HASH_BASE;
    }
    
    // open fasta
    std::ifstream infile;
    infile.open(fa);
    if (! infile.is_open()) {
        return 1;
    }
    
    // read fasta
    std::string line;
    std::string seq;
    while (std::getline(infile, line)) {
        if (line[0] == FASTA_HEADER_START) {
            if (! seq.empty()) {
                if (seq.size() >= seed_len) {
                    rolling_hash_bidirectional(seq.c_str(), seed_len, t, hvs); // hash
                }
                seq.clear();  // clear for next fasta entry
            }
        } else {
            transform(line.begin(), line.end(), line.begin(),
                [](unsigned char c) {return toupper(c);}); // convert to upper letters
            seq += line;
        }
    }
    if (seq.size() >= seed_len) {
        rolling_hash_bidirectional(seq.c_str(), seed_len, t, hvs); // hash
    }
    
    return 0;
}


/*
 This is the de facto main function.
 Reads a fasta file with multiple seqs and calc. hash values.
 Args:
    1) name of fasta file
    2) length to calc. hash value (seed_len)
 */
int hash_rep(char* fa, int seed_len) {
    // hash
    std::vector<int> hvs;  // stores hash values
    int ret=hash_a_fasta_file(fa, seed_len, hvs);
    if (! ret == 0) {
        return ret;
    }
    std::cout << hvs.size() << std::endl;
    
    // sort vec
    std::sort(hvs.begin(), hvs.end());
    hvs.erase(std::unique(hvs.begin(), hvs.end()), hvs.end());
    
    return 0;
}


/*
 This is just a wrapper of hash_rep().
 Usage: ./prog input.fa seed_len
 */
int _main(int argc, char* argv[]) {
    if (argc <= 2) {
        std::cerr << "Please specify a fasta file and seq length to be hashed." << std::endl;
        return 1;
    }
    char* fa=argv[1];
    int seed_len=atoi(argv[2]);
    std::cout << "Reading " << argv[1] << " ..." << std::endl;

    int ret=hash_rep(fa, seed_len);
    return ret;
}

