/*
 Author: Shohei Kojima @ RIKEN
 Description:
    This reads a BAM or CRAM file and extracts unmapped reads.
    This only outputs unmapped read when a read contains repeat k-mer(s).
    This also exports the number of unmapped reads.
 Compile:
    g++ -o extract_unmapped -I /path/to/htslib -L /path/to/htslib extract_unmapped.cpp -pthread -O2 -lhts
    g++ -shared -fPIC -o extract_unmapped.so -I /path/to/htslib -L /path/to/htslib extract_unmapped.cpp -pthread -O2
 Usage:
    usage: %prog input.bam/cram input.mk output_dir n_thread [reference.fa]
    (When CRAM file, it requires the reference fasta file.)
 Input:
    1) input.bam/cram : this must have an index file (.bai or .crai).
    2) input.mk : Repeat k-mer file. MEGAnE automatically generates this file.
    3) output_dir : This dir must be present.
    4) n_thread : 1 or more.
    5) (optional) reference.fa : must be specified when CRAM input.
 Output:
    1) unmapped.fa : repeat-derived unmapped reads
 Misc info:
    This requires ~1 GB RAM (for io buffering) in total.
    This depends on external C++ library, htslib. This program is compatible with at least htslib-1.14.
    This program is the complete re-write of `extract_discordant.py` which was used until MEGAnE v0.1.1.
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <unistd.h>
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "dna_to_2bit.hpp"
#include "extract_unmapped.hpp"
using namespace dna_to_2bit_hpp;

#define MAX_READ_LEN 512
#define REP_KMER_SIZE 11
#define SHIFT_16_TO_11 10

typedef unsigned long ul;
typedef unsigned long long ull;
const char* SAM_UNMAPPED_RNAME = "*";


/*
 This is a custom binary search similar to std::lower_bound
 */
inline ull custom_binary_search(const std::vector<uint32_t>& vec, const ull& vec_size, uint32_t& key) {
    if (vec[vec_size-1] < key) {
        return vec_size;
    }
    ull step= 1 << (63 - __builtin_clz(vec_size-1));
    ull pos= vec[step - 1] < key ? vec_size - step - 1 : -1;
    while ((step >>= 1) > 0) {
        pos= (vec[pos + step] < key ? pos + step : pos);
    }
    return pos + 1;
}


/*
 This converts DNA to 11-nt 2bit and judges whether the 11-mer is repeat-derived.
 Args:
    1) seq
    2) length of the seq
    3) vector containing repeat k-mer set
    4) vector size
 */
inline bool is_rep(std::string seq, uint32_t clip_len, const std::vector<uint32_t>& crepkmer, const ull& num_kmer) {
    ull pos=0;
    // first window_size (window_size = REP_KMER_SIZE)
    uint32_t bit2f=0;
    int nn=0;
    for (int i=0; i < REP_KMER_SIZE; i++) {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_32[seq[i]];
        if (seq[i] == 'N' || seq[i] == 'n') {  // ignore when N or n appears
            nn= REP_KMER_SIZE - 1;
        } else if (nn > 0) {  // within window_size-nt from N or n
            nn--;
        }
    }
    pos=custom_binary_search(crepkmer, num_kmer, bit2f);
    if (crepkmer[pos] == bit2f) { return true; }
    
    // rolling calc.
    for (int i=REP_KMER_SIZE; i < clip_len; i++) {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_32[seq[i]];
        bit2f <<= SHIFT_16_TO_11;
        bit2f >>= SHIFT_16_TO_11;
        if (seq[i] == 'N' || seq[i] == 'n') {  // ignore when N or n appears
            nn= REP_KMER_SIZE - 1;
        } else if (nn > 0) {  // within window_size-nt from N or n
            nn--;
        } else {
            pos=custom_binary_search(crepkmer, num_kmer, bit2f);
            if (crepkmer[pos] == bit2f) { return true; }
        }
    }
    return false;
}


/*
 This is a de facto main function.
 */
int extract_unmapped(std::string bam, std::string mk, std::string mi, std::string ref_fa,
                     std::string outdir, int n_thread, bool is_cram) {
    // prep to read file
    std::ifstream infile;
    std::string line;
    
    // load repeat .mk
    infile.open(mi.c_str());
    if (! infile.is_open()) { return 1; }
    if (! std::getline(infile, line)) { return 1; }
    ull num_kmer=std::stoull(line);
    infile.close();
    std::cout << "Number of k-mers loading from " << mk << ": " << num_kmer << std::endl;
    // load .mk
    infile.open(mk.c_str(), std::ios::binary);
    if (! infile.is_open()) { return 1; }
    std::vector<uint32_t> repkmer;
    repkmer.resize(num_kmer);
    infile.read((char*)&repkmer[0], sizeof(uint32_t) * num_kmer);
    infile.close();
    const std::vector<uint32_t>& crepkmer=repkmer;
    // for 2bit conversion
    const int window_size=init_dna_to_2bit_32();
    
    // open bam
    const char *f= bam.c_str();
    htsFile *fp=hts_open(f, "r");
    if (is_cram) {
        hts_set_opt(fp, CRAM_OPT_REFERENCE, ref_fa.c_str());
    }
    sam_hdr_t *h=sam_hdr_read(fp);
    bam1_t *b= bam_init1();
    
    // open bai
    hts_idx_t *idx = nullptr;
    idx=sam_index_load(fp, f);
    
    // threading
    htsThreadPool p = {nullptr, 0};
    p.pool = hts_tpool_init(n_thread);
    hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
    
    // ofstreams (C fopen)
    std::FILE* ofs_unmapped = fopen((outdir + std::string("/unmapped.fa")).c_str(), "w");
    if (ofs_unmapped == nullptr) { return 1; }
    
    // extract unmapped
    hts_itr_t *iter = sam_itr_querys(idx, h, SAM_UNMAPPED_RNAME);
    char* seq= new char[MAX_READ_LEN];    // seq of read, ATGCN
    std::string seqstr;
    int ret;
    while ((ret = sam_itr_next(fp, iter, b)) >= 0) {    // end file = -1
        // get read seq
        char* qname=bam_get_qname(b);      // read name
        int32_t &l_qseq = b->core.l_qseq;  // length of read
        uint8_t *tmp_s  = bam_get_seq(b);  // seq of read, nt16
        uint16_t &flag  = b->core.flag;    // SAM flag
        bool is_read2   = (flag & BAM_FREAD2) > 0;
        for (int i=0; i < l_qseq; i++) {
            seq[i]=seq_nt16_str[bam_seqi(tmp_s, i)];  // get nucleotide id and convert into IUPAC id
        }
        seq[l_qseq]='\0';
        seqstr=seq;
        // check rep k-mer
        if (is_rep(seqstr, l_qseq, crepkmer, num_kmer)) {
            // output
            std::fprintf(ofs_unmapped, "%s/%d\n%s\n", qname, ((int)is_read2 + 1), seq);
        }
    }
    
    // close file obj
    std::fclose(ofs_unmapped);
    
    // close bam
    sam_hdr_destroy(h);
    sam_close(fp);
    bam_destroy1(b);
    
    return 0;
}


/*
 This is a main func for direct use.
 usage: %prog input.bam/cram input.mk output_dir n_thread [reference.fa]
 */
int main(int argc, char *argv[]) {
    // file check (not decent)
    bool is_cram=false;
    if (argc == 5) {
        // bam
    } else if (argc == 6) {
        // cram
        is_cram=true;
    } else {
        std::cerr << "Please specify required files." << std::endl;
        return 1;
    }
    
    std::string bam=argv[1];
    std::string mk=argv[2];
    std::string mi=argv[2];
    mi.pop_back();
    mi += 'i';
    std::string outdir=argv[3];
    int n_thread= atoi(argv[4]);
    std::string ref_fa="";
    if (argc == 6) {
        ref_fa=argv[5];
    }
    
    std::cout << "bam " << bam << std::endl;
    std::cout << "mk " << mk << std::endl;
    std::cout << "mi " << mi << std::endl;
    std::cout << "outdir " << outdir << std::endl;
    std::cout << "n_thread " << n_thread << std::endl;
    std::cout << "ref_fa " << ref_fa << std::endl;
    std::cout << "is_cram " << is_cram << std::endl;
    
    int ret = extract_unmapped(bam, mk, mi, ref_fa, outdir, n_thread, is_cram);
    return ret;
}
