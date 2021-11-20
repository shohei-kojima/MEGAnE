/*
 Author: Shohei Kojima @ RIKEN
 This is a minimal DNA to 2bit converter.
 */

#ifndef DNA_TO_2BIT_HPP
#define DNA_TO_2BIT_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>
typedef unsigned long long ull;
typedef unsigned long ul;


extern "C" {
namespace dna_to_2bit_hpp {


/*
 uint32
 */
static uint32_t dna_to_2bitf_32[128];
static uint32_t dna_to_2bitr_32[128];

const int init_dna_to_2bit_32() {
    int shift=30;
    const int window_size=16;
    
    dna_to_2bitf_32['A']=0ul; dna_to_2bitf_32['a']=0ul;
    dna_to_2bitf_32['T']=1ul; dna_to_2bitf_32['t']=1ul;
    dna_to_2bitf_32['G']=2ul; dna_to_2bitf_32['g']=2ul;
    dna_to_2bitf_32['C']=3ul; dna_to_2bitf_32['c']=3ul;
    dna_to_2bitf_32['N']=0ul; dna_to_2bitf_32['n']=0ul;
    
    dna_to_2bitr_32['A']=1ul << shift; dna_to_2bitr_32['a']=1ul << shift;
    dna_to_2bitr_32['T']=0ul << shift; dna_to_2bitr_32['t']=0ul << shift;
    dna_to_2bitr_32['G']=3ul << shift; dna_to_2bitr_32['g']=3ul << shift;
    dna_to_2bitr_32['C']=2ul << shift; dna_to_2bitr_32['c']=2ul << shift;
    dna_to_2bitr_32['N']=0ul << shift; dna_to_2bitr_32['n']=0ul << shift;
    
    return window_size;
}

/*
 This converts 16-nt DNA to 2bit and stores as uint32_t in the vector.
 Args:
    1) pointer to seq
    2) seq len (ull)
    3) window_size (must be 16)
    4) vector to store the 2bit seq
 */
void dna_to_2bit_bidirectional_32(char* seq, ull& seqlen, const int& window_size, std::vector<uint32_t>& v) {
    if (seqlen < window_size) {
        return;
    }
    
    // first window_size
    uint32_t bit2f=0;
    uint32_t bit2r=0;
    int nn=0;
    for (int i=0; i < window_size; i++) {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_32[seq[i]];
        bit2r >>= 2;
        bit2r |= dna_to_2bitr_32[seq[i]];
        if (seq[i] == 'N' || seq[i] == 'n') {  // ignore when N or n appears
            nn= window_size - 1;
        } else if (nn > 0) {  // within window_size-nt from N or n
            nn--;
        }
    }
    if (nn == 0) {
        if (bit2f == bit2r) {
            v.push_back(bit2f);
        } else {
            v.push_back(bit2f);
            v.push_back(bit2r);
        }
    }
    
    // rolling calc.
    for (ull i=window_size; i < seqlen; i++) {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_32[seq[i]];
        bit2r >>= 2;
        bit2r |= dna_to_2bitr_32[seq[i]];
        if (seq[i] == 'N' || seq[i] == 'n') {  // ignore when N or n appears
            nn= window_size - 1;
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
 uint64
 */
static uint64_t dna_to_2bitf_64[128];
static uint64_t dna_to_2bitr_64[128];

const int init_dna_to_2bit_64() {
    int shift=62;
    const int window_size=32;
    
    dna_to_2bitf_64['A']=0ull; dna_to_2bitf_64['a']=0ull;
    dna_to_2bitf_64['T']=1ull; dna_to_2bitf_64['t']=1ull;
    dna_to_2bitf_64['G']=2ull; dna_to_2bitf_64['g']=2ull;
    dna_to_2bitf_64['C']=3ull; dna_to_2bitf_64['c']=3ull;
    dna_to_2bitf_64['N']=0ull; dna_to_2bitf_64['n']=0ull;
    
    dna_to_2bitr_64['A']=1ull << shift; dna_to_2bitr_64['a']=1ull << shift;
    dna_to_2bitr_64['T']=0ull << shift; dna_to_2bitr_64['t']=0ull << shift;
    dna_to_2bitr_64['G']=3ull << shift; dna_to_2bitr_64['g']=3ull << shift;
    dna_to_2bitr_64['C']=2ull << shift; dna_to_2bitr_64['c']=2ull << shift;
    dna_to_2bitr_64['N']=0ull << shift; dna_to_2bitr_64['n']=0ull << shift;
    
    return window_size;
}

/*
 This converts 32-nt DNA to 2bit and stores as uint64_t in the vector.
 Args:
    1) pointer to seq
    2) seq len (ull)
    3) window_size (must be 32)
    4) vector to store the 2bit seq
 */
void dna_to_2bit_bidirectional_64(char* seq, ull& seqlen, const int& window_size, std::vector<uint64_t>& v) {
    if (seqlen < window_size) {
        return;
    }
    
    // first window_size
    uint64_t bit2f=0;
    uint64_t bit2r=0;
    int nn=0;
    for (int i=0; i < window_size; i++) {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_64[seq[i]];
        bit2r >>= 2;
        bit2r |= dna_to_2bitr_64[seq[i]];
        if (seq[i] == 'N' || seq[i] == 'n') {
            nn= window_size - 1;
        } else if (nn > 0) {  // within window_size-nt from N or n
            nn--;
        }
    }
    if (nn == 0) {
        if (bit2f == bit2r) {
            v.push_back(bit2f);
        } else {
            v.push_back(bit2f);
            v.push_back(bit2r);
        }
    }
    
    // rolling calc.
    for (ull i=window_size; i < seqlen; i++) {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_64[seq[i]];
        bit2r >>= 2;
        bit2r |= dna_to_2bitr_64[seq[i]];
        if (seq[i] == 'N' || seq[i] == 'n') {
            nn= window_size - 1;
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
 This is the main for testing - this is also an example how to use this.
 g++ -o dna_to_2bit dna_to_2bit.cpp
 */
//int _main() {
//    const char* seq ="ATGCATCGACTAGCATCGACTAGCATGACnAGATGCATCGACTAGCATCGACTAGCATGACTAC";
//    const char* cseq="GTAGTCATGCTAGTCGATGCTAGTCGATGCATCTnGTCATGCTAGTCGATGCTAGTCGATGCAT";
//    const char* seq ="CGACTAGCATGACTAC";
//    const char* cseq="GTAGTCATGCTAGTCG";
//    ull seqlen=64;
//    const int window_size=init_dna_to_2bit_64();
//
//    std::vector<uint64_t> v;
//    dna_to_2bit_bidirectional_64(seq,  seqlen, window_size, v);
//    dna_to_2bit_bidirectional_64(cseq, seqlen, window_size, v);
//
//    std::sort(v.begin(), v.end());
//    for (uint64_t b : v) {
//        std::cout << b << std::endl;
//    }
//
//    return 0;
//}



}  // namespace
}  // extern C
#endif
