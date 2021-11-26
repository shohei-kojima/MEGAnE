/*
 Author: Shohei Kojima @ RIKEN
 Description:
    This reads a BAM or CRAM file and extracts discordantly mapped reads
    based on cigar, SA tag, and XA tag.
    Discordantly mapped reads are:
        1) chimeric reads derived from ME insertions
        2) chimeric reads derived from absent MEs
        3) hybrid reads derived from ME insertions
    This also exports the number of discordantly mapped reads.
 Compile:
    g++ -o extract_discordant -I /path/to/htslib/htslib-1.13 -L /path/to/htslib/htslib-1.13 extract_discordant.cpp -lhts -pthread -O2
    g++ -shared -fPIC -o extract_discordant.so -I /path/to/htslib/htslib-1.13 -L /path/to/htslib/htslib-1.13 extract_discordant.cpp -lhts -pthread -O2
 Usage:
    usage: %prog input.bam/cram main_chrs.txt input.mk output_dir n_thread [reference.fa]
    (When CRAM file, it requires the reference fasta file.)
 Input:
    1) input.bam/cram : this must have an index file (.bai or .crai).
    2) main_chrs.txt : file containing the names of main chrs.
    3) input.mk : Repeat k-mer file. MEGAnE automatically generates this file.
    4) output_dir : This dir must be present.
    5) n_thread : 1 or more.
    6) (optional) reference.fa : must be specified when CRAM input.
 Output:
    1) overhang_pA.txt[thread_id] : pA reads
    2) overhang.fa[thread_id] : all chimeric reads
    3) mapped.fa[thread_id] : mapped regions of pA and chimeric reads
    4) distant.txt[thread_id] : hybrid reads
    5) absent.txt[thread_id] : chimeric reads coming from absence
    6) stats.txt[thread_id] : discordant read stats
 Misc info:
    This requires ~2 GB RAM in total.
    This depends on external C++ library, htslib. This program is compatible with at least htslib-1.13.
    This program is the complete re-write of `extract_discordant.py` which was used until MEGAnE v0.1.1.
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <functional>
#include <unistd.h>
#include "htslib/sam.h"
#include "dna_to_2bit.hpp"
#include "complementary_seq.hpp"
#include "extract_discordant.hpp"
#include "ThreadPool.h"
using namespace dna_to_2bit_hpp;
using namespace complementary_seq_hpp;

#define MAX_CIGAR_LEN 128
#define TMP_BUF_SIZE  131072
#define READ_PAIR_GAP_LEN 2000
#define MAX_SEQ_LEN 512
#define DISCORDANT_READS_CLIP_LEN 20
#define REP_KMER_SIZE 11
#define SHIFT_16_TO_11 10
#define ABS_MIN_DIST 50
#define ABS_MAX_DIST 20000
#define MAPPED_REGION_LOW_COMPLEX_THRESHOLD 0.7
#define POLYA_OVERHANG_THRESHOLD 0.7

typedef unsigned long ul;
typedef unsigned long long ull;

const char ATGC[]="ATGC";


/*
 to do: implement cram reader; CRAM_OPT_REFERENCE;
 if (hts_set_opt(out, CRAM_OPT_REFERENCE, outref) < 0) {
     fail("setting reference %s for %s", outref, outfname);
     goto err;
 }
https://github.com/samtools/htslib/blob/9672589346459d675d62851d5b7b5f2e5c919076/test/sam.c
 */


/*
 Input files and options.
 */
struct options {
    std::string bam;
    std::string f_mainchr;
    std::string mk;
    std::string mi;
    std::string outdir;
    std::string ref_fa;
    int n_thread;
    bool is_cram;
    
    options(std::string bam, std::string f_mainchr, std::string mk, std::string mi, std::string ref_fa,
            std::string outdir, int n_thread, bool is_cram) {
        this->bam           = bam;
        this->f_mainchr     = f_mainchr;
        this->mk            = mk;
        this->mi            = mi;
        this->outdir        = outdir;
        this->ref_fa        = ref_fa;
        this->n_thread      = n_thread;
        this->is_cram       = is_cram;
    }
};


/*
 Struct to keep track of read stats.
 */
struct read_stats {
    int64_t pA=0;
    int64_t chimeric=0;
    int64_t distant=0;
    int64_t absent=0;
};


/*
 Struct to keep file streams during threading.
 */
struct fstr {
    std::FILE* ofs_pA;
    std::FILE* ofs_overhang;
    std::FILE* ofs_distant;
    std::FILE* ofs_mapped;
    std::FILE* ofs_abs;
    std::FILE* ofs_stat;
    bool is_occupied;
    
    fstr(std::FILE* ofs_pA, std::FILE* ofs_overhang, std::FILE* ofs_distant, std::FILE* ofs_mapped,
         std::FILE* ofs_abs, std::FILE* ofs_stat, bool is_occupied) {
        this->ofs_pA=ofs_pA;
        this->ofs_overhang=ofs_overhang;
        this->ofs_distant=ofs_distant;
        this->ofs_mapped=ofs_mapped;
        this->ofs_abs=ofs_abs;
        this->ofs_stat=ofs_stat;
        this->is_occupied=is_occupied;
    }
};


/*
 Manage file streams during threading.
 This will be accessed from multiple threads - do mutex when occupy() and release().
 */
class cfstrs {
    std::mutex _mutex;
    std::vector<fstr*> vec;  // stores fstr objects
public:
    void push_back_new(std::FILE* ofs_pA, std::FILE* ofs_overhang, std::FILE* ofs_distant, std::FILE* ofs_mapped,
                       std::FILE* ofs_abs, std::FILE* ofs_stat) {
        fstr* fstrp = new fstr(ofs_pA, ofs_overhang, ofs_distant, ofs_mapped, ofs_abs, ofs_stat, false);
        this->vec.push_back(fstrp);
    }
    
    fstr* occupy() {
        std::unique_lock<std::mutex> lock(this->_mutex);  // need mutex
        fstr* fstrp=nullptr;
        for (fstr* t : this->vec) {
            if (t->is_occupied == false) {  // if not the ofstream is used by the threads
                t->is_occupied=true;  // take this ofstream by a thread
                fstrp=t;
                break;
            }
        }
        if (fstrp == nullptr) {  // this happens if the num of ofstream object was fewer than the thread num.
            std::cerr << "ERR: no available ofstream found." << std::endl;
            std::exit(1);
        }
        return fstrp;
    }
    
    void release(fstr* fstrp) {
        std::unique_lock<std::mutex> lock(this->_mutex);  // need mutex
        fstrp->is_occupied=false;  // release this ofstream from a thread
    }
    
    void close_fileobjs() {
        for (fstr* t : this->vec) {
            fclose(t->ofs_pA);
            fclose(t->ofs_overhang);
            fclose(t->ofs_distant);
            fclose(t->ofs_mapped);
            fclose(t->ofs_abs);
            fclose(t->ofs_stat);
        }
    }
    
    ~cfstrs() {
        for (fstr* t : this->vec) {
            delete t;
        }
    }
};


/*
 Class cfstrs that will be accessed from multiple threads.
 */
static cfstrs fstrs;


/*
 Stores info of clipped reads, SA tags, XA tags.
 */
struct softclip_info {
    std::string chr;
    int64_t     pos;
    bool        is_reverse;
    char        breakpoint;
    int32_t     l_clip_len;
    int32_t     r_clip_len;
    uint64_t    clipstart;
    uint64_t    clipend;
    uint64_t    rlen;
    
    softclip_info(std::string chr, int64_t pos, bool is_reverse, char breakpoint,
                  int32_t l_clip_len, int32_t r_clip_len, uint64_t clipstart, uint64_t clipend, uint64_t rlen) {
        this->chr        = chr;
        this->pos        = pos;
        this->is_reverse = is_reverse;
        this->breakpoint = breakpoint;
        this->l_clip_len = l_clip_len;
        this->r_clip_len = r_clip_len;
        this->clipstart  = clipstart;
        this->clipend    = clipend;
        this->rlen       = rlen;
    };
    
    // for debug, just a print function
    void print() {
        std::cout << this->chr << " " << this->pos << " "
        << this->is_reverse << " " << this->breakpoint << " "
        << this->l_clip_len << " " << this->r_clip_len << " "
        << this->clipstart << " " << this->clipend << " "
        << this->rlen << std::endl;
    }
};


/*
 Stores absent read info.
 */
struct abs_info {
    std::vector<softclip_info*> R;
    std::vector<softclip_info*> L;
    int R_num;
    int L_num;
    
    abs_info() {
        this->R_num=0;
        this->L_num=0;
    }
};


/*
 Manages abs_info.
 */
class abs_info_manager {
public:
    std::unordered_map<std::string, abs_info*> sa_p;
    std::unordered_map<std::string, abs_info*> sa_m;
    std::unordered_map<std::string, abs_info*> xa_p;
    std::unordered_map<std::string, abs_info*> xa_m;
    std::unordered_map<std::string, abs_info*>* ums[4];
    
    void init() {
        this->ums[0]=&(this->sa_p);
        this->ums[1]=&(this->sa_m);
        this->ums[2]=&(this->xa_p);
        this->ums[3]=&(this->xa_m);
    }
    
    void emplace_chrs(const std::vector<std::string>& main_chrs) {
        abs_info* _abs_info;
        for (int i=0; i < 4; i++) {
            for (std::string chr : main_chrs){
                _abs_info = new abs_info();
                (*(this->ums[i])).emplace(chr, _abs_info);
            }
        }
    }
    
    void clean_up_by_chr(const std::string& chr) {
        for (int i=0; i < 4; i++) {
            if (this->ums[i]->at(chr)->R_num > 0) {
                this->ums[i]->at(chr)->R.clear();
                this->ums[i]->at(chr)->R_num = 0;
            }
            if (this->ums[i]->at(chr)->L_num > 0) {
                this->ums[i]->at(chr)->L.clear();
                this->ums[i]->at(chr)->L_num = 0;
            }
        }
    }
    
    ~abs_info_manager() {
        std::unordered_map<std::string, abs_info*>::iterator it;
        for (int i=0; i < 4; i++) {
            it=this->ums[i]->begin();
            while (it != this->ums[i]->end()) {
                delete it->second;
                it++;
            }
        }
    }
};


/*
 Returns chr names (tid) with at least one read.
 Chr names are in vector sorted by number of reads with descending manner.
 This is only available for BAM.
 */
void sort_chr_order(hts_idx_t *idx, sam_hdr_t *h, std::vector<int> &sorted_chr) {
    // retrieve mapped read counts from bam index
    int nseq = hts_idx_nseq(idx);    // number of chrs
    uint64_t mapped,umapped;
    std::vector<std::pair<int64_t, int32_t>> mapped_counts;  // pair<mapped_read_num, tid>
    mapped_counts.reserve(nseq);
    std::cout << nseq << " chrs found." << std::endl;
    for (int32_t tid=0; tid < nseq; tid++) {
        hts_idx_get_stat(idx, tid, &mapped, &umapped);  // return can be -1 when no read on chr
        if (mapped >= 1) {
            mapped_counts.push_back(std::make_pair(mapped, tid));
        }
    }
    
    // ascending by mapped read counts
    std::sort(mapped_counts.begin(), mapped_counts.end());
    int vec_size=mapped_counts.size();
    
    // convert to desceding order
    sorted_chr.reserve(vec_size);
    for (int i= vec_size - 1; i >= 0; i--) {
        sorted_chr.push_back(mapped_counts[i].second);
    }
}


/*
 Returns chr names (tid).
 Chr names are in vector sorted by the lengths of chrs.
 This is available for CRAM.
 */
void sort_chr_order(sam_hdr_t *h, std::vector<int> &sorted_chr) {
    // retrieve mapped read counts from bam index
    int32_t nseq = h->n_targets;    // number of chrs
    std::cout << nseq << " chrs found." << std::endl;
    std::vector<std::pair<int64_t, int32_t>> ref_lens;  // pair<ref_len, tid>
    ref_lens.reserve(nseq);
    for (int32_t tid=0; tid < nseq; tid++) {
        ref_lens.push_back(std::make_pair(h->target_len[tid], tid));
    }
    
    // ascending by mapped read counts
    std::sort(ref_lens.begin(), ref_lens.end());
    int vec_size=ref_lens.size();
    
    // convert to desceding order
    sorted_chr.reserve(vec_size);
    for (int i= vec_size - 1; i >= 0; i--) {
        sorted_chr.push_back(ref_lens[i].second);
    }
}


/*
 Parses cigar and stores in an array. Also checks whether the cigar contains H and S.
 Array type is std::pair<char, uint32_t>, which is <cigar_opchr, bam_cigar_oplen>.
 */
inline void parse_cigar(uint32_t* cigar, uint32_t& n_cigar,
                        std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                        bool& contains_H, bool& contains_S) {
    for (int i=0; i < n_cigar; ++i) {
        char opchr=bam_cigar_opchr(cigar[i]);
        cigar_arr[i]=std::make_pair(opchr, bam_cigar_oplen(cigar[i]));
        if (opchr == 'H') { contains_H=true; }
        if (opchr == 'S') { contains_S=true; }
    }
}


/*
 Parses cigar and stores in an array.
 Array type is std::pair<char, uint32_t>, which is <cigar_opchr, bam_cigar_oplen>.
 */
inline void parse_cigar(uint32_t* cigar, uint32_t& n_cigar,
                        std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN]) {
    for (int i=0; i < n_cigar; ++i) {
        cigar_arr[i]=std::make_pair(bam_cigar_opchr(cigar[i]),
                                    bam_cigar_oplen(cigar[i]));
    }
}


/*
 Defines breakpoint.
 */
inline void define_breakpoint(std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN], uint32_t& n_cigar,
                             char& breakpoint, int32_t& l_clip_len, int32_t& r_clip_len,
                              int64_t& clipstart, int64_t& clipend, int32_t& l_qseq) {
    // left length
    if (cigar_arr[0].first == 'S') {
        l_clip_len=cigar_arr[0].second;
    } else {
        l_clip_len=0;
    }
    // right length
    if (cigar_arr[n_cigar - 1].first == 'S') {
        r_clip_len=cigar_arr[n_cigar - 1].second;
    } else {
        r_clip_len=0;
    }
    // judge
    if (l_clip_len == r_clip_len) {
        breakpoint ='N';  // undetermined
        clipstart  = 0;
        clipend    = 0;
    }
    if (l_clip_len > r_clip_len) {
        breakpoint ='L';
        clipstart  = 0;
        clipend    = l_clip_len;
    } else {
        breakpoint ='R';
        clipstart  = l_qseq - r_clip_len;
        clipend    = l_qseq;
    }
}


/*
 Parses one SA entry in an SA tag.
 */
inline void parse_one_SA(std::string& tmp_str,
                         std::vector<softclip_info>& softclips,
                         std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                         uint32_t*& cig_buf, size_t& cig_m,
                         char& breakpoint, int32_t& l_clip_len, int32_t& r_clip_len,
                         int64_t& clipstart, int64_t& clipend, int32_t& l_qseq) {
    // split SA info by ','
    int comma1=0,comma2=0,comma3=0,comma4=0;
    int n=0, pos=0;
    for (char c : tmp_str) {
        if (c == ',') {
            if (n == 0) {
                comma1=pos;
            } else if (n == 1) {
                comma2=pos;
            } else if (n == 2) {
                comma3=pos;
            } else {
                comma4=pos;
                break;
            }
            n++;
        }
        pos++;
    }
    
    // parse cigar
    cig_m=0;
    uint32_t _n_cigar = sam_parse_cigar((tmp_str.substr(comma3 + 1, comma4 - comma3 - 1)).c_str(),
                                        nullptr, &cig_buf, &cig_m);
    int64_t _rlen     = bam_cigar2rlen(_n_cigar, cig_buf);
    parse_cigar(cig_buf, _n_cigar, cigar_arr);
    
    // save breakpoint info
    define_breakpoint(cigar_arr, _n_cigar, breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq);
    std::string _chr=tmp_str.substr(0, comma1);
    int64_t _pos=std::stoll(tmp_str.substr(comma1 + 1, comma2 - comma1 - 1)) - 1;  // 0-based
    char _strand=tmp_str[comma2 + 1];
    bool _is_reverse;
    if (_strand == '+') {
        _is_reverse=false;
    } else {
        _is_reverse=true;
    }
    softclips.push_back(softclip_info(_chr, _pos, _is_reverse, breakpoint,
                                      l_clip_len, r_clip_len, clipstart, clipend, _rlen));
}


/*
 Parses char SA tag and stores in std::vector<softclip_info>
 */
inline void SA_parser(char* sa_tag, std::string& tmp_str,
                      std::vector<softclip_info>& softclips,
                      std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                      uint32_t*& cig_buf, size_t& cig_m,
                      char& breakpoint, int32_t& l_clip_len, int32_t& r_clip_len,
                      int64_t& clipstart, int64_t& clipend, int32_t& l_qseq) {
//    std::cout << sa_tag << std::endl;
    // count SA entires
    int n_entry=0;
    int sa_len=0;
    tmp_str.clear();
    while (sa_tag[sa_len]) {
        if (sa_tag[sa_len] == ';') {
            parse_one_SA(tmp_str, softclips, cigar_arr, cig_buf, cig_m,
                         breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq);
            tmp_str.clear();
        } else {
            tmp_str += sa_tag[sa_len];
        }
        sa_len++;
    }
}


/*
 Parses one XA entry in an XA tag.
 */
inline void parse_one_XA(std::string& tmp_str,
                         std::vector<softclip_info>& softclips,
                         std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                         uint32_t*& cig_buf, size_t& cig_m,
                         char& breakpoint, int32_t& l_clip_len, int32_t& r_clip_len,
                         int64_t& clipstart, int64_t& clipend, int32_t& l_qseq) {
    // split XA info by ','
    int comma1=0,comma2=0,comma3=0;
    int n=0, pos=0;
    for (char c : tmp_str) {
        if (c == ',') {
            if (n == 0) {
                comma1=pos;
            } else if (n == 1) {
                comma2=pos;
            } else {
                comma3=pos;
                break;
            }
            n++;
        }
        pos++;
    }
    
    // parse cigar
    cig_m=0;
    uint32_t _n_cigar = sam_parse_cigar((tmp_str.substr(comma2 + 1, comma3 - comma2 - 1)).c_str(),
                                        nullptr, &cig_buf, &cig_m);
    int64_t _rlen     = bam_cigar2rlen(_n_cigar, cig_buf);
    parse_cigar(cig_buf, _n_cigar, cigar_arr);
    
    // save breakpoint info
    define_breakpoint(cigar_arr, _n_cigar, breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq);
    std::string _chr=tmp_str.substr(0, comma1);
    int64_t _pos=std::stoll(tmp_str.substr(comma1 + 2, comma2 - comma1 - 2)) - 1;  // 0-based
    char _strand=tmp_str[comma1 + 1];
    bool _is_reverse;
    if (_strand == '+') {
        _is_reverse=false;
    } else {
        _is_reverse=true;
    }
    softclips.push_back(softclip_info(_chr, _pos, _is_reverse, breakpoint,
                                      l_clip_len, r_clip_len, clipstart, clipend, _rlen));
}


/*
 Parses char XA tag and stores in std::vector<softclip_info>
 */
inline void XA_parser(char* sa_tag, std::string& tmp_str,
                      std::vector<softclip_info>& softclips,
                      std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                      uint32_t*& cig_buf, size_t& cig_m,
                      char& breakpoint, int32_t& l_clip_len, int32_t& r_clip_len,
                      int64_t& clipstart, int64_t& clipend, int32_t& l_qseq) {
    // count SA entires
    int n_entry=0;
    int sa_len=0;
    tmp_str.clear();
    while (sa_tag[sa_len]) {
        if (sa_tag[sa_len] == ';') {
            parse_one_XA(tmp_str, softclips, cigar_arr, cig_buf, cig_m,
                         breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq);
            tmp_str.clear();
        } else {
            tmp_str += sa_tag[sa_len];
        }
        sa_len++;
    }
}


/*
 Simple repeat checker.
 */
inline bool is_simple_repeat(const std::string& seq, int32_t seqlen) {
    int32_t count=0;
    for (char c: ATGC) {
        for (int32_t i=0; i < seqlen; i++) {
            if (seq[i] == c) { count++; }
        }
        if (((double) count / seqlen) >= MAPPED_REGION_LOW_COMPLEX_THRESHOLD) {
            return true;
        }
        count=0;
    }
    return false;
}


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
inline bool is_rep_overhang(std::string seq, uint64_t clip_len, const std::vector<uint32_t>& crepkmer, const ull& num_kmer) {
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
 This 1) output pA and 2) stores overhang info and mapped seq to unordered_maps.
 */
inline bool output_pA(std::vector<softclip_info>& softclips, const std::vector<uint32_t>& crepkmer, const ull& num_kmer,
                     char* chr, std::string& fseq, std::string& rseq, int64_t& start, int64_t& end,
                     char* qname, int32_t& l_qseq, bool& is_read2, bool& is_reverse, char& strand,
                     std::string& tmp_str, std::string& tmp_str1, std::string& tmp_str2, std::string& tmp_str3,
                     std::unordered_map<std::string, std::string>& um1,
                     std::unordered_map<std::string, std::string>& um2,
                     char* tmp_buf, fstr* fs, const std::unordered_map<std::string, bool>& is_mainchr) {
    bool written=false;
    for (softclip_info s : softclips) {
        if (s.breakpoint == 'N') { continue; }
        if (! is_mainchr.at(s.chr)) { continue; }
        if ((s.l_clip_len >= DISCORDANT_READS_CLIP_LEN) || (s.r_clip_len >= DISCORDANT_READS_CLIP_LEN)) {
            tmp_str.clear();  // seq of read
            if (s.is_reverse && is_reverse) {
                tmp_str=fseq;
            } else if ((! s.is_reverse) && (! is_reverse)) {
                tmp_str=fseq;
            } else {
                tmp_str=rseq;
            }
            int32_t mapped_len= l_qseq - s.l_clip_len - s.r_clip_len;
            tmp_str1.clear();  // mapped seq
            tmp_str1=tmp_str.substr(s.l_clip_len, mapped_len);
            if (is_simple_repeat(tmp_str1, mapped_len)) { continue; }
            tmp_str2.clear();  // clipped seq
            int32_t Acount=0;
            if (s.breakpoint == 'L') {
                tmp_str2=tmp_str.substr(s.clipstart, s.clipend);
                for (char c: tmp_str2) {
                    if (c == 'A') { Acount++; }
                }
            } else {
                tmp_str2=tmp_str.substr(s.clipstart);
                for (char c: tmp_str2) {
                    if (c == 'T') { Acount++; }
                }
            }
            // judge whether either pA read or ME overhang
            bool pA=false;
            bool ME_overhang=false;
            if (((double) Acount / (s.clipend - s.clipstart)) > POLYA_OVERHANG_THRESHOLD) {
                pA=true;
            } else if (is_rep_overhang(tmp_str2, (s.clipend - s.clipstart), crepkmer, num_kmer)) {  // proceed if repeat overhang
                ME_overhang=true;
            }
            if (pA || ME_overhang) {
                std::sprintf(tmp_buf, "%s:%ld-%ld/%c/%s/%d/%c", chr, s.pos, (s.pos + s.rlen), s.breakpoint, qname, (is_read2 + 1), strand);
                tmp_str3.clear();  // header name
                tmp_str3=tmp_buf;
                if (pA) {
                    std::fprintf(fs->ofs_pA, "%s\t%ld\n", tmp_buf, (s.clipend - s.clipstart));
                    written=true;
                } else {
                    auto itr=um2.find(tmp_str2);  // um2 = (std::string, std::string) = (clipped seq, header)
                    if (itr == um2.end()) {
                        um2.emplace(tmp_str2, tmp_str3);
                    } else {
                        itr->second = itr->second + tmp_str3;
                    }
                }
                auto itr=um1.find(tmp_str1);  // um1 = (std::string, std::string) = (mapped seq, header)
                if (itr == um1.end()) {
                    um1.emplace(tmp_str1, tmp_str3);
                } else {
                    itr->second = itr->second + tmp_str3;
                }
            }
        }
    }
    return written;
}


/*
 Summarize absent reads.
 */
inline void summarize_abs(std::vector<softclip_info>& softclips,
                          std::unordered_map<std::string, abs_info*>& ump,
                          std::unordered_map<std::string, abs_info*>& umm,
                          const std::unordered_map<std::string, bool>& is_mainchr,
                          std::set<std::string>& tmp_set, std::string& tmp_str) {
    int v_size=softclips.size();
    softclip_info* s;
    for (int i=0; i < v_size; i++) {
        s=&(softclips[i]);
        if (s->breakpoint == 'N') { continue; }
        if ((s->clipend - s->clipstart) < DISCORDANT_READS_CLIP_LEN) { continue; }
        if (! is_mainchr.at(s->chr)) { continue; }
        if (s->is_reverse && (s->breakpoint == 'R')) {
            umm.at(s->chr)->R.push_back(s);
            umm.at(s->chr)->R_num++;
        } else if (s->is_reverse && (s->breakpoint == 'L')) {
            umm.at(s->chr)->L.push_back(s);
            umm.at(s->chr)->L_num++;
        } else if (s->breakpoint == 'R') {
            ump.at(s->chr)->R.push_back(s);
            ump.at(s->chr)->R_num++;
        } else {
            ump.at(s->chr)->L.push_back(s);
            ump.at(s->chr)->L_num++;
        }
        tmp_set.insert(s->chr);
    }
}


/*
 Output absent reads.
 */
inline bool output_abs(std::string _chr, std::unordered_map<std::string, abs_info*>& umr, std::unordered_map<std::string, abs_info*>& uml,
                       fstr* fs, char* qname, int32_t& l_qseq, bool& is_read2) {
    int64_t leng;
    bool written=false;
    if ((umr.at(_chr)->R_num > 0) && (uml.at(_chr)->L_num > 0)) {
        for (softclip_info* p1 : umr.at(_chr)->R) {
            for (softclip_info* p2 : uml.at(_chr)->L) {
                leng= ((int64_t)(p2->pos) - ((int64_t)(p1->pos) + (int64_t)(p1->rlen)));
                if ((leng >= ABS_MIN_DIST) && (leng <= ABS_MAX_DIST)) {
                    const char* _chrp= _chr.c_str();
                    std::fprintf(fs->ofs_abs, "%s/%d\t%s\t%ld\t%ld\t%s:%ld-%ld\t%s:%ld-%ld\t%d-%d\t%d-%d\n",
                                 qname, (is_read2 + 1), _chrp, (p1->pos + p1->rlen), p2->pos,
                                 _chrp, p1->pos, (p1->pos + p1->rlen), _chrp, p2->pos, (p2->pos + p2->rlen),
                                 p1->l_clip_len, (l_qseq - p1->r_clip_len), p2->l_clip_len, (l_qseq - p2->r_clip_len));
                    written=true;
                }
            }
        }
    }
    return written;
}


/*
 This is the core function to retrieve discordantly mapped reads.
 This:
    1) reads one line of a bam file
    2) parses cigar, SA-tag, and XA-tag
    3) outputs information if the read is discordantly mapped.
 */
inline void process_aln(htsFile *fp, sam_hdr_t *h, bam1_t *b, const std::vector<std::string>& mainchrs,
                        const std::unordered_map<std::string, bool>& is_mainchr,
                        const std::vector<uint32_t>& crepkmer, const ull& num_kmer,
                        std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                        std::vector<softclip_info>& softclips_sa,
                        std::vector<softclip_info>& softclips_xa,
                        std::string& tmp_str,
                        uint32_t*& cig_buf, size_t& cig_m,
                        std::string& fseq, std::string& rseq,
                        std::string& tmp_str1, std::string& tmp_str2, std::string& tmp_str3,
                        std::unordered_map<std::string, std::string>& um1,
                        std::unordered_map<std::string, std::string>& um2,
                        std::unordered_map<std::string, std::string>& um3,
                        std::unordered_map<std::string, std::string>::iterator& it1,
                        char* tmp_buf, std::set<std::string>& tmp_set, std::set<std::string>::iterator& it2,
                        fstr* fs, abs_info_manager& abs, read_stats& stats) {
    // remove 1) supplementary alignments, 2) single-end reads, 3) unmapped reads
    uint16_t& flag = b->core.flag;
    if (((flag & BAM_FSUPPLEMENTARY) > 0) == true) { return; }
    if (((flag & BAM_FPAIRED) > 0) == false) { return; }
    if (((flag & BAM_FUNMAP) > 0) == true) { return; }
    
    // prep
    uint32_t* cigar   = bam_get_cigar(b);  // cigar, 32-bit
    uint32_t& n_cigar = b->core.n_cigar;   // number of CIGAR operations
    bool contains_H =false;
    bool contains_S =false;
    bool contains_SA=false;
    bool contains_XA=false;
    parse_cigar(cigar, n_cigar, cigar_arr, contains_H, contains_S);
    uint8_t *sa_p    = bam_aux_get(b, "SA");
    uint8_t *xa_p    = bam_aux_get(b, "XA");
    if (sa_p) { contains_SA=true; }
    if (xa_p) { contains_XA=true; }
    
    bool is_distant_read=false;
    int64_t& isize   = b->core.isize;  // insertion size (beween R1 and R2)
    if ((isize == 0) || (isize <= -READ_PAIR_GAP_LEN) || (isize >= READ_PAIR_GAP_LEN)) {
        if ((((flag & BAM_FUNMAP) > 0) == false) && (((flag & BAM_FMUNMAP) > 0) == false)) { // both mapped
            is_distant_read=true;
        }
    }
    if ((! is_distant_read) && (! contains_S) && (! contains_SA) && (! contains_XA)) {
        return;  // neither chimeric nor distant read
    }
    
    bool is_reverse=false;
    if (((flag & BAM_FREVERSE) > 0) == true) { is_reverse=true; }
    char* chr        = h->target_name[b->core.tid];    // chr name
    int32_t &l_qseq  = b->core.l_qseq;                 // length of read
    int32_t l_clip_len, r_clip_len;
    int64_t clipstart, clipend;
    char breakpoint;
    int ret;  // return
    
    int64_t& start = b->core.pos;  // read mapping start, 0-based
    int64_t  rlen  = bam_cigar2rlen(n_cigar, cigar);
    
    softclips_sa.clear();
    softclips_xa.clear();
    
    // detect chimeric, softclip -> save as XA tag
    if (contains_S) {
        define_breakpoint(cigar_arr, n_cigar, breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq);
        softclips_xa.push_back(softclip_info(std::string(chr), start, is_reverse, breakpoint,
                                             l_clip_len, r_clip_len, clipstart, clipend, rlen));
    }
    
    // detect chimeric, SA-tag
    bool is_short_deletion=false;
    if (contains_SA) {
        *sa_p++;
        char* sa_tag=(char*)sa_p;  // simpler version of `char *bam_aux2Z(const uint8_t *s)` in htslib sam.c
        SA_parser(sa_tag, tmp_str, softclips_sa, cigar_arr, cig_buf, cig_m,
                  breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq);
        int vec_size=softclips_sa.size();
        // stop when short deletion rather than insertion
        for (int i=1; i < vec_size; i++) {
            if (softclips_xa[0].chr == softclips_sa[i].chr) {
                if ((-50 < (softclips_xa[0].pos - softclips_sa[i].pos)) && ((softclips_xa[0].pos - softclips_sa[i].pos) < 50)) {
                    is_short_deletion=true;
                    break;
                }
            }
        }
    }
    
    // detect chimeric, XA-tag
    if (contains_XA) {
        *xa_p++;
        char* xa_tag=(char*)xa_p;  // simpler version of `char *bam_aux2Z(const uint8_t *s)` in htslib sam.c
        XA_parser(xa_tag, tmp_str, softclips_xa, cigar_arr, cig_buf, cig_m,
                  breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq);
    }
    
    // prep seq
    char nt;
    fseq.clear();
    rseq.clear();
    uint8_t* seq_p = bam_get_seq(b);
    for (int i=0; i < l_qseq; i++) {
        nt=seq_nt16_str[bam_seqi(seq_p, i)];
        fseq += nt;
        rseq += complement[nt];  // just complement, do not reverse
    }
    std::reverse(rseq.begin(), rseq.end());  // do reverse
    int64_t end   = bam_endpos(b);            // read mapping end, 0-based
    char*   qname = bam_get_qname(b);   // read name
    bool is_read2 = (flag & BAM_FREAD2) > 0;
    char strand;
    if (is_reverse) { strand='-'; }
    else { strand='+'; }
    
    // output overhangs
    bool written;
    if (! is_short_deletion) {
        um1.clear();
        um2.clear();
        if (contains_S || contains_XA) {
            bool ret=output_pA(softclips_xa, crepkmer, num_kmer, chr, fseq, rseq, start, end, qname, l_qseq, is_read2, is_reverse, strand,
                               tmp_str, tmp_str1, tmp_str2, tmp_str3, um1, um2, tmp_buf, fs, is_mainchr);
            if (ret) { written=true; }
        }
        if (contains_SA) {
            bool ret=output_pA(softclips_sa, crepkmer, num_kmer, chr, fseq, rseq, start, end, qname, l_qseq, is_read2, is_reverse, strand,
                               tmp_str, tmp_str1, tmp_str2, tmp_str3, um1, um2, tmp_buf, fs, is_mainchr);
            if (ret) { written=true; }
        }
        if (written) { stats.pA++; }
        
        // output mapped seq and overhangs
        if (! um1.empty()) {  // um1 = (std::string, std::string) = (mapped seq, header)
            it1=um1.begin();
            while (it1 != um1.end()) {
                std::fprintf(fs->ofs_mapped, "%s\n%s\n", (it1->second).c_str(), (it1->first).c_str());
                it1++;
            }
        }
        if (! um2.empty()) {  // um1 = (std::string, std::string) = (clipped seq, header)
            it1=um2.begin();
            while (it1 != um2.end()) {
                std::fprintf(fs->ofs_overhang, "%s\n%s\n", (it1->second).c_str(), (it1->first).c_str());
                it1++;
            }
            stats.chimeric++;
        }
    }
    
    // output distant reads
    bool proceed_next=false;
    if (is_distant_read) {
        if (! contains_S) {
            proceed_next=true;
        } else if (softclips_xa[0].breakpoint == 'N') {
            proceed_next=true;
        } else if ((softclips_xa[0].clipend - softclips_xa[0].clipstart) < DISCORDANT_READS_CLIP_LEN) {
            proceed_next=true;
        }
    }
    if (proceed_next) {
        tmp_str.clear();  // distant read info
        if (is_mainchr.at(std::string(chr))) {
            std::sprintf(tmp_buf, "%s:%ld-%ld/%c;", chr, softclips_xa[0].pos, (softclips_xa[0].pos + softclips_xa[0].rlen), strand);
            tmp_str=tmp_buf;
        }
        if (contains_XA) {
            char _strand;
            for (softclip_info s : softclips_xa) {
                if ((s.breakpoint == 'N') || ((s.clipend - s.clipstart) < DISCORDANT_READS_CLIP_LEN)) {
                    if (s.is_reverse) { _strand='-'; }
                    else { _strand='+'; }
                    std::sprintf(tmp_buf, "%s:%ld-%ld/%c;", (s.chr).c_str(), s.pos, (s.pos + s.rlen), _strand);
                    tmp_str1.clear();  // distant read info
                    tmp_str1=tmp_buf;
                    tmp_str += tmp_str1;
                }
            }
        }
        if (! tmp_str.empty()) {
            tmp_str.pop_back();  // delete last ';'
            std::fprintf(fs->ofs_distant, "%s/%d\t%s\n", qname, (is_read2 + 1), tmp_str.c_str());
            stats.distant++;
        }
    }
    
    // summarize and output absent reads
    if (! contains_SA) { return; }
    if (! contains_XA) { return; }
    tmp_set.clear();
    summarize_abs(softclips_sa, abs.sa_p, abs.sa_m, is_mainchr, tmp_set, tmp_str);
    summarize_abs(softclips_xa, abs.xa_p, abs.xa_m, is_mainchr, tmp_set, tmp_str);
    written=false;
    for (it2=tmp_set.begin(); it2 != tmp_set.end(); it2++) {
        if (output_abs(*it2, abs.sa_p, abs.xa_p, fs, qname, l_qseq, is_read2)) { written=true; }
        if (output_abs(*it2, abs.xa_p, abs.sa_p, fs, qname, l_qseq, is_read2)) { written=true; }
        if (output_abs(*it2, abs.sa_m, abs.xa_m, fs, qname, l_qseq, is_read2)) { written=true; }
        if (output_abs(*it2, abs.xa_m, abs.sa_m, fs, qname, l_qseq, is_read2)) { written=true; }
        abs.clean_up_by_chr(*it2);  // clean up for next
    }
    if (written) { stats.absent++; }
}


/*
 This is a core function to manage file streams during threading.
    Before processing alignments, this catches available file streams generated in extract_discordant().
    After finishing processing of all reads in a chromosome, it releases the file streams so the next job can catch it.
 This is a wrapper of process_aln() that extract discordantly mapped reads.
 This reads one chr and processes all alignment by process_aln().
 */
int extract_discordant_per_chr(const char* f, int tid, const options& opts,
                               const std::vector<uint32_t>& crepkmer, const ull& num_kmer,
                               const std::vector<std::string>& mainchrs, const std::unordered_map<std::string, bool>& is_mainchr) {
    // take ofstream
    fstr* fs=fstrs.occupy();
    
    // open bam
    htsFile *fp=hts_open(f, "r");
    if (opts.is_cram) {
        hts_set_opt(fp, CRAM_OPT_REFERENCE, (opts.ref_fa).c_str());
    }
    sam_hdr_t *h=sam_hdr_read(fp);
    bam1_t *b= bam_init1();
    hts_idx_t *idx = nullptr;
    idx=sam_index_load(fp, f);
        
    // determine iter region
    const hts_pos_t beg = 0;  // chr start pos
    hts_pos_t end = h->target_len[tid];  // chr end pos
    hts_itr_t *iter = sam_itr_queryi(idx, tid, beg, end);
    if (iter == NULL) {   // region invalid or reference name not found
        std::cout << "ERROR at sam_itr_queryi(idx, tid, beg, end)" << std::endl;
        return -1;
    }
    
    // prepare for process_aln()
    std::pair<char, uint32_t> cigar_arr[MAX_CIGAR_LEN];
    std::vector<softclip_info> softclips_sa;
    std::vector<softclip_info> softclips_xa;
    softclips_sa.reserve(256);
    softclips_xa.reserve(256);
    uint32_t* cig_buf = nullptr;
    size_t cig_m=0;
    std::string fseq;   // seq of read, ATGCN
    std::string rseq;   // seq of read, ATGCN
    std::string tmp_str;
    std::string tmp_str1;
    std::string tmp_str2;
    std::string tmp_str3;
    std::unordered_map<std::string, std::string> um1;
    std::unordered_map<std::string, std::string> um2;
    std::unordered_map<std::string, std::string> um3;
    std::unordered_map<std::string, std::string>::iterator it1;
    std::set<std::string> tmp_set;
    std::set<std::string>::iterator it2;
    char* tmp_buf= new char[TMP_BUF_SIZE];
    abs_info_manager abs;
    abs.init();
    abs.emplace_chrs(mainchrs);
    read_stats stats;
    
//    std::cout << "processing " << h->target_name[tid] << " ..." << std::endl;
    
    // read bam or cram
    int ret;
    int processed_cnt=0;
    kstring_t aux={0, 0, NULL};
    while ((ret = sam_itr_next(fp, iter, b)) >= 0) {
        process_aln(fp, h, b, mainchrs, is_mainchr, crepkmer, num_kmer, cigar_arr, softclips_sa, softclips_xa,
                    tmp_str, cig_buf, cig_m, fseq, rseq,
                    tmp_str1, tmp_str2, tmp_str3, um1, um2, um3, it1, tmp_buf, tmp_set, it2, fs, abs, stats);
    }
    
    // output read stats
    std::fprintf(fs->ofs_stat, "%s\t%ld\t%ld\t%ld\t%ld\n", h->target_name[tid],
                 stats.pA, stats.chimeric, stats.distant, stats.absent);
    
    // finish up
    fstrs.release(fs);  // always release!!!
    delete tmp_buf;
    hts_itr_destroy(iter);
    sam_hdr_destroy(h);
    sam_close(fp);
    bam_destroy1(b);
    return 0;
}


/*
 This is a de facto main function and a core function to manage threading.
    This generates multiple threads and throws jobs in the threads.
    This also generates multiple file streams for threading.
    One thread will process one chromosome during one job.
    After finishing the job, it proceeds to the next chromosome.
 */
int extract_discordant(std::string bam, std::string f_mainchr, std::string mk, std::string mi, std::string ref_fa,
                       std::string ourdir, int n_thread, bool is_cram) {
    // parse args and options
    const options opts(bam, f_mainchr, mk, mi, ref_fa, ourdir, n_thread, is_cram);
    
    // prep to read file
    std::ifstream infile;
    std::string line;
    
    // load repeat .mk
    infile.open((opts.mi).c_str());
    if (! infile.is_open()) { return 1; }
    if (! std::getline(infile, line)) { return 1; }
    ull num_kmer=std::stoull(line);
    infile.close();
    std::cout << "Number of k-mers loading from " << opts.mk << ": " << num_kmer << std::endl;
    // load .mk
    infile.open((opts.mk).c_str(), std::ios::binary);
    if (! infile.is_open()) { return 1; }
    std::vector<uint32_t> repkmer;
    repkmer.resize(num_kmer);
    infile.read((char*)&repkmer[0], sizeof(uint32_t) * num_kmer);
    infile.close();
    const std::vector<uint32_t>& crepkmer=repkmer;
    // for 2bit conversion
    const int window_size=init_dna_to_2bit_32();
    
    // open bam
    const char *f= (opts.bam).c_str();
    htsFile *fp=hts_open(f, "r");
    if (opts.is_cram) {
        hts_set_opt(fp, CRAM_OPT_REFERENCE, (opts.ref_fa).c_str());
    }
    sam_hdr_t *h=sam_hdr_read(fp);
    bam1_t *b= bam_init1();
    
    // open bai
    hts_idx_t *idx = nullptr;
    idx=sam_index_load(fp, f);
    
    // sort chr for multiprocessing
    std::vector<int> sorted_chr;
    if (! opts.is_cram) {
        sort_chr_order(idx, h, sorted_chr);  // bam
    } else {
        sort_chr_order(h, sorted_chr);  // cram
    }
    std::cout << "n = " << sorted_chr.size() << " chrs will be scanned." << std::endl;
    
    // read mainchrs
    std::vector<std::string> mainchrs;
    std::set<std::string> mainchrs_set;
    infile.open((opts.f_mainchr).c_str());
    if (! infile.is_open()) { return 1; }
    while (std::getline(infile, line)) {
        mainchrs.push_back(line);
        mainchrs_set.insert(line);
    }
    infile.close();
    // judge whether mainchrs
    std::unordered_map<std::string, bool> is_mainchr;
    int32_t nseq = h->n_targets;    // number of chrs, cannot use idx info, because .crai does not hold this info
    for (int tid=0; tid < nseq; tid++) {
        std::string _chr=h->target_name[tid];
        if (mainchrs_set.find(_chr) != mainchrs_set.end()) {
            is_mainchr.emplace(_chr, true);
        } else {
            is_mainchr.emplace(_chr, false);
        }
    }
    
    // close bam, keep bai
    sam_hdr_destroy(h);
    sam_close(fp);
    bam_destroy1(b);
        
    // threading
    megane_thread_pool::ThreadPool mpool(opts.n_thread);
    std::vector<std::future<int>> results;
    
    // ofstreams (C fopen)
    for (int i=0; i < opts.n_thread; i++) {
        std::FILE* ofs_pA       =fopen((opts.outdir + std::string("/overhang_pA.txt") + std::to_string(i)).c_str(), "w");
        std::FILE* ofs_overhang =fopen((opts.outdir + std::string("/overhang.fa") + std::to_string(i)).c_str(), "w");
        std::FILE* ofs_distant  =fopen((opts.outdir + std::string("/distant.txt") + std::to_string(i)).c_str(), "w");
        std::FILE* ofs_mapped   =fopen((opts.outdir + std::string("/mapped.fa") + std::to_string(i)).c_str(), "w");
        std::FILE* ofs_abs      =fopen((opts.outdir + std::string("/absent.txt") + std::to_string(i)).c_str(), "w");
        std::FILE* ofs_stat     =fopen((opts.outdir + std::string("/stats.txt") + std::to_string(i)).c_str(), "w");
        if (ofs_pA       == nullptr) { return 1; }
        if (ofs_overhang == nullptr) { return 1; }
        if (ofs_distant  == nullptr) { return 1; }
        if (ofs_mapped   == nullptr) { return 1; }
        if (ofs_abs      == nullptr) { return 1; }
        if (ofs_stat     == nullptr) { return 1; }
        fstrs.push_back_new(ofs_pA, ofs_overhang, ofs_distant, ofs_mapped, ofs_abs, ofs_stat);
    }
    
    // per chr processing
    comp_init();  // make complementary table
    for (int tid : sorted_chr) {
        results.emplace_back(
            mpool.enqueue([=] {
                return extract_discordant_per_chr(f, tid, opts, crepkmer, num_kmer,
                                                  (const std::vector<std::string>)mainchrs,
                                                  (const std::unordered_map<std::string, bool>)is_mainchr);  // core func
            })
        );
    }
    
    for (auto && result: results) {
        result.get();
    }
    
    // close all file obj
    fstrs.close_fileobjs();
    
    return 0;
}


/*
 This is a main func for direct use.
 usage: %prog input.bam/cram main_chrs.txt input.mk output_dir n_thread [reference.fa]
 */
int main(int argc, char *argv[]) {
    // file check (not decent)
    bool is_cram=false;
    if (argc == 6) {
        // bam
    } else if (argc == 7) {
        // cram
        is_cram=true;
    } else {
        std::cerr << "Please specify required files." << std::endl;
        return 1;
    }
    
    std::string bam=argv[1];
    std::string f_mainchr=argv[2];
    std::string mk=argv[3];
    std::string mi=argv[3];
    mi.pop_back();
    mi += 'i';
    std::string outdir=argv[4];
    int n_thread= atoi(argv[5]);
    std::string ref_fa="";
    if (argc == 7) {
        ref_fa=argv[6];
    }
    
    std::cout << "bam " << bam << std::endl;
    std::cout << "f_mainchr " << f_mainchr << std::endl;
    std::cout << "mk " << mk << std::endl;
    std::cout << "mi " << mi << std::endl;
    std::cout << "outdir " << outdir << std::endl;
    std::cout << "n_thread " << n_thread << std::endl;
    std::cout << "ref_fa " << ref_fa << std::endl;
    std::cout << "is_cram " << is_cram << std::endl;
    
    int ret = extract_discordant(bam, f_mainchr, mk, mi, ref_fa, outdir, n_thread, is_cram);
    return ret;
}
