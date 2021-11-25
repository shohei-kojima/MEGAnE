#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include "htslib/sam.h"
#include "dna_to_2bit.hpp"
#include "complementary_seq.hpp"
#include "ThreadPool.h"
using namespace dna_to_2bit_hpp;
using namespace complementary_seq_hpp;

#define MAX_CIGAR_LEN 128
#define TMP_BUF_SIZE  131072
#define READ_PAIR_GAP_LEN 2000
#define MAX_SEQ_LEN 512
#define DISCORDANT_READS_CLIP_LEN 20
#define REP_KMER_LEN 8
#define MAPPED_REGION_LOW_COMPLEX_THRESHOLD 0.7
#define POLYA_OVERHANG_THRESHOLD 0.7

typedef unsigned long ul;
typedef unsigned long long ull;

const char ATGC[]="ATGC";


/*
 g++ -o extract_discordant -I /home/kooojiii/Desktop/htslib/htslib-1.13 -L /home/kooojiii/Desktop/htslib/htslib-1.13 extract_discordant.cpp -lhts -pthread -O2
 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/kooojiii/Desktop/htslib/htslib-1.13
 time ./extract_discordant /home/kooojiii/Documents/testdata/bams/1kgp/GRCh38DH/NA12878.final.bam
 */



/*
 Struct to keep file streams during threading
 */
struct fstr {
    std::FILE* ofs_pA;
    std::FILE* ofs_overhang;
    std::FILE* ofs_distant;
    std::FILE* ofs_mapped;
    std::FILE* ofs_unmapped;
    std::FILE* ofs_abs;
    std::FILE* ofs_stat;
    bool is_occupied;
    
    fstr(std::FILE* ofs_pA, std::FILE* ofs_overhang, std::FILE* ofs_distant, std::FILE* ofs_mapped,
         std::FILE* ofs_unmapped, std::FILE* ofs_abs, std::FILE* ofs_stat, bool is_occupied) {
        this->ofs_pA=ofs_pA;
        this->ofs_overhang=ofs_overhang;
        this->ofs_distant=ofs_distant;
        this->ofs_mapped=ofs_mapped;
        this->ofs_unmapped=ofs_unmapped;
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
                       std::FILE* ofs_unmapped, std::FILE* ofs_abs, std::FILE* ofs_stat) {
        fstr* fstrp = new fstr(ofs_pA, ofs_overhang, ofs_distant, ofs_mapped,
                               ofs_unmapped, ofs_abs, ofs_stat, false);
        this->vec.push_back(fstrp);
    }
    
    fstr* occupy() { 
        std::unique_lock<std::mutex> lock(this->_mutex);  // need mutex
        fstr* fstrp=nullptr;
        for (fstr* t : this->vec) {
            if (t->is_occupied == false) {  // if not the ofstream is used by the threads
                t->is_occupied=true;  // take this ofstream by a thread
                fstrp=t;
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
            fclose(t->ofs_unmapped);
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
 Class cfstrs that will be accessed from multiple threads
 */
static cfstrs fstrs;


/*
 Stores info of clipped reads, SA tags, XA tags.
 */
struct softclip_info {
    char*    chr=nullptr;
    int64_t  pos;
    bool     is_reverse;
    char     breakpoint;
    int32_t  l_clip_len;
    int32_t  r_clip_len;
    uint64_t clipstart;
    uint64_t clipend;
    uint64_t rlen;
    
    softclip_info(char* chr, int64_t pos, bool is_reverse, char breakpoint,
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
        std::cout << chr << " " << pos << " "
        << is_reverse << " " << breakpoint << " "
        << l_clip_len << " " << r_clip_len << " "
        << clipstart << " " << clipend << " "
        << rlen << std::endl;
    }
};


/*
 Returns chr names with at least one read.
 Chr names are in vector sorted by number of reads with descending manner.
 */
void sort_chr_order(hts_idx_t *idx, sam_hdr_t *h, std::vector<int> &sorted_chr) {
    // retrieve mapped read counts from bam index
    int nseq = hts_idx_nseq(idx);    // number of chrs
    uint64_t mapped,umapped;
    std::vector<std::pair<int, int>> mapped_counts;  // pair<mapped_read_num, tid>
    mapped_counts.reserve(nseq);
    for (int tid=0; tid < nseq; tid++) {
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
 Parses cigar and stores in an array. Also checks whether the cigar contains H and S.
 Array type is std::pair<char, uint32_t>, which is <cigar_opchr, bam_cigar_oplen>
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
 Array type is std::pair<char, uint32_t>, which is <cigar_opchr, bam_cigar_oplen>
 */
inline void parse_cigar(uint32_t* cigar, uint32_t& n_cigar,
                        std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN]) {
    for (int i=0; i < n_cigar; ++i) {
        cigar_arr[i]=std::make_pair(bam_cigar_opchr(cigar[i]),
                                    bam_cigar_oplen(cigar[i]));
    }
}


/*
 Define breakpoint
 Returns 0 if breakpoint was determined.
 Returns 1 if breakpoint was not determined.
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
 Parse one SA entry in an SA tag.
 */
inline void parse_one_SA(std::string& tmp_str,
                         std::vector<softclip_info>& softclips,
                         std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                         uint32_t*& cig_buf, size_t& cig_m,
                         char& breakpoint, int32_t& l_clip_len, int32_t& r_clip_len,
                         int64_t& clipstart, int64_t& clipend,
                         int32_t& l_qseq, bool& is_reverse) {
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
    uint32_t _n_cigar = sam_parse_cigar((const char*)(tmp_str.substr(comma3 + 1, comma4 - comma3 - 1)).c_str(),
                                      nullptr, &cig_buf, &cig_m);
    int64_t _rlen     = bam_cigar2rlen(_n_cigar, cig_buf);
    parse_cigar(cig_buf, _n_cigar, cigar_arr);
    
    // save breakpoint info
    define_breakpoint(cigar_arr, _n_cigar, breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq);
    std::string _chr=tmp_str.substr(0, comma1);
    int64_t _pos=std::stoll(tmp_str.substr(comma1 + 1, comma2 - comma1 - 1)) - 1;  // 0-based
    char _strand=tmp_str[comma2 + 1];
    bool _is_reverse;
    if ((is_reverse == false) && (_strand == '+')) {
        _is_reverse=false;
    } else if ((is_reverse == false) && (_strand == '-')) {
        _is_reverse=true;
    } else if ((is_reverse == true) && (_strand == '+')) {
        _is_reverse=true;
    } else {
        _is_reverse=false;
    }
    softclips.push_back(softclip_info((char*)_chr.c_str(), _pos, _is_reverse, breakpoint,
                                      l_clip_len, r_clip_len, clipstart, clipend, _rlen));
//    for (softclip_info s : softclips) {
//        s.print();
//    }
//    std::exit(0);
}


/*
 Parses char SA tag and stores in std::vector<softclip_info>
 */
inline void SA_parser(char* sa_tag, std::string& tmp_str,
                      std::vector<softclip_info>& softclips,
                      std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                      uint32_t*& cig_buf, size_t& cig_m,
                      char& breakpoint, int32_t& l_clip_len, int32_t& r_clip_len,
                      int64_t& clipstart, int64_t& clipend,
                      int32_t& l_qseq, bool& is_reverse) {
//    std::cout << sa_tag << std::endl;
    // count SA entires
    int n_entry=0;
    int sa_len=0;
    tmp_str.clear();
    while (sa_tag[sa_len]) {
        if (sa_tag[sa_len] == ';') {
            parse_one_SA(tmp_str, softclips, cigar_arr, cig_buf, cig_m,
                         breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq, is_reverse);
            tmp_str.clear();
        } else {
            tmp_str += sa_tag[sa_len];
        }
        sa_len++;
    }
}


/*
 Parse one XA entry in an XA tag.
 */
inline void parse_one_XA(std::string& tmp_str,
                         std::vector<softclip_info>& softclips,
                         std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                         uint32_t*& cig_buf, size_t& cig_m,
                         char& breakpoint, int32_t& l_clip_len, int32_t& r_clip_len,
                         int64_t& clipstart, int64_t& clipend,
                         int32_t& l_qseq, bool& is_reverse) {
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
    uint32_t _n_cigar = sam_parse_cigar((const char*)(tmp_str.substr(comma2 + 1, comma3 - comma2 - 1)).c_str(),
                                        nullptr, &cig_buf, &cig_m);
    int64_t _rlen     = bam_cigar2rlen(_n_cigar, cig_buf);
    parse_cigar(cig_buf, _n_cigar, cigar_arr);
    
    // save breakpoint info
    define_breakpoint(cigar_arr, _n_cigar, breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq);
    std::string _chr=tmp_str.substr(0, comma1);
    int64_t _pos=std::stoll(tmp_str.substr(comma1 + 2, comma2 - comma1 - 2)) - 1;  // 0-based
    char _strand=tmp_str[comma1 + 1];
    bool _is_reverse;
    if ((is_reverse == false) && (_strand == '+')) {
        _is_reverse=false;
    } else if ((is_reverse == false) && (_strand == '-')) {
        _is_reverse=true;
    } else if ((is_reverse == true) && (_strand == '+')) {
        _is_reverse=true;
    } else {
        _is_reverse=false;
    }
    softclips.push_back(softclip_info((char*)_chr.c_str(), _pos, _is_reverse, breakpoint,
                                      l_clip_len, r_clip_len, clipstart, clipend, _rlen));
//    for (softclip_info s : softclips) {
//        s.print();
//    }
//    std::exit(0);
//    std::cout << softclips.size() << std::endl;
}


/*
 Parses char XA tag and stores in std::vector<softclip_info>
 */
inline void XA_parser(char* sa_tag, std::string& tmp_str,
                      std::vector<softclip_info>& softclips,
                      std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                      uint32_t*& cig_buf, size_t& cig_m,
                      char& breakpoint, int32_t& l_clip_len, int32_t& r_clip_len,
                      int64_t& clipstart, int64_t& clipend,
                      int32_t& l_qseq, bool& is_reverse) {
//    std::cout << sa_tag << std::endl;
    // count SA entires
    int n_entry=0;
    int sa_len=0;
    tmp_str.clear();
    while (sa_tag[sa_len]) {
        if (sa_tag[sa_len] == ';') {
            parse_one_XA(tmp_str, softclips, cigar_arr, cig_buf, cig_m,
                         breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq, is_reverse);
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
inline ull custom_binary_search(const std::vector<uint16_t>& vec, const ull& vec_size, uint16_t& key) {
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
 This converts DNA to 8-nt 2bit and judges whether the 8-mer is repeat-derived.
 Args:
    1) seq
    2) length of the seq
    3) vector containing repeat k-mer set
    4) vector size
 */
inline bool is_rep_overhang(std::string seq, uint64_t clip_len, const std::vector<uint16_t>& crepkmer, const ull& num_kmer) {
    ull pos=0;
    // first window_size (window_size = REP_KMER_LEN)
    uint16_t bit2f=0;
    int nn=0;
    for (int i=0; i < REP_KMER_LEN; i++) {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_16[seq[i]];
        if (seq[i] == 'N' || seq[i] == 'n') {  // ignore when N or n appears
            nn= REP_KMER_LEN - 1;
        } else if (nn > 0) {  // within window_size-nt from N or n
            nn--;
        }
    }
    pos=custom_binary_search(crepkmer, num_kmer, bit2f);
    if (crepkmer[pos] == bit2f) { return true; }
    
    // rolling calc.
    for (int i=REP_KMER_LEN; i < clip_len; i++) {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_16[seq[i]];
        if (seq[i] == 'N' || seq[i] == 'n') {  // ignore when N or n appears
            nn= REP_KMER_LEN - 1;
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
 This 1) output pA and 2) stores overhang info and mapped seq to unordered_maps
 */
inline void output_pA(std::vector<softclip_info>& softclips, const std::vector<uint16_t>& crepkmer, const ull& num_kmer,
                     char* chr, std::string& fseq, std::string& rseq, int64_t& start, int64_t& end,
                     char* qname, int32_t& l_qseq, bool& is_read2, bool& is_reverse, char& strand,
                     std::string& tmp_str, std::string& tmp_str1, std::string& tmp_str2, std::string& tmp_str3,
                     std::unordered_map<std::string, std::string>& um1,
                     std::unordered_map<std::string, std::string>& um2,
                     char* tmp_buf, fstr* fs) {
    for (softclip_info s : softclips) {
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
            if (! is_simple_repeat(tmp_str1, mapped_len)) {
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
    }
}


/*
 Core function to judge discordant reads.
 */
inline void process_aln(htsFile *fp, sam_hdr_t *h, bam1_t *b, const std::vector<uint16_t>& crepkmer, const ull& num_kmer,
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
                        char* tmp_buf,
                        fstr* fs) {
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
        is_distant_read=true;
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
        softclips_xa.push_back(softclip_info(chr, start, is_reverse, breakpoint,
                                             l_clip_len, r_clip_len, clipstart, clipend, rlen));
    }
    
    // detect chimeric, SA-tag
    bool is_short_deletion=false;
    if (contains_SA) {
        *sa_p++;
        char* sa_tag=(char*)sa_p;  // simpler version of `char *bam_aux2Z(const uint8_t *s)` in htslib sam.c
        SA_parser(sa_tag, tmp_str, softclips_sa, cigar_arr, cig_buf, cig_m,
                  breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq, is_reverse);
        int vec_size=softclips_sa.size();
        // stop when short deletion rather than insertion
        for (int i=1; i < vec_size; i++) {
            if (*(softclips_xa[0].chr) == *(softclips_sa[i].chr)) {
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
                  breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq, is_reverse);
    }
//    std::cout << softclips_sa.size() << " " << softclips_xa.size() << std::endl;
    
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
//    std::cout << fseq << std::endl;
//    std::cout << rseq << std::endl;
    int64_t end   = bam_endpos(b);            // read mapping end, 0-based
    char*   qname = bam_get_qname(b);   // read name
    bool is_read2 = (flag & BAM_FREAD2) > 0;
    char strand;
    if (is_reverse) { strand='-'; }
    else { strand='+'; }
    
    // output overhangs
    um1.clear();
    um2.clear();
    if (! contains_H) {
        if (contains_S || contains_XA) {
            output_pA(softclips_xa, crepkmer, num_kmer, chr, fseq, rseq, start, end, qname, l_qseq, is_read2, is_reverse, strand,
                      tmp_str, tmp_str1, tmp_str2, tmp_str3, um1, um2, tmp_buf, fs);
        }
        if (contains_SA) {
            output_pA(softclips_sa, crepkmer, num_kmer, chr, fseq, rseq, start, end, qname, l_qseq, is_read2, is_reverse, strand,
                      tmp_str, tmp_str1, tmp_str2, tmp_str3, um1, um2, tmp_buf, fs);
        }
    }
    
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
    }
}


int extract_discordant_per_chr(char* f, hts_idx_t *idx, int tid, const std::vector<uint16_t>& crepkmer, const ull& num_kmer) {
    // take ofstream
    fstr* fs=fstrs.occupy();
    
    // open bam
    htsFile *fp=hts_open(f, "r");
    sam_hdr_t *h=sam_hdr_read(fp);
    bam1_t *b= bam_init1();
    
    std::cout << "processing " << h->target_name[tid] << " ..." << std::endl;
    
    // determine iter region
    const hts_pos_t beg = 0;  // chr start pos
    hts_pos_t end = h->target_len[tid];  // chr end pos
    hts_itr_t *iter = sam_itr_queryi(idx, tid, beg, end);
    if (iter == NULL) {   // region invalid or reference name not found
        std::cout << "ERROR" << std::endl;
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
    char* tmp_buf= new char[TMP_BUF_SIZE];
    
    // read bam or cram
    int ret;
    int processed_cnt=0;
    kstring_t aux={0, 0, NULL};
    while ((ret = sam_itr_next(fp, iter, b)) >= 0) {
        process_aln(fp, h, b, crepkmer, num_kmer, cigar_arr, softclips_sa, softclips_xa, tmp_str, cig_buf, cig_m, fseq, rseq,
                    tmp_str1, tmp_str2, tmp_str3, um1, um2, um3, it1, tmp_buf,
                    fs);
//        if (++processed_cnt >= 1) {
//            break;
//        }
    }
    
    // finish up
    fstrs.release(fs);  // always release!!!
    delete tmp_buf;
    hts_itr_destroy(iter);
    return 0;
}


/*
 read bam or cram
 */
int extract_discordant(int argc, char *argv[]) {
    // load repeat .mk
    const char* f_mk="/home/kooojiii/results/2021/misc/MEGAnE_test/211122_1/reshaped_repbase.mk";
    const char* f_mi="/home/kooojiii/results/2021/misc/MEGAnE_test/211122_1/reshaped_repbase.mi";
    std::ifstream infile;
    std::string line;
    infile.open(f_mi);
    if (! infile.is_open()) { return 1; }
    if (! getline(infile, line)) { return 1; }
    ull num_kmer=std::stoull(line);
    infile.close();
    std::cout << "Number of k-mers loading from " << f_mk << ": " << num_kmer << std::endl;
    // load .mk
    infile.open(f_mk, std::ios::binary);
    if (! infile.is_open()) { return 1; }
    std::vector<uint16_t> repkmer;
    repkmer.resize(num_kmer);
    infile.read((char*)&repkmer[0], sizeof(uint16_t) * num_kmer);
    infile.close();
    const std::vector<uint16_t>& crepkmer=repkmer;
    // for 2bit conversion
    const int window_size=init_dna_to_2bit_16();
    
    // open bam
    char *f=argv[1];
    htsFile *fp=hts_open(f, "r");
    sam_hdr_t *h=sam_hdr_read(fp);
    
    // open bai
    hts_idx_t *idx = NULL;
    idx=sam_index_load(fp, f);
    
    // sort chr for multiprocessing
    std::vector<int> sorted_chr;
    sort_chr_order(idx, h, sorted_chr);
    std::cout << "n = " << sorted_chr.size() << " chrs will be scanned." << std::endl;
    
    // close bam, keep bai
    sam_hdr_destroy(h);
    sam_close(fp);
    
    // threading
    const int thread_n=1;
    ThreadPool pool(thread_n);
    std::vector<std::future<int>> results;
    
    // ofstreams (C fopen)
    for (int i=0; i < thread_n; i++) {
        std::FILE* ofs_pA       =fopen("_tmp_pA.txt", "w");
        std::FILE* ofs_overhang =fopen("_tmp_overhang.txt", "w");
        std::FILE* ofs_distant  =fopen("_tmp_distant.txt", "w");
        std::FILE* ofs_mapped   =fopen("_tmp_mapped.txt", "w");
        std::FILE* ofs_unmapped =fopen("_tmp_unmapped.txt", "w");
        std::FILE* ofs_abs      =fopen("_tmp_abs.txt", "w");
        std::FILE* ofs_stat     =fopen("_tmp_stat.txt", "w");
        if (ofs_pA       == nullptr) { return 1; }
        if (ofs_overhang == nullptr) { return 1; }
        if (ofs_distant  == nullptr) { return 1; }
        if (ofs_mapped   == nullptr) { return 1; }
        if (ofs_unmapped == nullptr) { return 1; }
        if (ofs_abs      == nullptr) { return 1; }
        if (ofs_stat     == nullptr) { return 1; }
        fstrs.push_back_new(ofs_pA, ofs_overhang, ofs_distant, ofs_mapped,
                            ofs_unmapped, ofs_abs, ofs_stat);
    }
    
    // per chr processing
    comp_init();  // make complementary table
    int i=0;
    for (int tid : sorted_chr) {
        results.emplace_back(
            pool.enqueue([=] {
                return extract_discordant_per_chr(f, idx, 21 /*tid*/, crepkmer, num_kmer);  // core func
            })
        );
        i++;
        if (i >= 1) {
            break;
        }
    }
    
    for (auto && result: results) {
        result.get();
    }
    
    // close all file obj
    fstrs.close_fileobjs();
    
    return 0;
}


/*
 main func for direct use
 usage: %prog input.bam
 */
int main(int argc, char *argv[]) {
    int ret = extract_discordant(argc, argv);
}
