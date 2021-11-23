#include <iostream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <functional>
#include "htslib/sam.h"
#include "ThreadPool.h"

#define MAX_CIGAR_LEN 128
#define N_SA_INFO 6

typedef unsigned long ul;
typedef unsigned long long ull;


/*
 g++ -o extract_discordant -I /home/kooojiii/Desktop/htslib/htslib-1.13 -L /home/kooojiii/Desktop/htslib/htslib-1.13 extract_discordant.cpp -lhts -pthread -O2
 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/kooojiii/Desktop/htslib/htslib-1.13
 ./extract_discordant /home/kooojiii/Documents/testdata/bams/1kgp/GRCh38DH/NA12878.final.bam
 */



struct softclip_info {
    char*    chr=nullptr;
    int64_t  pos;
    bool     is_reverse;
    char     breakpoint;
    uint64_t clipstart;
    uint64_t clipend;
    
    softclip_info(char* chr, int64_t pos, bool is_reverse, char breakpoint, uint64_t clipstart, uint64_t clipend) {
        this->chr        = chr;
        this->pos        = pos;
        this->is_reverse = is_reverse;
        this->breakpoint = breakpoint;
        this->clipstart  = clipstart;
        this->clipend    = clipend;
    };
    
    // for debug, just a print function
    void print() {
        std::cout << chr << " "
        << pos << " "
        << is_reverse << " "
        << breakpoint << " "
        << clipstart << " "
        << clipend << std::endl;
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
inline int define_breakpoint(std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN], uint32_t& n_cigar,
                             char& breakpoint, int64_t& clipstart, int64_t& clipend, int32_t& l_qseq) {
    uint32_t left=0,right=0;
    // left length
    if (cigar_arr[0].first == 'S') {
        left=cigar_arr[0].second;
    }
    // right length
    if (cigar_arr[n_cigar - 1].first == 'S') {
        right=cigar_arr[n_cigar - 1].second;
    }
    // judge
    if (left == right) { return 1; }
    if (left > right) {
        breakpoint ='L';
        clipstart  = 0;
        clipend    = left;
    } else {
        breakpoint ='R';
        clipstart  = l_qseq - right;
        clipend    = l_qseq;
    }
    return 0;
}


/*
 Parse one SA entry in an SA tag.
 */
inline void parse_one_SA(std::string& tmp_sa,
                         std::vector<softclip_info>& softclips,
                         std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                         uint32_t* cig_buf, size_t& cig_m,
                         char& breakpoint, int64_t& clipstart, int64_t& clipend,
                         int32_t& l_qseq, bool& is_reverse) {
    // split SA info by ','
    int comma1=0,comma2=0,comma3=0,comma4=0;
    int n=0, pos=0;
    for (char c : tmp_sa) {
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
    uint32_t _n_cigar=sam_parse_cigar((const char*)(tmp_sa.substr(comma3 + 1, comma4 - comma3 - 1)).c_str(),
                                      nullptr, &cig_buf, &cig_m);
    parse_cigar(cig_buf, _n_cigar, cigar_arr);
    
    // save breakpoint info
    int ret=define_breakpoint(cigar_arr, _n_cigar, breakpoint, clipstart, clipend, l_qseq);
    if (ret == 0) {
        std::string _chr=tmp_sa.substr(0, comma1);
        int64_t _pos=std::stoll(tmp_sa.substr(comma1 + 1, comma2 - comma1 - 1)) - 1;  // 0-based
        char strand=tmp_sa[comma2 + 1];
        bool _is_reverse;
        if ((is_reverse == false) && (strand == '+')) {
            _is_reverse=false;
        } else if ((is_reverse == false) && (strand == '-')) {
            _is_reverse=true;
        } else if ((is_reverse == true) && (strand == '+')) {
            _is_reverse=true;
        } else {
            _is_reverse=false;
        }
        softclips.push_back(softclip_info((char*)_chr.c_str(), _pos, _is_reverse, breakpoint, clipstart, clipend));
    }
//    for (softclip_info s : softclips) {
//        s.print();
//    }
//    std::exit(0);
}


/*
 Parses char SA tag and stores in std::vector<softclip_info>
 */
inline void SA_parser(char* sa_tag, std::string& tmp_sa,
                      std::vector<softclip_info>& softclips,
                      std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                      uint32_t* cig_buf, size_t& cig_m,
                      char& breakpoint, int64_t& clipstart, int64_t& clipend,
                      int32_t& l_qseq, bool& is_reverse) {
//    std::cout << sa_tag << std::endl;
    // count SA entires
    int n_entry=0;
    int sa_len=0;
    tmp_sa.clear();
    while (sa_tag[sa_len]) {
        if (sa_tag[sa_len] == ';') {
            parse_one_SA(tmp_sa, softclips, cigar_arr, cig_buf, cig_m,
                         breakpoint, clipstart, clipend, l_qseq, is_reverse);
            tmp_sa.clear();
        } else {
            tmp_sa += sa_tag[sa_len];
        }
        sa_len++;
    }
}


/*
 Parse one XA entry in an XA tag.
 */
inline void parse_one_XA(std::string& tmp_sa,
                         std::vector<softclip_info>& softclips,
                         std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                         uint32_t* cig_buf, size_t& cig_m,
                         char& breakpoint, int64_t& clipstart, int64_t& clipend,
                         int32_t& l_qseq, bool& is_reverse) {
    // split SA info by ','
    int comma1=0,comma2=0,comma3=0;
    int n=0, pos=0;
    for (char c : tmp_sa) {
        if (c == ',') {
            if (n == 0) {
                comma1=pos;
            } else if (n == 1) {
                comma2=pos;
            } else {
                comma3=pos;
            }
            n++;
        }
        pos++;
    }
    
    // parse cigar
    cig_m=0;
    uint32_t _n_cigar=sam_parse_cigar((const char*)(tmp_sa.substr(comma2 + 1, comma3 - comma2 - 1)).c_str(),
                                      nullptr, &cig_buf, &cig_m);
    parse_cigar(cig_buf, _n_cigar, cigar_arr);
    
    // save breakpoint info
    int ret=define_breakpoint(cigar_arr, _n_cigar, breakpoint, clipstart, clipend, l_qseq);
    if (ret == 0) {
        std::string _chr=tmp_sa.substr(0, comma1);
        int64_t _pos=std::stoll(tmp_sa.substr(comma1 + 2, comma2 - comma1 - 2)) - 1;  // 0-based
        char strand=tmp_sa[comma1 + 1];
        bool _is_reverse;
        if ((is_reverse == false) && (strand == '+')) {
            _is_reverse=false;
        } else if ((is_reverse == false) && (strand == '-')) {
            _is_reverse=true;
        } else if ((is_reverse == true) && (strand == '+')) {
            _is_reverse=true;
        } else {
            _is_reverse=false;
        }
        softclips.push_back(softclip_info((char*)_chr.c_str(), _pos, _is_reverse, breakpoint, clipstart, clipend));
    }
//    for (softclip_info s : softclips) {
//        s.print();
//    }
//    std::exit(0);
}


/*
 Parses char XA tag and stores in std::vector<softclip_info>
 */
inline void XA_parser(char* sa_tag, std::string& tmp_sa,
                      std::vector<softclip_info>& softclips,
                      std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                      uint32_t* cig_buf, size_t& cig_m,
                      char& breakpoint, int64_t& clipstart, int64_t& clipend,
                      int32_t& l_qseq, bool& is_reverse) {
//    std::cout << sa_tag << std::endl;
    // count SA entires
    int n_entry=0;
    int sa_len=0;
    tmp_sa.clear();
    while (sa_tag[sa_len]) {
        if (sa_tag[sa_len] == ';') {
            parse_one_XA(tmp_sa, softclips, cigar_arr, cig_buf, cig_m,
                         breakpoint, clipstart, clipend, l_qseq, is_reverse);
            tmp_sa.clear();
        } else {
            tmp_sa += sa_tag[sa_len];
        }
        sa_len++;
    }
}


/*
 Core function to judge discordant reads.
 */
void inline process_aln(htsFile *fp, sam_hdr_t *h, bam1_t *b,
                        std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                        std::vector<softclip_info>& softclips,
                        std::string& tmp_sa,
                        uint32_t* cig_buf, size_t& cig_m) {
    // remove 1) supplementary alignments, 2) single-end reads, 3) unmapped reads
    uint16_t& flag = b->core.flag;
    if (((flag & BAM_FSUPPLEMENTARY) > 0) == true) { return; }
    if (((flag & BAM_FPAIRED) > 0) == false) { return; }
    if (((flag & BAM_FUNMAP) > 0) == true) { return; }
    
    // prep for retrieving overhangs
    // if (('S' in ls[5]) or ('SA:Z:' in line) or ('XA:Z:' in line)) and not ('H' in ls[5]):
    uint32_t* cigar   = bam_get_cigar(b);  // cigar, 32-bit
    uint32_t& n_cigar = b->core.n_cigar;   // number of CIGAR operations
    bool contains_H=false;
    bool contains_S=false;
    bool contains_SA=false;
    bool contains_XA=false;
    parse_cigar(cigar, n_cigar, cigar_arr, contains_H, contains_S);
    if (contains_H) { return; }
    uint8_t *sa_p    = bam_aux_get(b, "SA");
    uint8_t *xa_p    = bam_aux_get(b, "XA");
    if (sa_p) { contains_SA=true; }
    if (xa_p) { contains_XA=true; }
    
    bool is_reverse=false;
    if (((flag & BAM_FREVERSE) > 0) == true) { is_reverse=true; }
    char* chr        = h->target_name[b->core.tid];    // chr name
    int32_t &l_qseq  = b->core.l_qseq;                 // length of read
    int64_t clipstart, clipend;
    char breakpoint;
    int ret;  // return
    
    // retrieve overhangs, softclip
    if (contains_S) {
        ret=define_breakpoint(cigar_arr, n_cigar, breakpoint, clipstart, clipend, l_qseq);
        if (ret == 0) {
            softclips.push_back(softclip_info(chr, b->core.pos, is_reverse, breakpoint, clipstart, clipend));
        }
    }
    
    // retrieve overhangs, SA-tag
    bool is_short_deletion=false;
    if (contains_SA) {
        *sa_p++;
        char* sa_tag=(char*)sa_p;  // simpler version of `char *bam_aux2Z(const uint8_t *s)` in htslib sam.c
        SA_parser(sa_tag, tmp_sa, softclips, cigar_arr, cig_buf, cig_m,
                  breakpoint, clipstart, clipend, l_qseq, is_reverse);
        int vec_size=softclips.size();
        // stop when short deletion rather than insertion
        if (vec_size >= 2) {
            for (int i=1; i < vec_size; i++) {
                if (*(softclips[0].chr) == *(softclips[i].chr)) {
                    if ((-50 < (softclips[0].pos - softclips[i].pos)) && ((softclips[0].pos - softclips[i].pos) < 50)) {
                        is_short_deletion=true;
                        break;
                    }
                }
            }
        }
    }
    
    // retrieve overhangs, XA-tag
    if ((! is_short_deletion) && (contains_XA)) {
        *xa_p++;
        char* xa_tag=(char*)xa_p;  // simpler version of `char *bam_aux2Z(const uint8_t *s)` in htslib sam.c
        XA_parser(xa_tag, tmp_sa, softclips, cigar_arr, cig_buf, cig_m,
                  breakpoint, clipstart, clipend, l_qseq, is_reverse);
    }
//    std::cout << softclips.size() << std::endl;
}



int extract_discordant_per_chr(char* f, hts_idx_t *idx, int tid) {
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
    std::vector<softclip_info> softclips;
    softclips.reserve(256);
    std::string tmp_sa;
    uint32_t* cig_buf = nullptr;
    size_t cig_m=0;
    
    // read bam or cram
    int ret;
    int processed_cnt=0;
    kstring_t aux={0, 0, NULL};
    while ((ret = sam_itr_next(fp, iter, b)) >= 0) {
        process_aln(fp, h, b, cigar_arr, softclips, tmp_sa, cig_buf, cig_m);
        softclips.clear();
//        if (++processed_cnt >= 1) {
//            break;
//        }
    }
    hts_itr_destroy(iter);
    return 0;
}



/*
 read bam or cram
 */
int extract_discordant(int argc, char *argv[]) {
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
    ThreadPool pool(1);
    std::vector<std::future<int>> results;
    
    // per chr processing
    int i=0;
    for (int tid : sorted_chr) {
        results.emplace_back(
            pool.enqueue([=] {
                return extract_discordant_per_chr(f, idx, tid);
            })
        );
//        i++;
//        if (i >= 1) {
//            break;
//        }
    }
    
    for (auto && result: results) {
        result.get();
    }
    
    return 0;
}



/*
 main func for direct use
 usage: %prog input.bam
 */
int main(int argc, char *argv[]) {
    int ret = extract_discordant(argc, argv);
}
