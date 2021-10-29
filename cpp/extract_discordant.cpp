#include <iostream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <functional>
#include "htslib/sam.h"
#include "ThreadPool.h"
using namespace std;


/*
 g++ -o extract_discordant -I /home/kooojiii/Desktop/htslib/htslib-1.13 -L /home/kooojiii/Desktop/htslib/htslib-1.13 extract_discordant.cpp -lhts -pthread
 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/kooojiii/Desktop/htslib/htslib-1.13
 ./extract_discordant /home/kooojiii/Documents/testdata/bams/1kgp/GRCh38DH/NA12878.final.bam
 */


/*
 Returns chr names with at least one read.
 Chr names are in vector sorted by number of reads with descending manner.
 */
void sort_chr_order(hts_idx_t *idx, sam_hdr_t *h, vector<int> &sorted_chr) {
    // retrieve mapped read counts from bam index
    int nseq = hts_idx_nseq(idx);    // number of chrs
    uint64_t mapped,umapped;
    vector<pair<int, int>> mapped_counts;  // pair<mapped_read_num, tid>
    mapped_counts.reserve(nseq);
    for (int tid=0; tid < nseq; tid++) {
        hts_idx_get_stat(idx, tid, &mapped, &umapped);  // return can be -1 when no read on chr
        if (mapped >= 1) {
            mapped_counts.push_back(make_pair(mapped, tid));
        }
    }
    
    // ascending by mapped read counts
    sort(mapped_counts.begin(), mapped_counts.end());
    int vec_size=mapped_counts.size();
    
    // convert to desceding order
    sorted_chr.reserve(vec_size);
    for (int i= vec_size - 1; i >= 0; i--) {
        sorted_chr.push_back(mapped_counts[i].second);
    }
}


/*
 Core function to judge discordant reads.
 */
int inline process_aln(htsFile *fp, sam_hdr_t *h, bam1_t *b) {
    hts_pos_t &start  = b->core.pos;                    // left position, 0-based
    hts_pos_t end     = bam_endpos(b);                  // right position, 0-based
    char     *qname   = bam_get_qname(b);               // read name

    // format seq
    int32_t &l_qseq   = b->core.l_qseq;                 // length of read
    uint8_t *tmp_s    = bam_get_seq(b);                 // seq of read, nt16
    char *seq=(char *)malloc(l_qseq + 1);               // seq of read, ATGCN
    for (int i=0; i < l_qseq; i++) {
        seq[i]=seq_nt16_str[bam_seqi(tmp_s, i)];        // get nucleotide id and convert into IUPAC id
    }
    seq[l_qseq]='\0';
    
//    cout
//        << qname << " "
//        << start << " "
//        << end << " "
//        << seq << " "
//        << endl;
    delete(seq);
    return 0;
}



int extract_discordant_per_chr(char* f, hts_idx_t *idx, int tid) {
    // open bam
    htsFile *fp=hts_open(f, "r");
    sam_hdr_t *h=sam_hdr_read(fp);
    bam1_t *b= bam_init1();
    
    cout << "processing " << h->target_name[tid] << " ..." << endl;
    
    // determine iter region
    const hts_pos_t beg = 0;  // chr start pos
    hts_pos_t end = h->target_len[tid];  // chr end pos
    hts_itr_t *iter = sam_itr_queryi(idx, tid, beg, end);
    if (iter == NULL) {   // region invalid or reference name not found
        cout << "ERROR" << endl;
        return -1;
    }
    
    // read bam or cram
    int ret;
    int processed_cnt=0;
    kstring_t aux={0, 0, NULL};
    while ((ret = sam_itr_next(fp, iter, b)) >= 0) {
        process_aln(fp, h, b);
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
    vector<int> sorted_chr;
    sort_chr_order(idx, h, sorted_chr);
    cout << "n = " << sorted_chr.size() << " chrs will be scanned." << endl;
    
    // close bam, keep bai
    sam_hdr_destroy(h);
    sam_close(fp);
    
    // threading
    ThreadPool pool(2);
    vector<future<int>> results;
    
    // per chr processing
    int i=0;
    for (int tid : sorted_chr) {
        results.emplace_back(
            pool.enqueue([=] {
                return extract_discordant_per_chr(f, idx, tid);
            })
        );
//        if (i >= 4) {
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
