#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

def find(params, blast_res, q_path, outfpath):
    min_A_count_per_bin= params.pA_scan_bin - params.max_non_pA_count
    overhang_len_for_full_analysis= params.pA_scan_bin + params.scan_loop_from_edge - 1

    def judge_pA(seq):  # more than 10 'A' in 12-nt
        pA_judge=False
        seqlen=len(seq)
        if seqlen >= params.pA_scan_bin:
            loop=params.scan_loop_from_edge
            if seqlen < overhang_len_for_full_analysis:
                loop= seqlen - params.pA_scan_bin + 1
            for n in range(loop):
                s=seq[-n - (1 + params.pA_scan_bin) : -n - 1]
                pA_c=s.count('A')
                if pA_c >= min_A_count_per_bin:
                    pA_judge=True
                    break
        return pA_judge

    def judge_pT(seq):  # more than 10 'A' in 12-nt
        pA_judge=False
        seqlen=len(seq)
        if seqlen >= params.pA_scan_bin:
            loop=params.scan_loop_from_edge
            if seqlen < overhang_len_for_full_analysis:
                loop= seqlen - params.pA_scan_bin + 1
            for n in range(loop):
                s=seq[n : n + params.pA_scan_bin]
                pA_c=s.count('T')
                if pA_c >= min_A_count_per_bin:
                    pA_judge=True
                    break
        return pA_judge

    # load blast results
    hits=set()
    with open(blast_res) as infile:
        for line in infile:
            ls=line.split()
            hits.add(ls[0])

    # judge overhangs
    with open(outfpath, 'w') as outfile:
        with open(q_path) as infile:
            for line in infile:
                h=line.strip().replace('>', '')
                seq=next(infile)
                seq=seq.strip().upper()
                if not h in hits:
                    if '/L/' in h:
                        pA_b=judge_pA(seq)
                        if pA_b is True:
                            seqlen=str(len(seq))
                            for i in h.split(';'):
                                if '/L/' in i:
                                    outfile.write(i +'\t'+ seqlen +'\n')
                    if '/R/' in h:
                        pT_b=judge_pT(seq)
                        if pT_b is True:
                            seqlen=str(len(seq))
                            for i in h.split(';'):
                                if '/R/' in i:
                                    outfile.write(i +'\t'+ seqlen +'\n')
