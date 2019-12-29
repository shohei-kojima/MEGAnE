#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


from utils import parse_fasta
from pybedtools import BedTool


def filter(args, params, filenames):
    nts=['A', 'T']
    total_read_threshold= round(args.cov * params.first_filter_total_read_num_coefficient)
    zero_hybrid_total_read_threshold= round((args.cov * params.first_filter_total_read_num_coefficient_for_zero_hybrid))

    def L1_filter(line, r_pos, l_pos):
        cand=True
        if pA_only is False:
            ls=line.split()
            bed=''
            if (r_pos - l_pos) >= params.L1_filter_min_TSD_len:
                bed += ls[0] +'\t'+ str(l_pos) +'\t'+ str(r_pos) +'\n'
                bed=BedTool(bed, from_string=True)
                fa=bed.sequence(fi=args.fa)
                fa=parse_fasta(fa.seqfn)
                for h in fa:
                    seq=fa[h].upper()
                    seqlen=len(seq)
                    total_AT=0
                    for nt in nts:
                        c=seq.count(nt)
                        if (100 * (c / seqlen)) >= params.L1_filter_A_or_T_perc:
                            cand=False
                        total_AT += c
                    if (100 * (total_AT / seqlen)) >= params.L1_filter_A_plus_T_perc:
                        cand=False
        if len(R_eval) >= 1:
            if not (min(R_eval) < params.L1_filter_eval_threshold):
                cand=False
        if len(L_eval) >= 1:
            if not (min(L_eval) < params.L1_filter_eval_threshold):
                cand=False
        return cand

    all=[]
    high=set()
    with open(filenames.bp_merged_all) as infile:
        for line in infile:
            ls=line.split('\t')
            r_pos,r_num=ls[3].split(':')
            l_pos,l_num=ls[4].split(':')
            r_num=0 if r_num == 'NA' else int(r_num)
            l_num=0 if l_num == 'NA' else int(l_num)
            r_pos=0 if r_pos == 'NA' else int(r_pos)
            l_pos=0 if l_pos == 'NA' else int(l_pos)
            total_read_count= r_num + l_num
            R_eval,L_eval=[],[]
            vs=ls[8].split(';')
            vs=[ float(v) for v in vs if not (v == 'NA') and not (v == '') ]
            if len(vs) >= 1:
                for v in vs:
                    R_eval.append(v)
            vs=ls[9].split(';')
            vs=[ float(v) for v in vs if not (v == 'NA') and not (v == '') ]
            if len(vs) >= 1:
                for v in vs:
                    L_eval.append(v)
            retain_count=False
            if total_read_count >= total_read_threshold:
                retain_count=True
            retain_eval=False
            for l in [R_eval, L_eval]:
                if len(l) >= 1:
                    vmin=min(l)
                    if vmin <= params.first_filter_eval_threshold:
                        retain_eval=True
            # first filter
            if (retain_count is True) and (retain_eval is True):
                if (int(ls[12]) + int(ls[13])) >= params.first_filter_total_hybrid_read_num:
                    all.append(line)
                elif total_read_count >= zero_hybrid_total_read_threshold:
                    all.append(line)
            # second filter, 90% accuracy
            if (retain_count is True) and (retain_eval is True):
                if (int(ls[12]) >= params.second_filter_hybrid_read_num) and (int(ls[13]) >= params.second_filter_hybrid_read_num):
                    if 'L1' in line:
                        L1_judge=L1_filter(line, r_pos, l_pos)
                        if L1_judge is True:
                            high.add(line)
                    else:
                        high.add(line)
                elif (len(R_eval) >= 1) and (len(L_eval) >= 1):
                    if (min(R_eval) < params.second_filter_eval_threshold_for_few_hybrid) and (min(L_eval) < params.second_filter_eval_threshold_for_few_hybrid):
                        if 'L1' in line:
                            L1_judge=L1_filter(line, r_pos, l_pos)
                            if L1_judge is True:
                                high.add(line)
                        else:
                            high.add(line)

    with open(filenames.bp_merged_filt, 'w') as outfile:
        for line in all:
            if line in high:
                line=line.strip()
                outfile.write(line +'\thigh\n')
            else:
                line=line.strip()
                outfile.write(line +'\tlow\n')


def grouping(args, filenames):
    def merge(one_l, other_l):
        orig_len=len(one_l)
        tmp=[]
        tmp.extend(one_l)
        unmerged=[]
        for a in other_l:
            merged=False
            for line in one_l:
                if line in a:
                    tmp.extend(a)
                    merged=True
                    break
            if merged is False:
                unmerged.append(a)
        tmp=list(set(tmp))
        tmp=sorted(tmp)
        after_len=len(tmp)
        if after_len > orig_len:
            merged=True
        else:
            merged=False
        return merged, tmp, unmerged

    # group non-singletons
    L,R={},{}
    for chr in args.main_chrs_set:
        L[chr]={}
        R[chr]={}
    with open(filenames.bp_merged_filt) as infile:
        for line in infile:
            line=line.strip()
            ls=line.split()
            r=ls[3].split(':')[0]
            l=ls[4].split(':')[0]
            L[ls[0]][l]=line
            R[ls[0]][r]=line

    all_0=[]
    with open(filenames.overhang_MEI) as infile:
        for line in infile:
            ls=line.split()
            poss=ls[0].split(';')[:-1]
            tmp_s=[]
            for p in poss:
                chr,tmp=p.split(':', 1)
                start,tmp=tmp.split('-', 1)
                end,dir,_=tmp.split('/', 2)
                if dir == 'L':
                    if start in L[chr]:
                        tmp_s.append(L[chr][start])
                else:
                    if end in R[chr]:
                        tmp_s.append(R[chr][end])
            all_0.append(tmp_s)

    all_1=[]
    for a in all_0:
        if not a in all_1:
            all_1.append(a)
    del(all_0)

    unmerged=all_1
    unmerged_len=len(unmerged)
    final=[]
    while unmerged_len >= 1:
        b=True
        l=unmerged[0]
        unmerged=unmerged[1:]
        while b is True:
            b,l,unmerged=merge(l, unmerged)
        final.append(l)
        unmerged_len=len(unmerged)

    singletons,multis={},{}
    me_clas=set()
    n=0
    for f in final:
        if len(f) == 1:
            ls=f[0].split()
            if not ls[7] in singletons:
                singletons[ls[7]]=''
                me_clas.add(ls[7])
            singletons[ls[7]] += f[0] +'\tsingleton\n'
        else:
            for line in f:
                ls=line.split()
                if not ls[7] in singletons:
                    multis[ls[7]]=''
                    me_clas.add(ls[7])
                multis[ls[7]] += line +'\tgroup%d\n' % n
            n += 1
    me_clas=sorted(list(me_clas))

    with open(filenames.bp_merged_group, 'w') as outfile:
        for m in me_clas:
            if m in singletons:
                bed=BedTool(singletons[m], from_string=True).sort()
                for line in bed:
                    outfile.write(str(line))
        for m in me_clas:
            if m in multis:
                bed=BedTool(multis[m], from_string=True).sort()
                for line in bed:
                    multis_sorted += str(line)

