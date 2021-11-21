#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os 
import pybedtools
from pybedtools import BedTool
from collections import Counter
from statistics import mean
from utils import load_me_classification, parse_fasta_transd
import log,traceback


def find_abs(args, params, filenames):
    log.logger.debug('started')
    try:
        pybedtools.set_tempdir(args.pybedtools_tmp)
        
        min_chimeric_num= round(args.cov * params.abs_min_chimeric_num_coeff)
        if min_chimeric_num < 2:
            min_chimeric_num=2
        mes,_=load_me_classification(filenames.reshaped_rep)
        
        # load chimeric reads
        d={}
        with open(filenames.abs_txt) as infile:
            for line in infile:
                ls=line.split()
                if ls[1] in args.main_chrs_set:
                    R= int(ls[6].split('-')[1])
                    L= int(ls[7].split('-')[0])
                    id='\t'.join(ls[1:4])
                    if not id in d:
                        d[id]=[[],[]]
                    d[id][0].append(R - L)
                    d[id][1].append(ls[0])

        # identify hotspots
        bed_high_cov=[]
        for id in d:
            if len(d[id][0]) >= min_chimeric_num:
                m=Counter(d[id][0]).most_common(1)[0][0]
                bed_high_cov.append('%s\t%s\t%d\n' % (id, ';'.join(d[id][1]), m))
        del(d)
        bed_high_cov=BedTool(''.join(bed_high_cov), from_string=True)

        # identify unfixed TEs
        bed_te=BedTool(filenames.repout_bed)
        bed_te_intersect=bed_te.intersect(bed_high_cov, wa=True, wb=True, nonamecheck=True)
        d={}
        for line in bed_te_intersect:
            ls=str(line).split()
            id='\t'.join(ls[5:8])
            if not id in d:
                d[id]=[]
            d[id].append([int(ls[1]), int(ls[2]), ls[3], ls[8], ls[4], ls[9]])
        pybedtools.cleanup()

        outfile_abs=open(filenames.abs_res, 'w')
        outfile_transd=open(filenames.transd_res, 'w')
        global abs_n
        abs_n=0
        transd_l_f,transd_l_r=[],[]
        bed_for_transd_f,bed_for_transd_r=[],[]
        for id in d:
            chr,start,end=id.split('\t')
            start,end=int(start),int(end)
            ss,es=[],[]
            for s,e,_,_,_,_ in d[id]:
                ss.append(abs(start - s))
                es.append(abs(end   - e))
            te_start=min(ss)
            te_end  =min(es)
            if (te_start <= params.breakpoint_annotation_gap) and (te_end <= params.breakpoint_annotation_gap):
                count=0
                pos_set=set()
                for s,e,_,_,_,_ in d[id]:
                    for pos in range(s, e):
                        pos_set.add(pos)
                for pos in pos_set:
                    if start <= pos < end:
                        count += 1
                te_ratio= count / (end - start)
                if te_ratio >= params.abs_len_to_te_ratio:
                    te_names=[]
                    non_ME_len=0
                    for s,e,n,r,_,t in d[id]:
                        if start <= s:
                            c=(e - s) if e <= end else (end - s)
                        else:
                            c=(e - start) if e <= end else (end - start)
                        if (c / (e-s)) >= params.len_te_for_abs_ratio:
                            te_names.append(n)
                        te_name=n.split(':')[0]
                        if not te_name in mes:
                            non_ME_len += c
                    if (len(te_names) >= 1) and ((non_ME_len / (end - start)) <= params.non_ME_len_ratio):
                        te_names=sorted(list(set(te_names)))
                        outfile_abs.write('%s\t%d\t%d\t%s\t%s\tTSD_len=%s\n' % (chr, start, end, ';'.join(te_names), r, t))
                        abs_n += 1
            elif (te_start <= params.breakpoint_annotation_gap) or (te_end <= params.breakpoint_annotation_gap):
                if te_start <= params.breakpoint_annotation_gap:
                    for s,e,n,r,strand,t in d[id]:
                        if (abs(start - s) == te_start) and (strand == '+') and (e < end):
                            te_name=n.split(':')[0]
                            if te_name in mes:
                                te_clas=mes[te_name]
                                if te_clas in args.rep_with_pA:
                                    transd_l_f.append([chr, start, end, n, r, t])
                                    bed_for_transd_f.append('%s\t%d\t%d\n' % (chr, end - params.transduction_pA_len, end))
                elif te_end <= params.breakpoint_annotation_gap:
                    for s,e,n,r,strand,t in d[id]:
                        if (abs(end - e) == te_end) and (strand == '-') and (start < s):
                            te_name=n.split(':')[0]
                            if te_name in mes:
                                te_clas=mes[te_name]
                                if te_clas in args.rep_with_pA:
                                    transd_l_r.append([chr, start, end, n, r, t])
                                    bed_for_transd_r.append('%s\t%d\t%d\n' % (chr, start, start + params.transduction_pA_len))
        # transd_f
        transd_threshold=params.transduction_pA_len * params.transduction_pA_ratio
        tail=BedTool(''.join(bed_for_transd_f), from_string=True)
        tail=tail.sequence(fi=args.fa)
        fa=parse_fasta_transd(tail.seqfn)
        for l,seq in zip(transd_l_f, fa):
            if seq.count('A') >= transd_threshold:
                chr,start,end,n,r,t=l
                outfile_transd.write('%s\t%d\t%d\t%s\t%s\tTSD_len=%s\n' % (chr, start, end, n, r, t))
        pybedtools.cleanup()
        # transd_r
        tail=BedTool(''.join(bed_for_transd_r), from_string=True)
        tail=tail.sequence(fi=args.fa)
        fa=parse_fasta_transd(tail.seqfn)
        for l,seq in zip(transd_l_r, fa):
            if seq.count('T') >= transd_threshold:
                chr,start,end,n,r,t=l
                outfile_transd.write('%s\t%d\t%d\t%s\t%s\tTSD_len=%s\n' % (chr, start, end, n, r, t))
        pybedtools.cleanup()
        outfile_abs.flush()
        outfile_transd.flush()
        os.fdatasync(outfile_abs.fileno())
        os.fdatasync(outfile_transd.fileno())
        outfile_abs.close()
        outfile_transd.close()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
