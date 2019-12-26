#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


# load MEI info from RepBase
import utils
mes,all_clas=utils.load_me_classification(filenames.reshaped_rep)


def pairing(args, params, filenames):
    # load blast results
    L_poss={}
    R_poss={}
    for m in all_clas:
        L_poss[m]={}
        R_poss[m]={}
        for chr in args.main_chrs_set:
            L_poss[m][chr]={}
            R_poss[m][chr]={}
    files=(filenames.overhang_MEI, filenames.unmapped_MEI)
    for f in files:
        with open(f) as infile:
            for line in infile:
                ls=line.split()
                cposs=ls[0].split(';')[:-1]
                mposs=ls[1].split(';')
                mposs=[ mes[m] for m in mposs ]
                mposs=set(mposs)
                for c in cposs:
                    chr,tmp=c.split(':', 1)
                    start,tmp=tmp.split('-', 1)
                    end,breakpoint,readname=tmp.split('/', 2)
                    if breakpoint == 'L':
                        start=int(start)
                        for m in mposs:
                            if not start in L_poss[m][chr]:
                                L_poss[m][chr][start]=[]
                            L_poss[m][chr][start].append(c)
                    else:
                        end=int(end)
                        for m in mposs:
                            if not end in R_poss[m][chr]:
                                R_poss[m][chr][end]=[]
                            R_poss[m][chr][end].append(c)
    # load breakpoints from pA file
    L_pA,R_pA={},{}
    for chr in args.main_chrs_set:
        L_pA[chr]={}
        R_pA[chr]={}
    files=(filenames.overhang_pA, filenames.additional_pA)
    for f_pA in files:
        with open(f_pA) as infile:
            for line in infile:
                line=line.split()[0]
                chr,tmp=line.split(':', 1)
                start,tmp=tmp.split('-', 1)
                end,breakpoint,readname=tmp.split('/', 2)
                start,end=int(start),int(end)
                if breakpoint == 'L':
                    if not start in L_pA[chr]:
                        L_pA[chr][start]=[]
                    L_pA[chr][start].append(line)
                else:
                    if not end in R_pA[chr]:
                        R_pA[chr][end]=[]
                    R_pA[chr][end].append(line)
    # sweep, gather pA-only overhang and normal overhang
    outfile=open(filenames.breakpoint_pairs, 'w')
    L_added={}
    R_added={}
    for m in all_clas:
        out[m]=[]
        L_added[m]={}
        R_added[m]={}
        for chr in args.main_chrs_set:
            L_added[m][chr]={}
            R_added[m][chr]={}
    for m in all_clas:
        pA_MEI=True if m in args.rep_with_pA else False
        for chr in args.main_chrs_set:
            L_keys=L_poss[m][chr].keys()
            R_keys=R_poss[m][chr].keys()
            if (len(L_keys) >= 1) and (len(R_keys) >= 1):
                L_keys=sorted(list(L_keys))
                R_keys=sorted(list(R_keys))
                iter_n=0
                for l in L_keys:
                    for r in R_keys[iter_n:]:
                        if -params.max_breakpoint_gap < (l - r) < params.max_TSD_len:
                            if pA_MEI is True:
                                L_pA_num,R_pA_num=0,0
                                L_pA_id,R_pA_id='NA','NA'
                                if l in L_pA[chr]:
                                    L_pA_num += len(L_pA[chr][l])
                                    L_pA_id= ';'.join(L_pA[chr][l])
                                if r in R_pA[chr]:
                                    R_pA_num += len(R_pA[chr][r])
                                    R_pA_id= ';'.join(R_pA[chr][r])
                                if (L_pA_num >= 1) and (R_pA_num >= 1):
                                    L_pA_num,R_pA_num=0,0
                                    L_pA_id,R_pA_id='NA','NA'
                                L_count= len(L_poss[m][chr][l]) + L_pA_num
                                R_count= len(R_poss[m][chr][r]) + R_pA_num
                                L_id=';'.join(L_poss[m][chr][l])
                                R_id=';'.join(R_poss[m][chr][r])
                                L_pA_num,R_pA_num=str(L_pA_num),str(R_pA_num)
                            else:
                                L_pA_num,R_pA_num='NA','NA'
                                L_count= len(L_poss[m][chr][l])
                                R_count= len(R_poss[m][chr][r])
                                L_id=';'.join(L_poss[m][chr][l])
                                R_id=';'.join(R_poss[m][chr][r])
                                R_pA_id,L_pA_id='NA','NA'
                            if (L_count >= params.min_read_num_per_breakpoint_edge) and (R_count >= params.min_read_num_per_breakpoint_edge):
                                outfile.write(chr +'\t'+ str(min(l, r)) +'\t'+ str(max(l, r)) +'\t'+ str(r)+':'+str(R_count) +'\t'+ str(l)+':'+str(L_count) +'\t'+ R_pA_num +'\t'+ L_pA_num +'\t'+ m +'\t'+ R_id +'\t'+ L_id +'\t'+ R_pA_id +'\t'+ L_pA_id +'\n')
                                if l in L_pA[chr]:
                                    del(L_pA[chr][l])
                                if r in R_pA[chr]:
                                    del(R_pA[chr][r])
                        elif (l - r) >= params.max_breakpoint_gap:
                            iter_n += 1
                        elif (r - l) >= params.max_TSD_len:
                            break
    # antisense-direction, long pA
    for m in args.rep_with_pA:
        if m in all_clas:
            for chr in args.main_chrs_set:
                L_keys=L_poss[m][chr].keys()
                R_keys=R_pA[chr].keys()
                if (len(L_keys) >= 1) and (len(R_keys) >= 1):
                    L_keys=sorted(list(L_keys))
                    R_keys=sorted(list(R_keys))
                    iter_n=0
                    for l in L_keys:
                        for r in R_keys[iter_n:]:
                            if -params.max_breakpoint_gap < (l - r) < params.max_TSD_len:
                                L_count= len(L_poss[m][chr][l])
                                R_count= len(R_pA[chr][r])
                                L_id=';'.join(L_poss[m][chr][l])
                                R_pA_id=';'.join(R_pA[chr][r])
                                R_id,L_pA_id='NA','NA'
                                if (L_count >= params.min_read_num_per_breakpoint_edge) and (R_count >= params.min_read_num_per_breakpoint_edge):
                                    outfile.write(chr +'\t'+ str(min(l, r)) +'\t'+ str(max(l, r)) +'\t'+ str(r)+':'+str(R_count) +'\t'+ str(l)+':'+str(L_count) +'\t'+ str(R_count) +'\t'+ 'NA' +'\t'+ m +'\t'+ R_id +'\t'+ L_id +'\t'+ R_pA_id +'\t'+ L_pA_id +'\n')
                            elif (l - r) >= params.max_breakpoint_gap:
                                iter_n += 1
                            elif (r - l) >= params.max_TSD_len:
                                break
    # antisense-direction, long pA
    for m in args.rep_with_pA:
        if m in all_clas:
            for chr in args.main_chrs_set:
                L_keys=L_pA[chr].keys()
                R_keys=R_poss[m][chr].keys()
                if (len(L_keys) >= 1) and (len(R_keys) >= 1):
                    L_keys=sorted(list(L_keys))
                    R_keys=sorted(list(R_keys))
                    iter_n=0
                    for r in R_keys:
                        for l in L_keys[iter_n:]:
                            if -params.max_breakpoint_gap < (l - r) < params.max_TSD_len:
                                L_count= len(L_pA[chr][l])
                                R_count= len(R_poss[m][chr][r])
                                L_pA_id=';'.join(L_pA[chr][l])
                                R_id=';'.join(R_poss[m][chr][r])
                                L_id,R_pA_id='NA','NA'
                                if (L_count >= params.min_read_num_per_breakpoint_edge) and (R_count >= params.min_read_num_per_breakpoint_edge):
                                    outfile.write(chr +'\t'+ str(min(l, r)) +'\t'+ str(max(l, r)) +'\t'+ str(r)+':'+str(R_count) +'\t'+ str(l)+':'+str(L_count) +'\t'+ 'NA' +'\t'+ str(L_count) +'\t'+ m +'\t'+ R_id +'\t'+ L_id +'\t'+ R_pA_id +'\t'+ L_pA_id +'\n')
                            elif (r - l) >= params.max_breakpoint_gap:
                                iter_n += 1
                            elif (l - r) >= params.max_TSD_len:
                                break
    outfile.close()


