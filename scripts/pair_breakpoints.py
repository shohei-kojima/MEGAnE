#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os
import utils
from pybedtools import BedTool
import log,traceback


# list up simple repeat in the ref genome
def load_ref_simple_rep_as_bed(filenames):
    mes,all_clas=utils.load_me_classification(filenames.reshaped_rep)
    simple=''
    with open(filenames.repout_bed) as infile:
        for line in infile:
            ls=line.split()
            name,clas=ls[3].split(':')
            if clas == 'Simple_repeat':
                simple += line
            elif name in mes:
                simple += line
    simple=BedTool(simple, from_string=True)
    return simple


def pairing(args, params, filenames):
    log.logger.debug('started')
    try:
        mes,all_clas=utils.load_me_classification(filenames.reshaped_rep)
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
        outfile.flush()
        os.fdatasync(outfile.fileno())
        outfile.close()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def add_TE_subclass(args, filenames, infpath, outfpath):
    log.logger.debug('started')
    try:
        mes,all_clas=utils.load_me_classification(filenames.reshaped_rep)
        # load bed file
        L,R={},{}
        for chr in args.main_chrs_set:
            L[chr]={}
            R[chr]={}
        with open(infpath) as infile:
            for line in infile:
                ls=line.split()
                r=ls[3].split(':')[0]
                l=ls[4].split(':')[0]
                if not l in L[ls[0]]:
                    L[ls[0]][l]={}
                L[ls[0]][l][ls[7]]=[[],[]]
                if not r in R[ls[0]]:
                    R[ls[0]][r]={}
                R[ls[0]][r][ls[7]]=[[],[]]

        # load e-values
        files=(filenames.overhang_MEI, filenames.unmapped_MEI)
        for f in files:
            with open(f) as infile:
                for line in infile:
                    ls=line.split()
                    poss=ls[0].split(';')
                    meis=ls[1].split(';')
                    meis_d={}
                    for m in meis:
                        c=mes[m]
                        if not c in meis_d:
                            meis_d[c]=''
                        meis_d[c] += m +';'
                    for p in poss[:-1]:
                        chr,tmp=p.split(':', 1)
                        start,tmp=tmp.split('-', 1)
                        end,dir,_=tmp.split('/', 2)
                        if dir == 'L':
                            pos=start
                            if pos in L[chr]:
                                for c in meis_d:
                                    if c in L[chr][pos]:
                                        L[chr][pos][c][0].append(float(ls[2]))
                                        L[chr][pos][c][1].append(meis_d[c])
                        else:
                            pos=end
                            if pos in R[chr]:
                                for c in meis_d:
                                    if c in R[chr][pos]:
                                        R[chr][pos][c][0].append(float(ls[2]))
                                        R[chr][pos][c][1].append(meis_d[c])

        # add TE names to the existing file
        with open(outfpath, 'w') as outfile:
            with open(infpath) as infile:
                for line in infile:
                    line=line.strip()
                    ls=line.split()
                    r=ls[3].split(':')[0]
                    l=ls[4].split(':')[0]
                    if len(L[ls[0]][l][ls[7]][0]) >= 1:
                        l_min=min(L[ls[0]][l][ls[7]][0])
                    else:
                        l_min=100
                    if len(R[ls[0]][r][ls[7]][0]) >= 1:
                        r_min=min(R[ls[0]][r][ls[7]][0])
                    else:
                        r_min=100
                    eval=min(l_min,r_min)
                    vsl,vsr=[],[]
                    msl,msr=[],[]
                    if len(L[ls[0]][l][ls[7]][0]) >= 1:
                        for v,m in zip(L[ls[0]][l][ls[7]][0], L[ls[0]][l][ls[7]][1]):
                            vsl.append(v)
                            msl.append(m)
                        vsl=[ str(v) for v in vsl ]
                    else:
                        vsl.append('NA')
                        msl.append('NA')
                    if len(R[ls[0]][r][ls[7]][0]) >= 1:
                        for v,m in zip(R[ls[0]][r][ls[7]][0], R[ls[0]][r][ls[7]][1]):
                            vsr.append(v)
                            msr.append(m)
                        vsr=[ str(v) for v in vsr ]
                    else:
                        vsr.append('NA')
                        msr.append('NA')
                    outfile.write('\t'.join(ls[:8]) +'\t'+ ';'.join(vsr) +'\t'+ ';'.join(vsl) +'\t'+ ';'.join(msr) +'\t'+ ';'.join(msl) +'\n')
            outfile.flush()
            os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def remove_cand_inside_TE(args, params, filenames):
    log.logger.debug('started')
    try:
        # remove breakpoints in simple repeats
        simple=load_ref_simple_rep_as_bed(filenames)
        breakpoints=set()
        with open(filenames.breakpoint_info) as infile:
            for line in infile:
                ls=line.split()
                tmp=set()
                for l in ls[10:12]:
                    for i in l.split(';'):
                        if not (i == '') and not (i == 'NA'):
                            tmp.add(i)
                tmp=sorted(list(tmp))
                tmp=';'.join(tmp)
                breakpoints.add('\t'.join(ls[:3]) +'\t'+ tmp)
        breakpoints=list(breakpoints)
        breakpoints='\n'.join(breakpoints) +'\n'
        breakpoints=BedTool(breakpoints, from_string=True).sort()
        breakpoints=breakpoints.intersect(simple, v=True)
        del(simple)

        # retrieve candidate hotspots overlaps with reference TEs
        ref_genome_tes=BedTool(filenames.repout_bed)
        if not params.ref_TE_slop_len == 0:
            ref_genome_tes=ref_genome_tes.slop(g=args.fai, b=slop_len)
        breakpoints_intersect=breakpoints.intersect(ref_genome_tes, wa=True, wb=True)
        del(ref_genome_tes)

        # load similar TEs
        similar_te_d={}
        with open(filenames.similar_rep_list) as infile:
            for line in infile:
                ls=line.split()
                similar_te_d[ls[0]]=ls[1].split(';')

        # remove breakpoints overlaps with similar TEs
        remove=set()
        for line in breakpoints_intersect:
            line=str(line)
            ls=line.strip().split('\t')
            ref_te,ref_te_class=ls[7].split(':')
            if not ref_te_class == 'Simple_repeat':
                hit_te=set(ls[3].split(';'))
                tmp=[]
                for m in hit_te:
                    tmp.extend(similar_te_d[m])
                tmp=set(tmp)
                if ref_te in tmp:
                    remove.add('\t'.join(ls[:3]))

        retain=set()
        for line in breakpoints:
            line=str(line)
            ls=line.split()
            id='\t'.join(ls[:3])
            if not id in remove:
                retain.add(id)

        with open(filenames.breakpoint_clean, 'w') as outfile:
            with open(filenames.breakpoint_info) as infile:
                for line in infile:
                    ls=line.split()
                    id='\t'.join(ls[:3])
                    if id in retain:
                        outfile.write(line)
            outfile.flush()
            os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
