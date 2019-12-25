#!/usr/bin/env python

"""
2019/12/05 15:56:32 +0900

# batch_191205_155632.py
# generated at /home/kooojiii/results/misc/1000_hgp/191203_1
# usage: python %prog
# python3.7

"""

# identify junctions by each TE class

import os,sys,glob

# params
minimum_read_num_per_breakpoint=1
minimum_read_num_per_single_breakpoint=3

hu_chrs=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
MEI_pA=['SINE', 'SINE1/7SL', 'SINE2/tRNA', 'SINE3/5S', 'L1', 'L2', 'RTE', 'RTEX', 'Non-LTR Retrotransposon']


# load MEI info from RepBase
f='/home/kooojiii/results/misc/1000_hgp/lib/humrepsub.fa'
meis={}
all_clas=set()
with open(f) as infile:
    for line in infile:
        if '>' in line:
            ls=line.strip().replace('>', '').split('\t')
            if len(ls) >= 2:
                clas=ls[1]
                if ls[0] == 'SVA2':
                    clas='SINE'
            elif ('Alu' in ls[0]) or ('FLA' in ls[0]):
                clas='SINE1/7SL'
            elif 'LTR26' in ls[0]:
                clas='ERV1'
            else:
                print('read error by kojima: %s' % f)
            clas=clas.replace(' ', '_')
            meis[ls[0]]=clas
            all_clas.add(clas)
        pass



# load blast results
L_poss={}
R_poss={}
for m in all_clas:
    L_poss[m]={}
    R_poss[m]={}
    for chr in hu_chrs:
        L_poss[m][chr]={}
        R_poss[m][chr]={}

files=['overhang_MEI.txt', 'overhang_MEI_unmapped.txt']
for f in files:
    with open(f) as infile:
        for line in infile:
            ls=line.split()
            cposs=ls[0].split(';')[:-1]
            mposs=ls[1].split(';')
            mposs=[ meis[m] for m in mposs ]
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
            pass


# load breakpoints from pA file
L_pA,R_pA={},{}
for chr in hu_chrs:
    L_pA[chr]={}
    R_pA[chr]={}
files=['overhang_pA.txt', 'overhang_pA_additional.txt']
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
            pass

# sweep, gather pA-only overhang and normal overhang
out={}
L_added={}
R_added={}
for m in all_clas:
    out[m]=[]
    L_added[m]={}
    R_added[m]={}
    for chr in hu_chrs:
        L_added[m][chr]={}
        R_added[m][chr]={}
for m in all_clas:
    pA_MEI=False
    if m in MEI_pA:
        pA_MEI=True
    for chr in hu_chrs:
        L_keys=L_poss[m][chr].keys()
        R_keys=R_poss[m][chr].keys()
        if (len(L_keys) >= 1) and (len(R_keys) >= 1):
            L_keys=sorted(list(L_keys))
            R_keys=sorted(list(R_keys))
            iter_n=0
            for l in L_keys:
                for r in R_keys[iter_n:]:
                    if -50 < (l - r) < 50:
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
                        if (L_count >= minimum_read_num_per_breakpoint) and (R_count >= minimum_read_num_per_breakpoint):
                            L_added[m][chr][l]=''
                            R_added[m][chr][r]=''
                            out[m].append(chr +'\t'+ str(min(l, r)) +'\t'+ str(max(l, r)) +'\t'+ str(r)+':'+str(R_count) +'\t'+ str(l)+':'+str(L_count) +'\t'+ R_pA_num +'\t'+ L_pA_num +'\t'+ m +'\t'+ R_id +'\t'+ L_id +'\t'+ R_pA_id +'\t'+ L_pA_id)
                            if l in L_pA[chr]:
                                del(L_pA[chr][l])
                            if r in R_pA[chr]:
                                del(R_pA[chr][r])
                    elif (l - r) >= 50:
                        iter_n += 1
                    elif (r - l) >= 50:
                        break
                    pass


# antisense-direction, long pA
for m in MEI_pA:
    if m in all_clas:
        for chr in hu_chrs:
            L_keys=L_poss[m][chr].keys()
            R_keys=R_pA[chr].keys()
            if (len(L_keys) >= 1) and (len(R_keys) >= 1):
                L_keys=sorted(list(L_keys))
                R_keys=sorted(list(R_keys))
                iter_n=0
                for l in L_keys:
                    for r in R_keys[iter_n:]:
                        if -50 < (l - r) < 50:
                            L_count= len(L_poss[m][chr][l])
                            R_count= len(R_pA[chr][r])
                            L_id=';'.join(L_poss[m][chr][l])
                            R_pA_id=';'.join(R_pA[chr][r])
                            R_id,L_pA_id='NA','NA'
                            if (L_count >= minimum_read_num_per_breakpoint) and (R_count >= minimum_read_num_per_breakpoint):
                                L_added[m][chr][l]=''
                                out[m].append(chr +'\t'+ str(min(l, r)) +'\t'+ str(max(l, r)) +'\t'+ str(r)+':'+str(R_count) +'\t'+ str(l)+':'+str(L_count) +'\t'+ str(R_count) +'\t'+ 'NA' +'\t'+ m +'\t'+ R_id +'\t'+ L_id +'\t'+ R_pA_id +'\t'+ L_pA_id)
                        elif (l - r) >= 50:
                            iter_n += 1
                        elif (r - l) >= 50:
                            break
                        pass


# antisense-direction, long pA
for m in MEI_pA:
    if m in all_clas:
        for chr in hu_chrs:
            L_keys=L_pA[chr].keys()
            R_keys=R_poss[m][chr].keys()
            if (len(L_keys) >= 1) and (len(R_keys) >= 1):
                L_keys=sorted(list(L_keys))
                R_keys=sorted(list(R_keys))
                iter_n=0
                for r in R_keys:
                    for l in L_keys[iter_n:]:
                        if -50 < (l - r) < 50:
                            L_count= len(L_pA[chr][l])
                            R_count= len(R_poss[m][chr][r])
                            L_pA_id=';'.join(L_pA[chr][l])
                            R_id=';'.join(R_poss[m][chr][r])
                            L_id,R_pA_id='NA','NA'
                            if (L_count >= minimum_read_num_per_breakpoint) and (R_count >= minimum_read_num_per_breakpoint):
                                R_added[m][chr][r]=''
                                out[m].append(chr +'\t'+ str(min(l, r)) +'\t'+ str(max(l, r)) +'\t'+ str(r)+':'+str(R_count) +'\t'+ str(l)+':'+str(L_count) +'\t'+ 'NA' +'\t'+ str(L_count) +'\t'+ m +'\t'+ R_id +'\t'+ L_id +'\t'+ R_pA_id +'\t'+ L_pA_id)
                        elif (r - l) >= 50:
                            iter_n += 1
                        elif (l - r) >= 50:
                            break
                    pass

# output paired breakpoints
tmp=[]
for m in out:
    for i in out[m]:
        tmp.append(i)
print(len(tmp))

with open('junction_breakpoint_pairs.txt', 'w') as outfile:
    outfile.write('\n'.join(tmp) +'\n')


# list up single breakpoints
out={}
for m in all_clas:
    out[m]=[]
for m in all_clas:
    for chr in hu_chrs:
        for l in L_poss[m][chr]:
            L_count= len(L_poss[m][chr][l])
            if L_count >= minimum_read_num_per_single_breakpoint:
                if not l in L_added[m][chr]:
                    L_id=';'.join(L_poss[m][chr][l])
                    out[m].append(chr +'\t'+ str(l) +'\t'+ str(l) +'\t'+ 'NA:NA' +'\t'+ str(l)+':'+str(L_count) +'\t'+ 'NA' +'\t'+ 'NA' +'\t'+ m +'\t'+ 'NA' +'\t'+ L_id +'\t'+ 'NA' +'\t'+ 'NA')
        for r in R_poss[m][chr]:
            R_count= len(R_poss[m][chr][r])
            if R_count >= minimum_read_num_per_single_breakpoint:
                if not r in R_added[m][chr]:
                    R_id=';'.join(R_poss[m][chr][r])
                    out[m].append(chr +'\t'+ str(r) +'\t'+ str(r) +'\t'+ str(r)+':'+str(R_count) +'\t'+ 'NA:NA' +'\t'+ 'NA' +'\t'+ 'NA' +'\t'+ m +'\t'+ R_id +'\t'+ 'NA' +'\t'+ 'NA' +'\t'+ 'NA')


tmp=[]
for m in out:
    for i in out[m]:
        tmp.append(i)
print(len(tmp))

with open('junction_breakpoint_single.txt', 'w') as outfile:
    outfile.write('\n'.join(tmp) +'\n')
