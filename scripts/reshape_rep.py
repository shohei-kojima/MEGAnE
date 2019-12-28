#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os
from utils import parse_fasta


def reshape(args, filenames):
    fa=parse_fasta(args.rep)
    fa_keep={}
    for header in fa:
        hs=header.split('\t')
        if len(hs) == 3:
            if not hs[1] in args.rep_headers_to_be_removed:
                fa_keep[header]=fa[header]
    with open(filenames.reshaped_rep, 'w') as outfile:
        for header in fa_keep:
            outfile.write(header +'\n'+ fa_keep[header] +'\n')


def slide_rep_file(args, params, filenames):
    
    def bin_slide(string, bin, interval):
        l=[]
        if len(string) >= bin:
            for i in range(0, len(string)-bin+1, interval):
                l.append(string[i:i+bin])
        else:
            l.append(string)
        return l

    fa=parse_fasta(filenames.reshaped_rep)
    with open(filenames.rep_slide_file, 'w') as outfile:
        for h in fa:
            seqs=bin_slide(fa[h], args.readlen, params.repbase_seq_slide_bin)
            tmp=[]
            n=0
            for seq in seqs:
                tmp.append(h.strip() +'/'+ str(n))
                tmp.append(seq)
                n += 1
            outfile.write('\n'.join(tmp) +'\n')


def parse_slide_rep_blastn_res(args, filenames):
    d={}
    with open(filenames.blast0_res) as infile:
        for line in infile:
            ls=line.split()
            id=ls[0].split('/')[0]
            if not id in d:
                d[id]=set()
            d[id].add(ls[1])
    for id in d:
        d[id]=sorted(list(d[id]))
    with open(filenames.similar_rep_list, 'w') as outfile:
        for id in d:
            outfile.write(id +'\t'+ ';'.join(d[id]) +'\n')


def reshape_repout_to_bed(args, filenames):
    with open(filenames.repout_bed, 'w') as outfile:
        with open(args.repout) as infile:
            for _ in range(3):
                next(infile)
            for line in infile:
                ls=line.split()
                start= int(ls[5]) - 1  # 0-based
                out= ls[4] +'\t'+ str(start) +'\t'+ ls[6] +'\t'+ ls[9]+':'+ls[10] +'\n'
                outfile.write(out)

