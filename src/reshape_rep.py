#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os
from utils import parse_fasta

def reshape(args, outfpath):
    remove_headers={'HSFAU', 'KER', 'Pseudogene', 'SAT', 'Satellite', 'satellite', 'snRNA'}
    fa=parse_fasta(args.rep)
    fa_keep={}
    for header in fa:
        hs=header.split('\t')
        if len(hs) == 3:
            if not hs[1] in remove_headers:
                fa_keep[header]=fa[header]
    with open(outfpath, 'w') as outfile:
        for header in fa_keep:
            outfile.write(header +'\n'+ fa_keep[header] +'\n')


def slide_rep_file(args, params, reshaped_rep, outfpath):
    
    def bin_slide(string, bin, interval):
        l=[]
        if len(string) >= bin:
            for i in range(0, len(string)-bin+1, interval):
                l.append(string[i:i+bin])
        else:
            l.append(string)
        return l

    fa=parse_fasta(reshaped_rep)
    with open(outfpath, 'w') as outfile:
        for h in fa:
            seqs=bin_slide(fa[h], args.readlen, params.repbase_seq_slide_bin)
            tmp=[]
            n=0
            for seq in seqs:
                tmp.append(h.strip() +'/'+ str(n))
                tmp.append(seq)
                n += 1
            outfile.write('\n'.join(tmp) +'\n')


def parse_slide_rep_blastn_res(args, blast_res, outfpath):
    d={}
    with open(blast_res) as infile:
        for line in infile:
            ls=line.split()
            id=ls[0].split('/')[0]
            if not id in d:
                d[id]=set()
            d[id].add(ls[1])
    for id in d:
        d[id]=sorted(list(d[id]))
    with open(outfpath, 'w') as outfile:
        for id in d:
            outfile.write(id +'\t'+ ';'.join(d[id]) +'\n')

