#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,gzip,shutil

def gzip_or_del(args, params, file):
    if args.keep is None:
        os.remove(file)
    else:
        with open(file, 'rt') as f_in:
            with gzip.open(file +'.gz', mode='wt', compresslevel=params.gzip_compresslevel) as f_out:
                shutil.copyfileobj(f_in, f_out)


def parse_fasta(path_to_file):
    tmp={}
    seq=''
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line and seq:
                tmp[header]=seq
                header=line.strip().replace(' ', '_')
                seq=''
            elif '>' in line and not seq:
                header=line.strip().replace(' ', '_')
            else:
                seq += line.strip()
        tmp[header]=seq
    return tmp

