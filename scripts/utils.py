#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,gzip,shutil


class empclass:
    pass


def gzip_or_del(args, params, file):
    if args.keep is True:
        with open(file, 'rt') as f_in:
            with gzip.open(file +'.gz', 'wt', compresslevel=params.gzip_compresslevel) as f_out:
                shutil.copyfileobj(f_in, f_out)
    os.remove(file)


def gzip_file(params, file):
    with open(file, 'rt') as f_in:
        with gzip.open(file +'.gz', 'wt', compresslevel=params.gzip_compresslevel) as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file)


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


def load_me_classification(path_to_file):
    mes,all_clas={},set()
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line:
                ls=line.strip().replace('>', '').split('\t')
                if len(ls) >= 2:
                    clas=ls[1]
                    if ls[0] == 'SVA2':  # exist in humrep.ref of RepBase24.01
                        clas='SINE'
                clas=clas.replace(' ', '_')
                mes[ls[0]]=clas
                all_clas.add(clas)
    return mes, all_clas

