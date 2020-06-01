#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.

Prerequisite:
Python3 or later.

Usage: python calc_mapping_depth.py -i samtools_out.txt -o out.txt

Exmaple:
samtools coverage --reference ref.fa in.cram > in.txt
python calc_mapping_depth.py -i in.txt -chr /path/to/prog/lib/human_autosomes_ucsc_style.txt -o out.txt

Usage details:
Before using this script, you need to run 'samtools coverage' to calculate depth of each chromosome. You need to use samtools 1.10 or later to use 'coverage' function. The output from 'samtools coverage' contains mean depth of each chromosome. This program takes the output file from 'samtools coverage' to calculate mean depth of autosomes. You can specify the output file with '-i' option. You can specify names of autosomes by specifying a file containing names of autosomes with -chr option.

Show help message:
python calc_mapping_depth.py -h
'''

import os,sys,argparse

parser=argparse.ArgumentParser(description='')
parser.add_argument('-i', metavar='str', type=str, help='Required. Specify output file from samtools coverage.')
parser.add_argument('-chr', metavar='str', type=str, help='Required. Specify a file containing autosomes. Default=/path/to/prog/lib/human_autosomes_ucsc_style.txt')
parser.add_argument('-o', metavar='str', type=str, help='Optional. Specify output file name. Otherwise, output result as stdout.')
parser.add_argument('-overwrite', help='Optional. Specify if you overwrite previous results.', action='store_true')
args=parser.parse_args()


# initial check
if args.i is None:
    print('Error: You need to specify input file with "-i" option.')
    exit(1)
else:
    if os.path.exists(args.i) is not True:
        print('Error: Input file does not exist.')
        exit(1)

if args.chr is None:
    base=os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    args.chr=os.path.join(base, '../lib/human_autosomes_ucsc_style.txt')
else:
    if os.path.exists(args.chr) is not True:
        print('Error: Chromosome file do not exists.')
        exit(1)
    
if args.o is not None:
    if os.path.exists(args.o) is True:
        if args.overwrite is False:
            print('Error: Output file already exists. Please specify "-overwrite" option if you overwrite the existing file.')
            exit(1)
        else:
            print('Info: Output file already existed, but overwrote it according to "-overwrite" option.')


# main
chrs=set()
with open(args.chr) as infile:
    for line in infile:
        chrs.add(line.strip())

genome_len=0
cov_dep=0
chrs_used=[]
with open(args.i) as infile:
    for line in infile:
        if not line[0] == '#':
            ls=line.strip().split('\t')
            if ls[0] in chrs:
                genome_len += int(ls[4])
                cov_dep += int(ls[4]) * float(ls[6])
                chrs_used.append(ls[0])
mean_dep= cov_dep / genome_len
mean_dep= round(mean_dep)

if args.o is not None:
    with open(args.o, 'w') as outfile:
        outfile.write('%d\t%s\t%s\n' % (mean_dep, os.path.abspath(args.i), ';'.join(chrs_used)))
else:
    print('%d\t%s\t%s' % (mean_dep, os.path.abspath(args.i), ';'.join(chrs_used)))
