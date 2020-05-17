#!/usr/bin/env python

""'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.

Before using this, you need to run samtools coverage to calculate depth of each chromosome.
$ samtools coverage --reference ref.fa in.cram > chr_coverage.txt

Usage: python calc_mapping_depth.py chr_coverage.txt
This outputs 'chr_coverage_genome_depth.txt' which contains genome coverage.
'''

import os,sys

f=sys.argv[1]

genome_len=0
cov_dep=0
with open(f) as infile:
    for line in infile:
        if not line[0] == '#':
            ls=line.strip().split('\t')
            genome_len += int(ls[4])
            cov_dep += int(ls[4]) * float(ls[6])
mean_dep= cov_dep / genome_len
mean_dep= round(mean_dep)

base,_=os.path.splitext(f)
with open('%s_genome_depth.txt' % base, 'w') as outfile:
    outfile.write('%d\t%s\n' % (mean_dep, os.path.abspath(f)))
