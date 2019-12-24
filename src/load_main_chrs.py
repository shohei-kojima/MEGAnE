#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

def load(f):
    chrs=[]
    with open(f) as infile:
        for line in infile:
            chr=line.strip()
            if not chr == '':
                chrs.append(chr)
    return chrs
