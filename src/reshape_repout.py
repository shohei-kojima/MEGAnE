#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os

def reshape(f, outdir):
    outpath=os.path.join(outdir, 'reshaped.bed')
    with open(outpath, 'w') as outfile:
        with open(f) as infile:
            for _ in range(3):
                next(infile)
            for line in infile:
                ls=line.split()
                start= int(ls[5]) - 1  # 0-based
                if ls[8] == '+':
                    strand='+'
                else:
                    strand='-'
                out= ls[4] +'\t'+ str(start) +'\t'+ ls[6] +'\t'+ ls[9]+':'+ls[10] +'\t'+ strand +'\n'
                outfile.write(out)


