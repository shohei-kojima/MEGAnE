#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.

Prerequisite:
    Python3.7 or later.

Usage:
    python calc_mapping_depth.py -i Dfam.fa -o Dfam.rep

Usage details:
    Please see the pdf of MEGAnE instruction for further details.

Show help message:
    python Dfam_fa_to_rep.py -h
'''

import os,sys,argparse

parser=argparse.ArgumentParser(description='')
parser.add_argument('-i', metavar='str', type=str, help='Required. Specify fasta file containing mobile elements. This must be generated from Dfam.')
parser.add_argument('-o', metavar='str', type=str, help='Optional. Specify output file name. Otherwise, output result as stdout.')
parser.add_argument('-do_not_overwrite', help='Optional. Specify if you do NOT want to overwrite previous results.', action='store_true')
args=parser.parse_args()

# initial check
if args.i is None:
    print('Error: You need to specify input file with "-i" option.', file=sys.stderr)
    exit(1)
else:
    if os.path.exists(args.i) is not True:
        print('Error: Input file does not exist.', file=sys.stderr)
        exit(1)

if args.o is not None:
    if os.path.exists(args.o) is True:
        if args.do_not_overwrite is True:
            print('Error: Output file already exists. Please remove "-do_not_overwrite" option if you overwrite the existing file.', file=sys.stderr)
            exit(1)
        else:
            print('Info: Output file already existed, but overwrote it according to the absence of "-do_not_overwrite" option.')


# main
out=[]
header_n=0
clas_set=set()
with open(args.i) as infile:
    for line in infile:
        if '>' in line:
            if not '#' in line:
                print('Error: fasta header of your input does not look like the usual fasta file conveted from Dfam.', file=sys.stderr)
                exit(1)
            ls=line.split()
            if len(ls) <= 1:
                print('Error: fasta header of your input does not look like the usual fasta file conveted from Dfam.', file=sys.stderr)
                exit(1)
            if not '@' in ls[1]:
                print('Error: fasta header of your input does not look like the usual fasta file conveted from Dfam.', file=sys.stderr)
                exit(1)
            species=ls[1].replace('@', '')
            me=ls[0].replace('>', '').split('#')[0]
            clas=ls[0].split('#')[1]
            clas_set.add(clas)
            out.append('>%s\t%s\t%s\n' % (me, clas, species))
            header_n += 1
        else:
            out.append(line)


if args.o is not None:
    with open(args.o, 'w') as outfile:
        outfile.write(''.join(out))
else:
    print(''.join(out), end='')

print('%d entry(s) processed.' % header_n, file=sys.stderr)

clas=sorted(list(clas_set))
print('%d subclass found:\n\n%s\n\nWhen running MEGAnE with this output, please specify 1) subclass names of non-ME subclasses with the "-repremove" flag and 2) subclass names of MEs with polyA tail with the "-pA_ME" flag.' % (len(clas), '\n'.join(clas)), file=sys.stderr)
