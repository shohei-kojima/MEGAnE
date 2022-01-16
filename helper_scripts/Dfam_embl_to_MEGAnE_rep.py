#!/usr/bin/env python

'''
Copyright (c) 2020-2022 RIKEN
All Rights Reserved
See file LICENSE for details.

Prerequisite:
    Python3.7 or later.

Usage:
    Dfam_embl_to_MEGAnE_rep.py -i Dfam_target_species.embl
    
Details:
    First, please convert the "Dfam_curatedonly.h5" file to the embl format using the "famdb.py".
    e.g. in the case of human: python famdb.py -i /path/to/Dfam_curatedonly.h5 families -f embl -ad 9606 > Dfam_9606.embl
    "famdb.py" can be downloaded from https://github.com/Dfam-consortium/FamDB/blob/master/famdb.py
    You can convert the embl format file to a MEGAnE rep file.
    e.g. python Dfam_embl_to_MEGAnE_rep.py -i Dfam_9606.embl
    This will generate "Dfam_9606.rep" file which can be used in MEGAnE step 1.

Show help message:
    python Dfam_fa_to_rep.py -h
'''


import os,sys,collections,argparse

parser=argparse.ArgumentParser(description='')
parser.add_argument('-i', metavar='str', type=str, help='Required. Specify fasta file containing mobile elements. This must be generated from Dfam.')
parser.add_argument('-o', metavar='str', type=str, help='Optional. Specify output file name. Otherwise, the basename of the input file will be used as the output file basename.')
parser.add_argument('-do_not_overwrite', help='Optional. Specify if you do NOT want to overwrite previous results.', action='store_true')
args=parser.parse_args()


# initial check
if args.i is None:
    print('Error: You need to specify input .embl file with "-i" option.', file=sys.stderr)
    exit(1)
else:
    if os.path.exists(args.i) is not True:
        print('Error: Input file does not exist.', file=sys.stderr)
        exit(1)

if args.o is None:
    if args.i[-5:] == '.embl':
        args.o= args.i[:-5] + '.rep'
    else:
        args.o= args.i + '.rep'

if os.path.exists(args.o) is True:
    if args.do_not_overwrite is True:
        print('Error: Output file already exists. Please remove "-do_not_overwrite" option if you overwrite the existing file.', file=sys.stderr)
        exit(1)
    else:
        print('Info: Output file already existed, but overwrote it according to the absence of "-do_not_overwrite" option.')


remove_type={'tRNA', 'Unknown', 'snRNA', 'scRNA', 'Satellite', 'rRNA'}


def embl_parser(f):
    d={}
    # initialize
    seqname=None
    repbase=None
    type=None
    subtype=None
    seqline=False
    seq=[]
    # read
    with open(f) as infile:
        for line in infile:
            key=line[:2]
            if key == '//':
                if not seqname is None:
                    if type is None:
                        print('type is None')
                        exit()
                    if not type in remove_type:
                        if subtype is None:
                            family=type
                        else:
                            family='%s/%s' % (type, subtype)
                        d[seqname]=[repbase, family, ''.join(seq)]
                seqname=None
                repbase=None
                type=None
                subtype=None
                seqline=False
                seq=[]
                continue
            if seqline is True:
                _seq=''.join(line.split()[:-1])
                seq.append(_seq)
                continue
            if key == 'NM':
                seqname=line.split()[1]
            elif key == 'DR':
                repbase=line.strip().replace('DR   Repbase; ', '')[:-1]
            elif key == 'CC':
                if 'CC        Type: ' in line:
                    ls=line.split()
                    if len(ls) == 3:
                        type=ls[2]
                elif 'CC        SubType: ' in line:
                    ls=line.split()
                    if len(ls) == 3:
                        subtype=ls[2]
            elif key == 'SQ':
                seqline=True
    return d


def check_L1_seq_overlap(end5, orf2, end3):
    second_overlap_orf2=orf2[-150:]
    second_overlap_end3=end3[:150]
    first_overlap=0
    for n in range(380):
        pos=400 - n
        first_overlap_end5=end5[-pos:]
        first_overlap_orf2=orf2[:pos]
        if first_overlap_end5 == first_overlap_orf2:
            first_overlap=pos
            break
    if first_overlap == 0:
        return False, 0, 0
    second_overlap=0
    for n in range(380):
        pos=400 - n
        second_overlap_orf2=orf2[-pos:]
        second_overlap_end3=end3[:pos]
        if second_overlap_orf2 == second_overlap_end3:
            second_overlap=pos
            break
    if second_overlap == 0:
        return False, 0, 0
    return True, first_overlap, second_overlap


def reconstruct_L1(embl, out, headers):
    # pick up L1 with 5end, orf2, 3end
    L1_checker={}
    reconstructed=set()
    for seqname in embl:
        repbase,family,seq=embl[seqname]
        if family == 'LINE/L1':
            if '_5end' in seqname:
                L1_name=seqname.replace('_5end', '')
                if not L1_name in L1_checker:
                    L1_checker[L1_name]=[[False, False, False], ['', '', '']]
                L1_checker[L1_name][0][0]=True
                L1_checker[L1_name][1][0]=seq
            elif '_orf2' in seqname:
                L1_name=seqname.replace('_orf2', '')
                if not L1_name in L1_checker:
                    L1_checker[L1_name]=[[False, False, False], ['', '', '']]
                L1_checker[L1_name][0][1]=True
                L1_checker[L1_name][1][1]=seq
            elif '_3end' in seqname:
                L1_name=seqname.replace('_3end', '')
                if not L1_name in L1_checker:
                    L1_checker[L1_name]=[[False, False, False], ['', '', '']]
                L1_checker[L1_name][0][2]=True
                L1_checker[L1_name][1][2]=seq
    to_be_reconstructed=set()
    for L1_name in L1_checker:
        if sum(L1_checker[L1_name][0]) == 3:
            judge,first_overlap,second_overlap=check_L1_seq_overlap(*L1_checker[L1_name][1])
            if judge is False:
                continue
            seq=L1_checker[L1_name][1][0] + L1_checker[L1_name][1][1][first_overlap:] + L1_checker[L1_name][1][2][second_overlap:]
            out.append('>%s\t%s\t.\n%s\n' % (L1_name, 'LINE/L1', seq))
            headers[L1_name] += 1
            reconstructed.add(L1_name)
    return out, headers, reconstructed


def format(embl):
    out=[]
    headers=collections.Counter()
    
    # check redundant repbase name
    counter=collections.Counter()
    for seqname in embl:
        repbase,family,seq=embl[seqname]
        counter[repbase] += 1
    redundant=set()
    for repbase in counter:
        if counter[repbase] >= 2:
            redundant.add(repbase)
    
    # reconstruct L1
    out,headers,reconstructed=reconstruct_L1(embl, out, headers)
    print('Reconstructed L1: ', len(reconstructed))
    
    # others
    for seqname in embl:
        repbase,family,seq=embl[seqname]
        if repbase in reconstructed:
            continue
        if family is None:
            family=''
        else:
            family= '\t' + family
        out.append('>%s%s\t.\n%s\n' % (seqname, family, seq))
        headers[seqname] += 1
    
    # redundancy check
    for seqname in headers:
        if headers[seqname] >= 2:
            print('ERR: redundant: ', seqname)
            exit()
    return out


embl=embl_parser(args.i)
out=format(embl)
print('Repeats found: ', len(out))
with open(args.o, 'w') as outfile:
    outfile.write(''.join(out))

print('When running MEGAnE with this output, please specify %s with "-rep" flag.' % args.o)
