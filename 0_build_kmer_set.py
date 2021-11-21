#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020-2021 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,sys,datetime,argparse,glob,shutil,logging
import ctypes as ct

# version
version='v1.0.0 2021/11/21'

# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-fa', metavar='str', type=str, help='Fasta file of the reference genome (e.g. GRCh38.fa). Please make sure there is a fasta index (.fa.fai)', required=True)
parser.add_argument('-prefix', metavar='str', type=str, help='Optional. Specify prefix of the k-mer file. If this was not specified, prefix will be the file name of the input reference genome. Two files, [prefix].mk and [prefix].mi, will be generated.')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./megane_kmer_set', default='./megane_kmer_set')
parser.add_argument('-v', '--version', action='version', version='MEGAnE %s %s' % (os.path.basename(__file__), version))
args=parser.parse_args()
args.version=version


# start
import init
init.init_build_kmer(args, version)


# logging
import log
args.logfilename='megane_build_kmer_set.log'
if os.path.exists(os.path.join(args.outdir, args.logfilename)) is True:
    os.remove(os.path.join(args.outdir, args.logfilename))
log.start_log(args)
log.logger.debug('Logging started.')


# initial check
import initial_check
log.logger.debug('This is %s version %s' % (__file__, version))
print()
log.logger.info('Initial check started.')
initial_check.check_build_kmer(args, sys.argv)


# output file names
import utils
filenames=utils.empclass()
filenames.mk =os.path.join(args.outdir, args.prefix + '.mk')
filenames.mi =os.path.join(args.outdir, args.prefix + '.mi')


# build k-mer set
stdout=os.dup(1)  # copy stdout stream
os.dup2(logging._handlerList[2]().stream.fileno(), 1)  # change stdout stream to logger
so='%s/cpp/save_redundant_kmers.so' % init.base
cpp=ct.CDLL(so)

def to_char_p(char):
    return ct.c_char_p(char.encode('utf-8'))

res=cpp.find_and_save_red_kmers(to_char_p(args.fa), to_char_p(args.fa + '.fai'), to_char_p(filenames.mk), to_char_p(filenames.mi))
os.dup2(stdout, 1)  # reset stdout stream
log.logger.debug('exit status of cpp.find_and_save_red_kmers: %s' % res)

