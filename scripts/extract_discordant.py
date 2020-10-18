#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,pysam,itertools,math,shutil,string
import log,traceback


# flagstat
def flagstat(args):
    log.logger.debug('started')
    try:
        flag=pysam.flagstat(args.b, '-@ %d' % (args.p - 1))  # multi process
        count= int(flag.split()[0])
        interval= math.ceil(count / args.p)
        return count,interval
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


# concatenate result files
def concat_for_ins(args, filenames):
    log.logger.debug('started')
    try:
        outfiles=[filenames.overhang_fa, filenames.overhang_pA, filenames.distant_txt, filenames.unmapped_fa, filenames.mapped_fa, filenames.blast1_res]
        for file_base in outfiles:
            shutil.move(file_base +'0.txt', file_base)
            with open(file_base, 'a') as outfile:
                for n in range(1, args.p):
                    with open(file_base + str(n) +'.txt') as infile:
                        shutil.copyfileobj(infile, outfile)
                    os.remove(file_base + str(n) +'.txt')
                outfile.flush()
                os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def concat_for_abs(args, filenames):
    log.logger.debug('started')
    try:
        outfiles=[filenames.abs_txt]
        for file_base in outfiles:
            shutil.move(file_base +'0.txt', file_base)
            with open(file_base, 'a') as outfile:
                for n in range(1, args.p):
                    with open(file_base + str(n) +'.txt') as infile:
                        shutil.copyfileobj(infile, outfile)
                    os.remove(file_base + str(n) +'.txt')
                outfile.flush()
                os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

