#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020-2021 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,shutil
import log,traceback

def char_p(char):
    return ct.c_char_p(char.encode('utf-8'))

def extract_discordant(args, params, filenames, init_base):
    log.logger.debug('started')
    try:
        stdout=os.dup(1)  # copy stdout stream
        os.dup2(logging._handlerList[2]().stream.fileno(), 1)  # change stdout stream to logger
        so='%s/cpp/extract_discordant.so' % init_base
        cpp=ct.CDLL(so)
        if args.b is not None:
            res=cpp.main(char_p(args.b),
                         char_p(args.mainchr),
                         char_p(args.mk),
                         char_p(args.outdir),
                         char_p(args.p))
        else:
            res=cpp.main(char_p(args.c),
                         char_p(args.mainchr),
                         char_p(args.mk),
                         char_p(args.outdir),
                         char_p(args.p),
                         char_p(args.fa))
        os.dup2(stdout, 1)  # reset stdout stream
        log.logger.debug('exit status of cpp.main: %s' % res)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def extract_unmapped(args, params, filenames, init_base):
    log.logger.debug('started')
    try:
        stdout=os.dup(1)  # copy stdout stream
        os.dup2(logging._handlerList[2]().stream.fileno(), 1)  # change stdout stream to logger
        so='%s/cpp/extract_unmapped.so' % init_base
        cpp=ct.CDLL(so)
        if args.b is not None:
            res=cpp.main(char_p(args.b),
                         char_p(args.mk),
                         char_p(args.outdir),
                         char_p(args.p))
        else:
            res=cpp.main(char_p(args.c),
                         char_p(args.mk),
                         char_p(args.outdir),
                         char_p(args.p),
                         char_p(args.fa))
        os.dup2(stdout, 1)  # reset stdout stream
        log.logger.debug('exit status of cpp.main: %s' % res)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

