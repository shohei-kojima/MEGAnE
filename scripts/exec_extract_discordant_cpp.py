#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020-2021 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,shutil,logging
import log,traceback
import ctypes as ct

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
            argv=(ct.c_char_p * 6)('dummy'.encode('utf-8'),
                                   args.b.encode('utf-8'),
                                   args.mainchr.encode('utf-8'),
                                   filenames.reshaped_rep_mk.encode('utf-8'),
                                   args.outdir.encode('utf-8'),
                                   str(args.p).encode('utf-8'))
            res=cpp.main(len(argv), argv)
        else:
            argv=(ct.c_char_p * 7)('dummy'.encode('utf-8'),
                                   args.b.encode('utf-8'),
                                   args.mainchr.encode('utf-8'),
                                   filenames.reshaped_rep_mk.encode('utf-8'),
                                   args.outdir.encode('utf-8'),
                                   str(args.p).encode('utf-8'),
                                   args.fa.encode('utf-8'))
            res=cpp.main(len(argv), argv)
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
            argv=(ct.c_char_p * 5)('dummy'.encode('utf-8'),
                                   args.b.encode('utf-8'),
                                   filenames.reshaped_rep_mk.encode('utf-8'),
                                   args.outdir.encode('utf-8'),
                                   str(args.p).encode('utf-8'))
            res=cpp.main(len(argv), argv)
        else:
            argv=(ct.c_char_p * 6)('dummy'.encode('utf-8'),
                                   args.b.encode('utf-8'),
                                   filenames.reshaped_rep_mk.encode('utf-8'),
                                   args.outdir.encode('utf-8'),
                                   str(args.p).encode('utf-8'),
                                   args.fa.encode('utf-8'))
            res=cpp.main(len(argv), argv)
        os.dup2(stdout, 1)  # reset stdout stream
        log.logger.debug('exit status of cpp.main: %s' % res)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

