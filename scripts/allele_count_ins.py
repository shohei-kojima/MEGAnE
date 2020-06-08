#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,gzip
import pysam
from pybedtools import BedTool
import log,traceback


def limit(args, params, filenames):
    log.logger.debug('started')
    try:
        if args.ins_bed is True and args.abs_bed is True:
            insbed=BedTool(args.ins_bed).slop(b=params.ins_slop_len, g=args.fai)
            absbed=BedTool(args.abs_bed).slop(b=params.abs_slop_len, g=args.fai)
            slopbed= insbed + absbed
            slopbed=slopbed.merge()
        elif args.ins_bed is True:
            slopbed=BedTool(args.ins_bed).slop(b=params.ins_slop_len, g=args.fai).merge()
        else:
            slopbed=BedTool(args.abs_bed).slop(b=params.abs_slop_len, g=args.fai).merge()
        if args.b is not None:
            pysam.view('-@', '%d' % args.p, '-bh', '-f', '2', '-M', '-L', slopbed, '-o', filenames.limited, args.b, catch_stdout=False)
        else:
            pysam.view('-@', '%d' % args.p, '-Ch', '-f', '2', '-M', '-L', slopbed, '-T', args.fa, '-o', filenames.limited, args.b, catch_stdout=False)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


# main
def genotype_ins(args, params, filenames, n):
    log.logger.debug('started')
    try:
        insbed=BedTool(args.ins_bed).slop(b=params.ins_slop_len, g=args.fai).merge()
        # convert to ins-only
        if args.b is not None:
            pysam.view('-@', '%d' % args.p, '-bh', '-f', '2', '-M', '-L', insbed, '-o', filenames.limited_ins, filenames.limited, catch_stdout=False)
        else:
            pysam.view('-@', '%d' % args.p, '-Ch', '-f', '2', '-M', '-L', insbed, '-T', args.fa, '-o', filenames.limited_ins, filenames.limited, catch_stdout=False)
        # convert to depth
        pysam.depth(filenames.limited_ins, '-o', filenames.depth_ins)

