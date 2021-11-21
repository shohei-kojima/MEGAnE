#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os
import pysam
import numpy as np
import log,traceback


def init(args):
    if args.sex.lower().replace('"', '') == 'unknown' and args.make_sex_auto is True:
        args.sex='auto'
        log.logger.warning('Sex is not specified. It will be automatically estimated. Please make sure the input must be a human sample.')
    elif args.sex.lower().replace('"', '') == 'unknown':
        log.logger.warning('Sex is not specified. MEGAnE consider all sex chromosomes as diploid. This option is NOT recommended. Please use this option at your own risk. Please specify sex whenever possible.')
    return {'auto', 'Auto', 'AUTO'}


def estimate_readlen(args):
    log.logger.debug('started')
    try:
        if args.c is not None:
            infile=pysam.AlignmentFile(args.c, 'rc', reference_filename=args.fa)
        else:
            infile=pysam.AlignmentFile(args.b, 'rb')
        lens=[]
        for read in infile:
            if read.is_secondary is False and read.is_supplementary is False and read.is_duplicate is False:
                lens.append(read.query_length)
                if len(lens) == 200:
                    break
        avelen=round(np.mean(lens))
        args.readlen=int(avelen)
        log.logger.debug('avelen=%d;lens=%s' % (avelen, ','.join([ str(v) for v in lens ])))
        log.logger.info('estimated read lenth = %d' % avelen)
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def estimate_depth_sorted(infile, chr, start, end):
    lens=0
    for read in infile.fetch(chr, start, end):
        if read.is_secondary is False and read.is_supplementary is False and read.is_duplicate is False:
            lens += read.query_length
    depth= round(lens / (end - start))
    return depth


def estimate_depth_unsorted(args, infile, chrx, chry):
    fai=0
    with open(args.fai) as faifile:
        for line in faifile:
            ls=line.split()
            fai += int(ls[1])
    readnum=0
    xnum=0
    ynum=0
    for read in infile:
        if read.is_secondary is False and read.is_supplementary is False and read.is_duplicate is False and read.is_unmapped is False:
            readnum += 1
            if read.reference_name == chrx:
                xnum += 1
            elif read.reference_name == chry:
                ynum += 1
    depth= round(1.1 * args.readlen * readnum / fai)
    return depth,xnum,ynum


def estimate_depth_sex(args, params, auto, read_cnts):
    log.logger.debug('started')
    try:
        # load fai
        nonXY_len=0
        X_len=0
        Y_len=0
        with open(args.fai) as infile:
            for line in infile:
                ls=line.split()
                if ls[0] in params.chrX:
                    X_len += int(ls[1])
                elif ls[0] in params.chrY:
                    Y_len += int(ls[1])
                else:
                    nonXY_len += int(ls[1])
        if args.nochrY is False:
            if X_len == 0 or Y_len == 0:
                log.logger.error('Sex chromosome not found in the reference fasta index. Please check again if the reference fasta contains sex chromosomes.')
                exit(1)
        else:
            if X_len == 0:
                log.logger.error('Sex chromosome not found in the reference fasta index. Please check again if the reference fasta contains sex chromosomes.')
                exit(1)
        if nonXY_len == 0:
            log.logger.error('Autosomes not found in the reference fasta index. Please check again if the reference fasta contains autosomes.')
            exit(1)
        depth= (args.readlen * read_cnts[0]) / nonXY_len
        x_depth= (args.readlen * read_cnts[1]) / X_len
        if args.nochrY is False:
            y_depth= (args.readlen * read_cnts[2]) / Y_len
            sex='female' if (y_depth / x_depth) < params.sex_est_XY_ratio_threshold else 'male'
            log.logger.debug('depth=%d;x_depth=%f;y_depth=%f;sex=%s;nonXY_len=%d' % (depth, x_depth, y_depth, sex, nonXY_len))
        else:
            sex='female' if (x_depth / depth) > params.sex_est_XY_ratio_threshold_for_nochrY else 'male'
            log.logger.debug('depth=%d;x_depth=%f;y_depth=%s;sex=%s;nonXY_len=%d' % (depth, x_depth, 'NA', sex, nonXY_len))
        if args.cov in auto:
            args.cov= int(round(depth))
            log.logger.info('estimated autosome depth = %d' % int(round(depth)))
        if args.sex in auto:
            args.sex=sex
            log.logger.info('estimated sex = %s' % sex)
            if args.nochrY is True:
                log.logger.warning('Sex estimation may not correct, due to absence of male sex chromosome in the input BAM/CRAM file. Please use this result at your own risk. Specifying sex with "-sex" flag is highly recommended.')
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def estimate_depth_sex_only_human(args, params, auto):  # deprecated
    log.logger.debug('started')
    try:
        if args.c is not None:
            infile=pysam.AlignmentFile(args.c, 'rc', reference_filename=args.fa)
        else:
            infile=pysam.AlignmentFile(args.b, 'rb')
        header=infile.header.to_dict()
        # determine chrX and chrY names
        chr_names=set()
        for d in header['SQ']:
            chr_names.add(d['SN'])
        chrx,chry=0,0
        if 'chrX' in chr_names:
            chrx='chrX'
        elif 'X' in chr_names:
            chrx='X'
        if 'chrY' in chr_names:
            chry='chrY'
        elif 'Y' in chr_names:
            chry='Y'
        if chrx == 0:
            log.logger.error('chrX or X not found in BAM/CRAM header. Cannot estimate sex of this sample. Exit.')
            exit(1)
        if chry == 0:
            log.logger.error('chrY or Y not found in BAM/CRAM header. Cannot estimate sex of this sample. Exit.')
            exit(1)
        # determine depth and sex
        if header['HD']['SO'] == 'coordinate':
            log.logger.debug('Input BAM/CRAM is coordinate-sorted.')
            if args.b is not None:
                if os.path.exists(args.b + '.bai') is False:
                    base,_=os.path.splitext(args.b)
                    if os.path.exists(base + '.bai') is False:
                        log.logger.info('Generating BAM index...')
                        pysam.index(args.b, '-@ %s' % args.p)
            else:
                if os.path.exists(args.c + '.crai') is False:
                    base,_=os.path.splitext(args.c)
                    if os.path.exists(base + '.crai') is False:
                        log.logger.info('Generating CRAM index...')
                        pysam.index(args.c, '-@ %s' % args.p)
            chr1=header['SQ'][0]['SN']
            depth=estimate_depth_sorted(infile, chr1, params.chr1_start_depth_est, params.chr1_end_depth_est)
            x_depth=estimate_depth_sorted(infile, chrx, params.chrX_start_depth_est, params.chrX_end_depth_est)
            y_depth=estimate_depth_sorted(infile, chry, params.chrY_start_depth_est, params.chrY_end_depth_est)
        else:
            log.logger.warning('Input BAM/CRAM is not coordinate-sorted. MEGAnE will start autosome depth and sex prediction, but this will take a longer time. Please consider to sort input alignment file.')
            depth,xnum,ynum=estimate_depth_unsorted(args, infile, chrx, chry)
            for d in header['SQ']:
                if d['SN'] == chrx:
                    xlen=d['LN']
                elif d['SN'] == chry:
                    ylen=d['LN']
            log.logger.debug('xnum=%d;ynum=%d;xlen=%d;ylen=%d' % (xnum, ynum, xlen, ylen))
            x_depth= args.readlen * xnum / xlen
            y_depth= args.readlen * ynum / ylen
        if x_depth > 0:
            sex='female' if (y_depth / x_depth) < params.sex_est_XY_ratio_threshold else 'male'
        else:
            log.logger.debug('depth=%d;x_depth=%f;y_depth=%f' % (depth, x_depth, y_depth))
            log.logger.error('Depth of %s:%d-%d is 0. Cannot estimate sex. Please try to determine the sex by different way. Exit.' % (chrx, params.chrX_start_depth_est, params.chrX_end_depth_est))
            exit(1)
        log.logger.debug('depth=%d;x_depth=%f;y_depth=%f;sex=%s' % (depth, x_depth, y_depth, sex))
        if args.cov in auto:
            args.cov=depth
            log.logger.info('estimated autosome depth = %d' % depth)
        if args.sex in auto:
            args.sex=sex
            log.logger.info('estimated sex = %s' % sex)
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

