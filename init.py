#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,glob


def init(args, version):
    # pythonpath
    global base
    base=os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    sys.path.insert(0, os.path.join(base, 'scripts'))

    # make output dir
    if args.do_not_overwrite is True:
        if os.path.exists(args.outdir) is True:
            if args.outdir[-1] == '/':
                dir=args.outdir[:-1]
            else:
                dir=args.outdir
            if len(glob.glob('%s/*' % dir)) >= 1:
                print('Error: %s already exists. Please specify another directory name.' % args.outdir, file=sys.stderr)
                exit(1)
        else:
            os.makedirs(args.outdir, exist_ok=True)
    else:
        if os.path.exists(args.outdir) is False:
            os.makedirs(args.outdir, exist_ok=True)


def init_geno(args, version):
    # pythonpath
    global base
    base=os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    sys.path.insert(0, os.path.join(base, 'scripts'))

    # make output dir
    if args.do_not_overwrite is True:
        if os.path.exists(args.outdir) is True:
            if args.outdir[-1] == '/':
                dir=args.outdir[:-1]
            else:
                dir=args.outdir
            if len(glob.glob('%s/*' % dir)) >= 1:
                print('Error: %s already exists. Please specify another directory name.' % args.outdir, file=sys.stderr)
                exit(1)
        else:
            os.makedirs(args.outdir, exist_ok=True)
    else:
        if os.path.exists(args.outdir) is False:
            os.makedirs(args.outdir, exist_ok=True)


def init_jointcall(args, version):
    # pythonpath
    global base
    base=os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    sys.path.insert(0, os.path.join(base, 'scripts'))
    
    # make output dir
    if args.chunk_vcf_list is not None:
        args.outdir=os.path.join(args.outdir, 'chunks_merged')
    if args.input_scaffold is not None and args.chunk_f is not None:
        with open(args.chunk_f) as infile:
            first_sample=next(infile).split()[0]
            for line in infile:
                pass
        last_sample=line.split()[0]
        dir_name='%s_to_%s' % (first_sample, last_sample)
        args.outdir=os.path.join(args.outdir, dir_name)
    if args.chr is not None:
        args.outdir=os.path.join(args.outdir, args.chr.replace(',', '_'))
    if args.do_not_overwrite is True:
        if os.path.exists(args.outdir) is True:
            if args.outdir[-1] == '/':
                dir=args.outdir[:-1]
            else:
                dir=args.outdir
            if len(glob.glob('%s/*' % dir)) >= 1:
                print('Error: %s already exists. Please specify another directory name.' % args.outdir, file=sys.stderr)
                exit(1)
        else:
            os.makedirs(args.outdir, exist_ok=True)
    else:
        if os.path.exists(args.outdir) is False:
            os.makedirs(args.outdir, exist_ok=True)


def init_reshape_vcf(args, version):
    # pythonpath
    global base
    base=os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    sys.path.insert(0, os.path.join(base, 'scripts'))

    # make output dir
    if args.do_not_overwrite is True:
        if os.path.exists(args.outdir) is True:
            if args.outdir[-1] == '/':
                dir=args.outdir[:-1]
            else:
                dir=args.outdir
            if len(glob.glob('%s/*' % dir)) >= 1:
                print('Error: %s already exists. Please specify another directory name.' % args.outdir, file=sys.stderr)
                exit(1)
        else:
            os.makedirs(args.outdir, exist_ok=True)
    else:
        if os.path.exists(args.outdir) is False:
            os.makedirs(args.outdir, exist_ok=True)


def init_build_kmer(args, version):
    # pythonpath
    global base
    base=os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    sys.path.insert(0, os.path.join(base, 'scripts'))
    
    # make output dir
    os.makedirs(args.outdir, exist_ok=True)
