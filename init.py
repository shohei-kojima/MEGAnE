#!/usr/bin/env python

'''
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
            os.mkdir(args.outdir)
    else:
        if os.path.exists(args.outdir) is False:
            os.mkdir(args.outdir)


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
            os.mkdir(args.outdir)
    else:
        if os.path.exists(args.outdir) is False:
            os.mkdir(args.outdir)


def init_jointcall(args, version):
    # pythonpath
    global base
    base=os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    sys.path.insert(0, os.path.join(base, 'scripts'))
    
    # make output dir
    if args.chr is not None and args.outdir == './jointcall_out':
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
            os.mkdir(args.outdir)
    else:
        if os.path.exists(args.outdir) is False:
            os.mkdir(args.outdir)


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
            os.mkdir(args.outdir)
    else:
        if os.path.exists(args.outdir) is False:
            os.mkdir(args.outdir)
