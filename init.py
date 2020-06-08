#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys


def init(args, version):
    # version check
    if args.version is True:
        print('MEI pipeline\nVersion = %s' % version)
        exit(0)
    
    # pythonpath
    global base
    base=os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    sys.path.insert(0, os.path.join(base, 'scripts'))

    # make output dir
    if args.overwrite is False:
        if os.path.exists(args.outdir) is True:
            print('Error: %s already exists. Please specify another directory name.' % args.outdir)
            exit(1)
        else:
            os.mkdir(args.outdir)
    else:
        if os.path.exists(args.outdir) is False:
            os.mkdir(args.outdir)


def init_geno(args, version):
    # version check
    if args.version is True:
        print('MEI genotyping pipeline\nVersion = %s' % version)
        exit(0)

    # pythonpath
    global base
    base=os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    sys.path.insert(0, os.path.join(base, 'scripts'))

    # make output dir
    if args.overwrite is False:
        if os.path.exists(args.outdir) is True:
            print('Error: %s already exists. Please specify another directory name.' % args.outdir)
            exit(1)
        else:
            os.mkdir(args.outdir)
    else:
        if os.path.exists(args.outdir) is False:
            os.mkdir(args.outdir)
