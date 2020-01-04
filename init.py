#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys


def init(args):
    global base
    base=os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    sys.path.insert(0, os.path.join(base, 'scripts'))

    # make output dir
    if args.overwrite is False:
        if os.path.exists(args.outdir) is True:
            print('Error: %s already exists. Please specify another directory name.' % args.outdir)
            exit()
        else:
            os.mkdir(args.outdir)
    else:
        if os.path.exists(args.outdir) is False:
            os.mkdir(args.outdir)
