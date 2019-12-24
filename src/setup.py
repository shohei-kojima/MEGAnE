#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,multiprocessing
from os.path import abspath,dirname,realpath,join


def setup(args, base):
    # start
    print('Setup started: '+ datetime.datetime.now().isoformat())
    
    # make output dir
    if args.overwrite is False:
        if os.path.exists(args.outdir) is True:
            print('Error: %s already exists. Please specify another directory name.' % outdir)
            exit()
        else:
            os.mkdir(args.outdir)
    else:
        if os.path.exists(args.outdir) is False:
            os.mkdir(args.outdir)

    # load main chrs
    global main_chrs
    if args.mainchr is not None:
        main_chr_path=args.mainchr
    else:
        main_chr_path=join(base, 'lib/human_main_chrs.txt')
    import load_main_chrs
    main_chrs=load_main_chrs.load(main_chr_path)
    if args.mainchr is None:
        if args.gender == 'female':
            main_chrs=main_chrs[:-1]
            print('You specified "female" as gender, chrY was excluded.')

    # load parameter settings
    global params
    if args.setting is not None:
        param_path=args.setting
    else:
        param_path=join(base, 'lib/parameter_settings.txt')
    import load_parameters
    params=load_parameters.load(param_path)
