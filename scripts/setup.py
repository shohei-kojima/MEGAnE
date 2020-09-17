#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,multiprocessing
from os.path import abspath,dirname,realpath,join
import log,traceback


def setup(args, base):
    log.logger.debug('started')
    try:
        sys.setrecursionlimit(10000)
        
        # tmp dir for pybedtools
        args.pybedtools_tmp=args.outdir
        
        # load main chrs
        global main_chrs
        if args.mainchr is not None:
            main_chr_path=args.mainchr
        else:
            main_chr_path=join(base, 'lib/hg38_human_main_chrs_plus_alt_ucsc_style.txt')
        main_chrs=[]
        with open(main_chr_path) as infile:
            for line in infile:
                chr=line.strip()
                if not chr == '':
                    main_chrs.append(chr)
        global female, male, chrX, chrY
        female={'female', 'Female', 'F', 'f'}
        male={'male', 'Male', 'M', 'm'}
        chrX={'chrX', 'X', '23'}
        chrY={'chrY', 'Y', '24'}
#        if not args.sex == 'male':
#            if args.sex in female:
#                new_main_chr=[]
#                for chr in main_chrs:
#                    if not chr in chrY:
#                        new_main_chr.append(chr)
#                main_chrs=[]
#                main_chrs.extend(new_main_chr)
#                log.logger.info('You specified "female" as sex, chrY is going to be excluded from analysis.')
#            else:
#                log.logger.error('Unknown value was specified with -sex. Please check again.')
#                exit(1)
        
        # fai
        global fai_path
        if os.path.exists(args.fa + '.fai') is False:
            import pysam
            pysam.faidx(args.fa)
            log.logger.info('Fasta index not found. Generated %s.fai.' % args.fa)
        else:
            diff= os.stat(args.fa + '.fai').st_mtime - os.stat(args.fa).st_mtime
            if not diff > 0:
                log.logger.error('Fasta index %s.fai was older than fasta. Please generate fasta index again.' % args.fa)
                exit(1)
        fai_path=args.fa + '.fai'
        
        # load parameter settings
        global params
        if args.setting is not None:
            param_path=args.setting
        elif args.verylowdep is True:
            param_path=join(base, 'lib/parameter_settings_lowdep.txt')
        elif args.lowdep is True:
            param_path=join(base, 'lib/parameter_settings_lowdep.txt')
        else:
            param_path=join(base, 'lib/parameter_settings.txt')
        import load_parameters
        params=load_parameters.load(args, param_path)

        # load rep headers to be removed
        global rep_headers_to_be_removed
        if args.repremove is not None:
            param_path=args.repremove
        else:
            param_path=join(base, 'lib/human_non_ME_rep_headers.txt')
        rep_headers_to_be_removed=set()
        with open(param_path) as infile:
            for line in infile:
                line=line.strip().replace(' ', '_')
                if not line == '':
                    rep_headers_to_be_removed.add(line)

        # load rep with polyA tail
        global rep_with_pA
        if args.pA_ME is not None:
            param_path=args.pA_ME
        else:
            param_path=join(base, 'lib/human_ME_with_polyA_tail.txt')
        rep_with_pA=set()
        with open(param_path) as infile:
            for line in infile:
                line=line.strip().replace(' ', '_')
                if not line == '':
                    rep_with_pA.add(line)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def setup_geno_only_load_params(args, base):
    log.logger.debug('started')
    try:
        # load parameter settings
        global params
        param_path=join(base, 'lib/parameter_settings_geno.txt')
        import load_parameters
        params=load_parameters.load_geno(args, param_path)
        global female, male, chrX, chrY
        female={'female', 'Female', 'F', 'f'}
        male={'male', 'Male', 'M', 'm'}
        chrX={'chrX', 'X', '23'}
        chrY={'chrY', 'Y', '24'}
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def setup_geno(args, base):
    log.logger.debug('started')
    try:
        # potential 3' transduction
        if not args.abs_3t_bed is None:
            if args.abs_bed is None:
                log.logger.error('Cannot run genotyping of %s without -abs_bed option.' % args.abs_3t_bed)
                exit(1)
        
        # load main chrs
        global main_chrs
        if args.mainchr is not None:
            main_chr_path=args.mainchr
        else:
            main_chr_path=join(base, 'lib/hg38_human_main_chrs_plus_alt_ucsc_style.txt')
        main_chrs=[]
        with open(main_chr_path) as infile:
            for line in infile:
                chr=line.strip()
                if not chr == '':
                    main_chrs.append(chr)
        global female, male, chrX, chrY
        female={'female', 'Female', 'F', 'f'}
        male={'male', 'Male', 'M', 'm'}
        chrX={'chrX', 'X', '23'}
        chrY={'chrY', 'Y', '24'}
#        if not args.sex == 'male':
#            if args.sex in female:
#                new_main_chr=[]
#                for chr in main_chrs:
#                    if not chr in chrY:
#                        new_main_chr.append(chr)
#                main_chrs=[]
#                main_chrs.extend(new_main_chr)
#                log.logger.info('You specified "female" as sex, chrY is going to be excluded from analysis.')
#            else:
#                log.logger.error('Unknown value was specified with -sex. Please check again.')
#                exit(1)
        
        # fai
        global fai_path
        if os.path.exists(args.fa + '.fai') is False:
            import pysam
            pysam.faidx(args.fa)
            log.logger.info('Fasta index not found. Generated %s.fai.' % args.fa)
        fai_path=args.fa + '.fai'
        
        # load parameter settings
        global params
        param_path=join(base, 'lib/parameter_settings_geno.txt')
        import load_parameters
        params=load_parameters.load_geno(args, param_path)
                
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
