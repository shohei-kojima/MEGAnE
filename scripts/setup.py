#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,multiprocessing,time,random
from os.path import abspath,dirname,realpath,join
import log,traceback


def setup(args, base):
    log.logger.debug('started')
    try:
        sys.setrecursionlimit(10000)
        
        # tmp dir for pybedtools
        if args.pybedtools_tmp is None:
            args.pybedtools_tmp=args.outdir
        else:
            if os.path.exists(args.pybedtools_tmp) is False:
                log.logger.info('%s does not exist. Please check if %s exists.' % (args.pybedtools_tmp, args.pybedtools_tmp))
                exit(1)
            nowtime=str(time.time()).replace('.', '_')
            randval=random.randint(0, 1000000000)
            nowtime='%s_%s' % (nowtime, str(randval))
            args.pybedtools_tmp=os.path.join(args.pybedtools_tmp, nowtime)
            if os.path.exists(args.pybedtools_tmp) is True:
                log.logger.error('Exited due to a quite minor and corner reason. Please run again.')
                exit(1)
            else:
                os.mkdir(args.pybedtools_tmp)
        
        # load main chrs
        global main_chrs
        if args.mainchr is not None:
            main_chr_path=args.mainchr
        else:
            log.logger.error('Please specify a file containing chromosomes to be analyzed with "-mainchr" flag.')
            exit(1)
        main_chrs=[]
        with open(main_chr_path) as infile:
            for line in infile:
                chr=line.strip()
                if not chr == '':
                    main_chrs.append(chr)
        
        # check sex, cov, readlen
        global female, male, chrX, chrY, auto
        auto=args.auto
        female={'female', 'Female', 'FEMALE', 'F', 'f'}
        male={'male', 'Male', 'MALE', 'M', 'm'}
        chrX=set([ chr for chr in args.female_sex_chr.split(',') ])
        chrY=set([ chr for chr in args.male_sex_chr.split(',') ])
        args.sex=args.sex.replace('"', '')
        if not args.sex in auto:
            if not args.sex in female:
                if not args.sex in male:
                    if not args.sex.lower().replace('"', '') == 'unknown':
                        log.logger.error('Unknown value was specified with -sex. Please check again.')
                        exit(1)
        if not args.cov in auto:
            if isinstance(args.cov, int) is True:
                pass
            else:
                args.cov=args.cov.replace(',', '')
                tmp=args.cov.replace('.', '').replace('-', '')
                if tmp.isnumeric() is False:
                    log.logger.error('Unknown value was specified with -cov. Please check again.')
                    exit(1)
                if '.' in args.cov:
                    if args.cov.count('.') == 1:
                        args.cov=float(args.cov)
                        log.logger.info('-cov %f was rounded to %s' % (args.cov, round(args.cov)))
                        args.cov=round(args.cov)
                    else:
                        log.logger.error('Unknown value was specified with -cov. Please check again.')
                        exit(1)
                else:
                    args.cov=int(args.cov)
            if args.cov < 1:
                log.logger.error('Zero or negative value was specified with -cov. Please check again.')
                exit(1)
            elif args.cov < 5:
                log.logger.warning('< 5 was specified with -cov. Efficiency of polymorphic ME detection may not be good.')
        if not args.readlen in auto:
            if isinstance(args.readlen, int) is True:
                pass
            else:
                args.readlen=args.readlen.replace(',', '')
                tmp=args.readlen.replace('.', '').replace('-', '')
                if tmp.isnumeric() is False:
                    log.logger.error('Unknown value was specified with -readlen. Please check again.')
                    exit(1)
                if '.' in args.readlen:
                    if args.readlen.count('.') == 1:
                        args.readlen=float(args.readlen)
                        log.logger.info('-readlen %f was rounded to %s' % (args.readlen, round(args.readlen)))
                        args.readlen=round(args.readlen)
                    else:
                        log.logger.error('Unknown value was specified with -readlen. Please check again.')
                        exit(1)
                else:
                    args.readlen=int(args.readlen)
            if args.readlen < 1:
                log.logger.error('Zero or negative value was specified with -readlen. Please check again.')
                exit(1)
            elif args.readlen < 100:
                log.logger.warning('< 100 was specified with -readlen. Efficiency of polymorphic ME detection may not be good.')
        
        # chromosome check
        import pysam
        if args.c is not None:
            header=pysam.view(args.c, '-H')
            samfilename=args.c
        elif args.b is not None:
            header=pysam.view(args.b, '-H')
            samfilename=args.b
        all_chrs_in_bam_header=set()
        for line in header.split('\n'):
            if '@SQ\tSN:' in line:
                seqname=line.split()[1][3:]
                all_chrs_in_bam_header.add(seqname)
        log.logger.debug('%d chromosome(s) were found in %s.' % (len(all_chrs_in_bam_header), samfilename))
        # main chr
        found_n=0
        not_found=[]
        for chr in main_chrs:
            if chr in all_chrs_in_bam_header:
                found_n += 1
            else:
                not_found.append(chr)
        if len(not_found) >= 1:
            log.logger.error('Chromosomes %s were NOT found in %s. Please check again if you are specifying correct chrosomes with "-mainchr" flag.' % (','.join(not_found), samfilename))
            exit(1)
        else:
            log.logger.info('All %d main chromosome(s) were found in %s.' % (found_n, samfilename))
        # sex chr
        if args.no_sex_chr is False:
            found_n=0
            for chr in chrX:
                if chr in all_chrs_in_bam_header:
                    log.logger.info('"%s" was found in %s. "%s" will be considered as a female sex chromosome.' % (chr, samfilename, chr))
                    found_n += 1
            if found_n == 0:
                log.logger.error('Sex chromosome %s was NOT found in %s. Please check again if you are specifying correct chrosomes with "-female_sex_chr" flag.' % (chrX, samfilename))
                exit(1)
            found_n=0
            for chr in chrY:
                if chr in all_chrs_in_bam_header:
                    log.logger.info('"%s" was found in %s. "%s" will be considered as a male sex chromosome.' % (chr, samfilename, chr))
                    found_n += 1
                    args.nochrY=False
            if found_n == 0:
                args.nochrY=True
                log.logger.warning('Sex chromosome %s was NOT found in %s. Male sex chromosome will not be analyzed. Will continue anyway.' % (chrY, samfilename))
        else:
            args.nochrY=True  # this is needed because of auto_setting.estimate_depth_sex
            log.logger.info('"-no_sex_chr" option was specified. All chromosomes will be treated as autosomes.')
        
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
            param_path=join(base, 'docs/parameter_settings_lowdep.txt')
        elif args.lowdep is True:
            param_path=join(base, 'docs/parameter_settings_lowdep.txt')
        else:
            param_path=join(base, 'docs/parameter_settings.txt')
        import load_parameters
        params=load_parameters.load(args, param_path)

        # load rep headers to be removed
        global rep_headers_to_be_removed
        if args.repremove is not None:
            param_path=args.repremove
        else:
            log.logger.error('Please specify a file containing non-ME repeats with "-repremove" flag.')
            exit(1)
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
            log.logger.error('Please specify a file containing MEs with poly-A tail with "-pA_ME" flag.')
            exit(1)
        rep_with_pA=set()
        with open(param_path) as infile:
            for line in infile:
                line=line.strip().replace(' ', '_')
                if not line == '':
                    rep_with_pA.add(line)
        
        # homozygous
        if args.homozygous is True:
            args.monoallelic=True
        
        # L1 filter
        args.L1_filt_off=True  # default from MEGAnE 1.0.0
        if args.L1_filt_on is True:
            args.L1_filt_off=False
        
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def setup_geno_only_load_params(args, base):
    log.logger.debug('started')
    try:
        # load parameter settings
        global params
        param_path=join(base, 'docs/parameter_settings_geno.txt')
        import load_parameters
        params=load_parameters.load_geno(args, param_path)
        # check sex, cov, readlen
        global female, male, chrX, chrY, auto
        auto=args.auto
        female={'female', 'Female', 'FEMALE', 'F', 'f'}
        male={'male', 'Male', 'MALE', 'M', 'm'}
        chrX=set([ chr for chr in args.female_sex_chr.split(',') ])
        chrY=set([ chr for chr in args.male_sex_chr.split(',') ])
        if args.no_sex_chr is True:
            chrX=set()
            chrY=set()
        args.sex=args.sex.replace('"', '')
        if not args.sex in auto:
            if not args.sex in female:
                if not args.sex in male:
                    if not args.sex.lower().replace('"', '') == 'unknown':
                        log.logger.error('Unknown value was specified with -sex. Please check again.')
                        exit(1)
        if not args.cov in auto:
            if isinstance(args.cov, int) is True:
                pass
            else:
                args.cov=args.cov.replace(',', '')
                tmp=args.cov.replace('.', '').replace('-', '')
                if tmp.isnumeric() is False:
                    log.logger.error('Unknown value was specified with -cov. Please check again.')
                    exit(1)
                if '.' in args.cov:
                    if args.cov.count('.') == 1:
                        args.cov=float(args.cov)
                        log.logger.info('-cov %f was rounded to %s' % (args.cov, round(args.cov)))
                        args.cov=round(args.cov)
                    else:
                        log.logger.error('Unknown value was specified with -cov. Please check again.')
                        exit(1)
                else:
                    args.cov=int(args.cov)
            if args.cov < 1:
                log.logger.error('Zero or negative value was specified with -cov. Please check again.')
                exit(1)
            elif args.cov < 5:
                log.logger.warning('< 5 was specified with -cov. Efficiency of polymorphic ME detection may not be good.')
        if not args.readlen in auto:
            if isinstance(args.readlen, int) is True:
                pass
            else:
                args.readlen=args.readlen.replace(',', '')
                tmp=args.readlen.replace('.', '').replace('-', '')
                if tmp.isnumeric() is False:
                    log.logger.error('Unknown value was specified with -readlen. Please check again.')
                    exit(1)
                if '.' in args.readlen:
                    if args.readlen.count('.') == 1:
                        args.readlen=float(args.readlen)
                        log.logger.info('-readlen %f was rounded to %s' % (args.readlen, round(args.readlen)))
                        args.readlen=round(args.readlen)
                    else:
                        log.logger.error('Unknown value was specified with -readlen. Please check again.')
                        exit(1)
                else:
                    args.readlen=int(args.readlen)
            if args.readlen < 1:
                log.logger.error('Zero or negative value was specified with -readlen. Please check again.')
                exit(1)
            elif args.readlen < 100:
                log.logger.warning('< 100 was specified with -readlen. Efficiency of polymorphic ME detection may not be good.')
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
            log.logger.error('Please specify a file containing chromosomes to be analyzed with "-mainchr" flag.')
            exit(1)
        main_chrs=[]
        with open(main_chr_path) as infile:
            for line in infile:
                chr=line.strip()
                if not chr == '':
                    main_chrs.append(chr)
        global female, male, chrX, chrY
        female={'female', 'Female', 'F', 'f'}
        male={'male', 'Male', 'M', 'm'}
        chrX=set([ chr for chr in args.female_sex_chr.split(',') ])
        chrY=set([ chr for chr in args.male_sex_chr.split(',') ])
        args.sex=args.sex.replace('"', '')
#        if not args.sex in auto:
#            if not args.sex in female:
#                if not args.sex in male:
#                    if not args.sex.lower().replace('"', '') == 'unknown':
#                        log.logger.error('Unknown value was specified with -sex. Please check again.')
#                        exit(1)

        # fai
        global fai_path
        if os.path.exists(args.fa + '.fai') is False:
            import pysam
            pysam.faidx(args.fa)
            log.logger.info('Fasta index not found. Generated %s.fai.' % args.fa)
        fai_path=args.fa + '.fai'
        
        # load parameter settings
        global params
        param_path=join(base, 'docs/parameter_settings_geno.txt')
        import load_parameters
        params=load_parameters.load_geno(args, param_path)
                
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def setup_merge_vcf(args, base):
    log.logger.debug('started')
    try:
        # check options
        if args.input_scaffold is not None:
            if args.chunk_f is None:
                log.logger.error('Chunk file was not specified. Please specify a file with the "-chunk_f" flag.')
                exit(1)
            if args.f is not None:
                log.logger.warning('A file was specified with the "-f" flag. This will be ignored.')
        else:
            if args.chunk_f is not None:
                log.logger.warning('Chunk file was specified with the "-chunk_f" flag. This will be ignored.')
            if args.chunk_vcf_list is None and args.f is None:
                log.logger.error('Necessary file was not specified. Please specify a file with the "-f" flag.')
                exit(1)
        
        # tmp dir for pybedtools
        if args.pybedtools_tmp is None:
            args.pybedtools_tmp=args.outdir
        else:
            if os.path.exists(args.pybedtools_tmp) is False:
                log.logger.error('%s does not exist. Please check if %s exists.' % (args.pybedtools_tmp, args.pybedtools_tmp))
                exit(1)
            nowtime=str(time.time()).replace('.', '_')
            randval=random.randint(0, 1000000000)
            nowtime='%s_%s' % (nowtime, str(randval))
            args.pybedtools_tmp=os.path.join(args.pybedtools_tmp, nowtime)
            if os.path.exists(args.pybedtools_tmp) is True:
                log.logger.error('Exited due to a quite minor and corner reason. Please run again.')
                exit(1)
            else:
                os.mkdir(args.pybedtools_tmp)
        
        # load main chrs
        global female, male, chrX, chrY
        female={'female', 'Female', 'F', 'f'}
        male={'male', 'Male', 'M', 'm'}
        if args.no_sex_chr is False:
            chrX=set([ chr for chr in args.female_sex_chr.split(',') ])
            chrY=set([ chr for chr in args.male_sex_chr.split(',') ])
        else:
            chrX=set()
            chrY=set()
        if args.chr is not None:
            args.chr=set([ chr for chr in args.chr.split(',') ])
        
        # load rep headers to be removed
        global rep_headers_to_be_removed
        if args.repremove is not None:
            param_path=args.repremove
        else:
            param_path=join(base, 'docs/human_non_ME_rep_headers.txt')
        rep_headers_to_be_removed=set()
        with open(param_path) as infile:
            for line in infile:
                line=line.strip().replace(' ', '_')
                if not line == '':
                    rep_headers_to_be_removed.add(line)
        
        # fai
        global fai_path
        if args.chunk_vcf_list is None:
            if os.path.exists(args.fa + '.fai') is False:
                import pysam
                pysam.faidx(args.fa)
                log.logger.info('Fasta index not found. Generated %s.fai.' % args.fa)
            fai_path=args.fa + '.fai'
        else:
            fai_path='dummy'
        
        # load parameter settings
        global params
        import load_parameters
        params=load_parameters.load_merge_vcf(args)
                
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
