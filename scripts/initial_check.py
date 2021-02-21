#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,multiprocessing
from os.path import abspath,dirname,realpath,join
import log,traceback

# http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    elif 'PATH' in os.environ:
        for path in os.environ['PATH'].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                log.logger.debug('%s found: %s' % (program, exe_file))
                return exe_file
    return None

def check(args, argv):
    log.logger.debug('started')
    try:
        log.logger.debug('command line:\n'+ ' '.join(argv))
        # check python version
        version=sys.version_info
        if (version[0] >= 3) and (version[1] >= 7):
            log.logger.debug('Python version=%d.%d.%d' % (version[0], version[1], version[2]))
        else:
            log.logger.error('Please use Python 3.7 or later. Your Python is version %d.%d.' % (version[0], version[1]))
            exit(1)
        
        # check cpu num
        cpu_num=multiprocessing.cpu_count()
        if args.p > cpu_num:
            log.logger.error('Too many thread number. Please specify the number less than your cpu cores. You specified = %d, cpu cores = %d.' % (args.p, cpu_num))
            exit(1)
        
        # check PATH
        for i in ['blastn', 'bedtools', 'samtools']:
            if which(i) is None:
                log.logger.error('%s not found in $PATH. Please check %s is installed and added to PATH.' % (i, i))
                exit(1)
        
        # check file path
        for f_check,f_name in zip([args.fa, args.rep, args.repout, args.mainchr, args.repremove, args.pA_ME], ['Reference genome (-fa flag)', 'RepBase library (-rep flag)', 'RepeatMasker output file (-repout flag)', 'chr name file (-mainchr flag)', 'Non-ME rep list (-repremove flag)', 'ME with pA list (-pA_ME flag)']):
            if f_check is None:
                log.logger.error('%s was not specified.' % f_name)
                exit(1)
            elif os.path.exists(f_check) is False:
                log.logger.error('%s was not found.' % f_check)
                exit(1)
            else:
                log.logger.debug('%s found.' % f_check)
        
        import pysam
        if args.c is not None:
            if os.path.exists(args.c) is False:
                log.logger.error('%s was not found.' % args.c)
                exit(1)
            infile=pysam.AlignmentFile(args.c, 'rc', reference_filename=args.fa)
            for line in infile:
                line=line.tostring()
                break
            log.logger.debug('%s was able to open.' % args.c)
        elif args.b is not None:
            if os.path.exists(args.b) is False:
                log.logger.error('%s was not found.' % args.b)
                exit(1)
            infile=pysam.AlignmentFile(args.b, 'rb')
            for line in infile:
                line=line.tostring()
                break
            log.logger.debug('%s was able to open.' % args.b)
        else:
            log.logger.error('Please specify BAM or CRAM file with -b or -c flag.')
            exit(1)
        
        # check prerequisite modules
        from Bio.Seq import Seq
        import gzip
        from pybedtools import BedTool
        import matplotlib
        from Bio.Blast.Applications import NcbiblastnCommandline
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def check_merge_vcf(args, argv):
    log.logger.debug('started')
    try:
        log.logger.debug('command line:\n'+ ' '.join(argv))
        # check python version
        version=sys.version_info
        if (version[0] >= 3) and (version[1] >= 7):
            log.logger.debug('Python version=%d.%d.%d' % (version[0], version[1], version[2]))
        else:
            log.logger.error('Please use Python 3.7 or later. Your Python is version %d.%d.' % (version[0], version[1]))
            exit(1)
        
        # check cpu num
        cpu_num=multiprocessing.cpu_count()
        if args.p > cpu_num:
            log.logger.error('Too many thread number. Please specify the number less than your cpu cores. You specified = %d, cpu cores = %d.' % (args.p, cpu_num))
            exit(1)
        
        # check PATH
        for i in ['blastn', 'bedtools']:
            if which(i) is None:
                log.logger.error('%s not found in $PATH. Please check %s is installed and added to PATH.' % (i, i))
                exit(1)
        
        # check file path
        for f_check,f_name in zip([args.fa, args.rep, args.f], ['Reference genome (-fa flag)', 'RepBase library (-rep flag)', 'List of vcf files (-f flag)']):
            if f_check is None:
                log.logger.error('%s was not specified.' % f_name)
                exit(1)
            elif os.path.exists(f_check) is False:
                log.logger.error('%s was not found.' % f_check)
                exit(1)
            else:
                log.logger.debug('%s found.' % f_check)
        
        if args.merge_mei is True and args.merge_absent_me is True:
            log.logger.error('Please specify either -merge_mei or -merge_absent_me.')
            exit(1)
        elif args.merge_mei is False and args.merge_absent_me is False:
            log.logger.error('Please specify either -merge_mei or -merge_absent_me.')
            exit(1)
        
        if args.cohort_name is None:
            dt_now=datetime.datetime.now()
            args.cohort_name='%d-%d-%d-%d%d%d' % (dt_now.year, dt_now.month, dt_now.day, dt_now.hour, dt_now.minute, dt_now.second)
            log.logger.info('-cohort_name was not specified. Will use %s as a cohort name' % args.cohort_name)

        # check prerequisite modules
        from pybedtools import BedTool
        from Bio.Blast.Applications import NcbiblastnCommandline
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def check_reshape_vcf(args, argv):
    log.logger.debug('started')
    try:
        log.logger.debug('command line:\n'+ ' '.join(argv))
        # check python version
        version=sys.version_info
        if (version[0] >= 3) and (version[1] >= 7):
            log.logger.debug('Python version=%d.%d.%d' % (version[0], version[1], version[2]))
        else:
            log.logger.error('Please use Python 3.7 or later. Your Python is version %d.%d.' % (version[0], version[1]))
            exit(1)
                        
        for f_check,f_name in zip([args.i, args.a], ['VCF of MEI (-i flag)', 'VCF of absent ME (-a flag)']):
            if f_check is None:
                log.logger.error('%s was not specified.' % f_name)
                exit(1)
            elif os.path.exists(f_check) is False:
                log.logger.error('%s was not found.' % f_check)
                exit(1)
            else:
                log.logger.debug('%s found.' % f_check)
        
        if not '.vcf' in args.i:
            log.logger.error('please specify .vcf file with -i flag.')
            exit(1)
        if not '.vcf' in args.a:
            log.logger.error('please specify .vcf file with -a flag.')
            exit(1)
        
        if args.cohort_name is None:
            dt_now=datetime.datetime.now()
            args.cohort_name='%d-%d-%d-%d%d%d' % (dt_now.year, dt_now.month, dt_now.day, dt_now.hour, dt_now.minute, dt_now.second)
            log.logger.info('-cohort_name was not specified. Will use %s as a cohort name' % args.cohort_name)
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
