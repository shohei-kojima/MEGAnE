#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,shutil,logging


# version
version='2021/01/06'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-merge_mei', help='Specify if you merge MEIs (e.g. MEI_final_gaussian_genotyped.vcf).', action='store_true')
parser.add_argument('-merge_absent_me', help='Specify if you merge absent MEs (e.g. absent_ME_final_genotyped.vcf).', action='store_true')
parser.add_argument('-f', metavar='str', type=str, help='Required. Specify paths to the vcf files to be merged. One line should contain one path to a vcf file.')
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: GRCh38DH.fa')
parser.add_argument('-rep', metavar='str', type=str, help='Required. Specify RepBase file used for repeatmasking. Example: humrep.ref')
parser.add_argument('-repremove', metavar='str', type=str, help='Optional. Specify full path to a file containing the names of non-ME repeat class. Default: /path/to/prog/lib/human_non_ME_rep_headers.txt')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./jointcall_out', default='./jointcall_out')
parser.add_argument('-cohort_name', metavar='str', type=str, help='Optional. Specify a cohort name. This will be used for the variant names as well the output file name. Default: YYYY-MM-DD-HHMMSS')
parser.add_argument('-overwrite', help='Optional. Specify if you overwrite previous results.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Optional. Number of threads. 3 or more is recommended. Default: 2', default=2)
parser.add_argument('-v', '--version', help='Print version.', action='store_true')
args=parser.parse_args()
args.version=version


# start
import init
init.init_jointcall(args, version)


# logging
import log
if args.merge_mei is True:
    args.logfilename='for_debug_jointcall_ins.log'
else:
    args.logfilename='for_debug_jointcall_abs.log'
if os.path.exists(os.path.join(args.outdir, args.logfilename)) is True:
    os.remove(os.path.join(args.outdir, args.logfilename))
log.start_log(args)
log.logger.debug('Logging started.')


# initial check
import initial_check
log.logger.debug('This is make_joint_call.py version %s' % version)
print()
log.logger.info('Initial check started.')
initial_check.check_merge_vcf(args, sys.argv)


# set up
import setup
setup.setup_merge_vcf(args, init.base)
params=setup.params
args.fai=setup.fai_path
args.rep_headers_to_be_removed=setup.rep_headers_to_be_removed


# output file names
import utils
filenames=utils.empclass()

filenames.repdb           =os.path.join(args.outdir, 'repdb')
filenames.rep_unknown_fa  =os.path.join(args.outdir, 'rep_unknown.fa')
filenames.blast_tmp_res   =os.path.join(args.outdir, 'blastn_tmp.txt')
filenames.reshaped_rep    =os.path.join(args.outdir, 'reshaped_repbase.fa')

filenames.merged_vcf_ins  =os.path.join(args.outdir, '%s_MEI_jointcall.vcf' % args.cohort_name)
filenames.merged_vcf_abs  =os.path.join(args.outdir, '%s_MEA_jointcall.vcf' % args.cohort_name)


# 0. preprocess repbase file
if args.merge_absent_me is True:
    import reshape_rep, blastn
    print()
    log.logger.info('Preprocess started.')
    reshape_rep.reshape(args, params, filenames)
    for ext in ['nhr', 'nin', 'nog', 'nsd', 'nsi', 'nsq']:
        os.remove('%s.%s' % (filenames.repdb, ext))


# 1. joint calling
import merge_vcf_strict
if args.merge_mei is True:
    merge_vcf_strict.merge_vcf_ins(args, params, filenames)
elif args.merge_absent_me is True:
    merge_vcf_strict.merge_vcf_abs(args, params, filenames)
    os.remove(filenames.reshaped_rep)


# all finish!
utils.output_finish_comment_merge_vcf(args, filenames)
