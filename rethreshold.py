#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,logging,gzip

'''
python rethreshold.py -b test_data/NA12878.chr22.bam -fa /home/kooojiii/Documents/genomes/hg38/hg38.fa -fadb /home/kooojiii/Documents/genomes/hg38/hg38 -rep test_data/humrepsub.fa -repout /home/kooojiii/Documents/genomes/hg38/ucsc/masked_using_RepBase24.01_humrep_humsub/hg38.fa.out -cov 35 -threshold 10
'''


# version
version='2021/02/12'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-b', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end BAM file.')  # , required=True
parser.add_argument('-c', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end CRAM file.')  # , required=True
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: GRCh38DH.fa')
parser.add_argument('-fadb', metavar='str', type=str, help='Required. Specify blastdb of reference genome. Example: GRCh38DH.fa.db')
parser.add_argument('-rep', metavar='str', type=str, help='Required. Specify RepBase file used for repeatmasking. Example: humrep.ref')
parser.add_argument('-repout', metavar='str', type=str, help='Required. Specify RepeatMasker output. Must be masked using the input RepBase file. Example: GRCh38DH.fa.out')
parser.add_argument('-cov', metavar='int or auto', type=str, help='Optional. Specify mapping depth. If "auto" was specified, it estimate autosome depth. Default: auto', default='auto')
parser.add_argument('-readlen', metavar='int or auto', type=str, help='Optional. Specify read length. If "auto" was specified, it estimate read length. Default: auto', default='auto')
parser.add_argument('-threshold', metavar='int', type=int, help='Optional. Specify user-defined threshold.')
parser.add_argument('-sample_name', metavar='str', type=str, help='Optional. Specify sample name which will be labeled in the VCF output. If not specified, BAM/CRAM filename will be output.')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./result_out', default='./result_out')
parser.add_argument('-mainchr', metavar='str', type=str, help='Optional. Specify full path if you analyze non-human sample. Default: /path/to/prog/lib/human_main_chrs.txt')
parser.add_argument('-monoallelic', help='Optional. Specify if you use monoalellic sample, such as inbread mouse strains or HAP1 cells.', action='store_true')
parser.add_argument('-sex', metavar='str', type=str, help='Optional. Available only when input has chrX and chrY (e.g. human, mouse). Specify "female" or "male" or "auto". If "auto" is specified, it automatically estimate the sex. If your input foes not have chrX and chrY, please remove this option. Default: auto', default='auto')
parser.add_argument('-verylowdep', help='Optional. Specify if you use parameter settings for low depth (generally less than 10x).', action='store_true')
parser.add_argument('-lowdep', help='Optional. Specify if you use parameter settings for low depth (generally less than 15x).', action='store_true')
parser.add_argument('-setting', metavar='str', type=str, help='Optional. Specify full path to the parameter setting file. Default: /path/to/prog/lib/parameter_settings.txt')
parser.add_argument('-repremove', metavar='str', type=str, help='Optional. Specify full path to a file containing the names of non-ME repeat class. Default: /path/to/prog/lib/human_non_ME_rep_headers.txt')
parser.add_argument('-pA_ME', metavar='str', type=str, help='Optional. Specify full path to a file containing repat class with polyA tail. Default: /path/to/prog/lib/ME_with_polyA_tail.txt')
parser.add_argument('-only_ins', help='Optional. Specify if you only analyze non-reference MEI insertions.', action='store_true')
parser.add_argument('-only_abs', help='Optional. Specify if you only analyze absence of reference MEI.', action='store_true')
parser.add_argument('-only_geno', help='Optional. Specify if you only genotype using previously analyzed data.', action='store_true')
parser.add_argument('-unsorted', help='Optional. Specify if an input BAM/CRAM is not position sorted.', action='store_true')
parser.add_argument('-pybedtools_tmp', metavar='str', type=str, help='Optional. Specify directory for temporary bedtools files, e.g. /dev/shm')
parser.add_argument('-do_not_overwrite', help='Optional. Specify if you do NOT overwrite previous results.', action='store_true')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-no_pdf', help='Optional. Specify if you do not want to output pdf summary files.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Optional. Number of threads. 3 or more is recommended. Default: 2', default=2)
parser.add_argument('-v', '--version', action='version', version='MEGAnE %s %s' % (os.path.basename(__file__), version))
parser.add_argument('-only_geno_precall', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('-skip_unmapped', action='store_true', help=argparse.SUPPRESS)
args=parser.parse_args()
args.version=version


# always overwrite
args.overwrite=True


# start
import init
init.init(args, version)


# logging
import log
args.logfilename='rethreshold.log'
if os.path.exists(os.path.join(args.outdir, args.logfilename)) is True:
    os.remove(os.path.join(args.outdir, args.logfilename))
log.start_log(args)
log.logger.debug('Logging started.')


# initial check
import initial_check
log.logger.debug('This is version %s' % version)
print()
log.logger.info('Initial check started.')
initial_check.check(args, sys.argv)
if args.p >= 2:
    log.logger.info('You specified %d cores (-p option), but rethresholding will use single core. Proceed anyway.' % args.p)


# set up
import setup,auto_setting
args.auto=auto_setting.init()
if args.readlen in args.auto:
    auto_setting.estimate_readlen(args)
setup.setup(args, init.base)
params=setup.params
args.main_chrs=setup.main_chrs
args.main_chrs_set=set(args.main_chrs)
args.rep_headers_to_be_removed=setup.rep_headers_to_be_removed
args.rep_with_pA=setup.rep_with_pA
args.fai=setup.fai_path
params.chrX=setup.chrX
params.chrY=setup.chrY
params.female=setup.female
params.male=setup.male
do_ins=False if args.only_abs is True else True
do_abs=False if args.only_ins is True else True
if args.cov in args.auto or args.sex in args.auto:
    auto_setting.estimate_depth_sex(args, params, args.auto)


# output file names
import utils
filenames=utils.empclass()

filenames.breakpoint_pairs=os.path.join(args.outdir, 'breakpoint_pairs.txt')
filenames.bp_merged_all   =os.path.join(args.outdir, 'breakpoint_pairs_pooled_all.txt')
filenames.overhang_MEI    =os.path.join(args.outdir, 'overhang_to_MEI_list.txt')

filenames.bp_merged_filt_g=os.path.join(args.outdir, 'breakpoint_pairs_pooled_filtered_gaussian_ret.txt')
filenames.bp_merged_filt_p=os.path.join(args.outdir, 'breakpoint_pairs_pooled_filtered_percentile_ret.txt')
filenames.bp_merged_filt_f=os.path.join(args.outdir, 'breakpoint_pairs_pooled_filtered_failed_ret.txt')
filenames.bp_merged_filt_u=os.path.join(args.outdir, 'breakpoint_pairs_pooled_filtered_user_ret.txt')
filenames.bp_merged_groupg=os.path.join(args.outdir, 'breakpoint_pairs_pooled_grouped_gaussian_ret.txt')
filenames.bp_merged_groupp=os.path.join(args.outdir, 'breakpoint_pairs_pooled_grouped_percentile_ret.txt')
filenames.bp_merged_groupf=os.path.join(args.outdir, 'breakpoint_pairs_pooled_grouped_failed_ret.txt')
filenames.bp_merged_groupu=os.path.join(args.outdir, 'breakpoint_pairs_pooled_grouped_user_ret.txt')
filenames.gaussian_plot   =os.path.join(args.outdir, 'plot_gaussian_fitting_ret.pdf')
filenames.abs_res         =os.path.join(args.outdir, 'absent_MEs.bed')

filenames.bp_final_g      =os.path.join(args.outdir, 'MEI_final_gaussian_ret.bed')
filenames.bp_final_p      =os.path.join(args.outdir, 'MEI_final_percentile_ret.bed')
filenames.bp_final_f      =os.path.join(args.outdir, 'MEI_final_failed_ret.bed')
filenames.bp_final_u      =os.path.join(args.outdir, 'MEI_final_user_ret.bed')


# unzip
if os.path.exists(filenames.breakpoint_pairs +'.gz') is True:
    utils.gzip_d(filenames.breakpoint_pairs +'.gz')
elif os.path.exists(filenames.breakpoint_pairs) is True:
    pass
else:
    log.logger.error('%s.gz does not exist.' % filenames.breakpoint_pairs)
    exit(1)

if os.path.exists(filenames.bp_merged_all +'.gz') is True:
    utils.gzip_d(filenames.bp_merged_all +'.gz')
elif os.path.exists(filenames.bp_merged_all) is True:
    pass
else:
    log.logger.error('%s.gz does not exist.' % filenames.bp_merged_all)
    exit(1)

if os.path.exists(filenames.overhang_MEI +'.gz') is True:
    utils.gzip_d(filenames.overhang_MEI +'.gz')
elif os.path.exists(filenames.overhang_MEI) is True:
    pass
else:
    log.logger.error('%s.gz does not exist.' % filenames.overhang_MEI)
    exit(1)


# main rethresholding
import filter_candidates
filter_candidates.filter(args, params, filenames)
args.gaussian_executed=filter_candidates.gaussian_executed
filter_candidates.grouping(args, filenames)
log.logger.info('%s ME insertion candidates found.' % filter_candidates.ins_ns)

import after_processing
after_processing.grouped_mei_to_bed(args, params, filenames)

utils.gzip_file(params, filenames.breakpoint_pairs)
utils.gzip_file(params, filenames.bp_merged_all)
utils.gzip_file(params, filenames.overhang_MEI)

# monoallelic, output VCF
if args.monoallelic is True:
    import output_genotyped_vcf_mono
    for f in [filenames.bp_final_g, filenames.bp_final_p, filenames.bp_final_f, filenames.bp_final_u]:
        if os.path.exists(f) is True:
            args.ins_bed=f
            log.logger.info('Will skip genotyping, insertion: %s' % args.ins_bed)
            base=os.path.splitext(os.path.basename(args.ins_bed))[0]
            filenames.ins_out_bed   =os.path.join(args.outdir, '%s_genotyped.bed' % base)
            filenames.ins_out_vcf   =os.path.join(args.outdir, '%s_genotyped.vcf' % base)
            output_genotyped_vcf_mono.output_ins_bed_vcf(args, params, filenames)
            log.logger.info('Did output VCF, insertion: %s' % args.ins_bed)
    if do_abs is True:
        log.logger.info('Will skip genotyping, absent MEs.')
        args.abs_bed=filenames.abs_res
        args.abs_3t_bed=None
        base=os.path.splitext(os.path.basename(args.abs_bed))[0]
        filenames.abs_out_bed     =os.path.join(args.outdir, '%s_genotyped.bed' % base)
        filenames.abs_out_vcf     =os.path.join(args.outdir, '%s_genotyped.vcf' % base)
        output_genotyped_vcf_mono.output_abs_bed_vcf(args, params, filenames)
        log.logger.info('Did output VCF, absent MEs.')
    # cleanup
    if not args.pybedtools_tmp == args.outdir:
        shutil.rmtree(args.pybedtools_tmp)
    # all finish!
    utils.output_finish_comment(args, do_ins, do_abs, filenames)
    exit(0)

log.logger.info('Re-thresholding finished!\n')
