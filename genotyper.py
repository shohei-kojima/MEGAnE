#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,logging

'''
python /home/kooojiii/results/2020/prog_develop/koji_mei_pipeline/genotyper.py -c NA12878.final.cram -fa /home/kooojiii/Documents/genomes/hg38/1kGP/GRCh38_full_analysis_set_plus_decoy_hla.fa -fai /home/kooojiii/Documents/genomes/hg38/1kGP/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai -ins_bed result_out/MEI_final_percentile.bed -abs_bed result_out/absent_MEs.bed -abs_3t_bed result_out/absent_MEs_transduction.bed  -p 4 -overwrite -cov 35
'''


# version
version='2020/07/06'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-b', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end BAM file.')
parser.add_argument('-c', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end CRAM file.')
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: GRCh38DH.fa')
parser.add_argument('-fai', metavar='str', type=str, help='Required. Specify fasta index of the reference genome. Example: GRCh38DH.fa.fai')
parser.add_argument('-ins_bed', metavar='str', type=str, help='Specify input bed file containing information of MEIs.')
parser.add_argument('-abs_bed', metavar='str', type=str, help='Specify input bed file containing information of absent MEs.')
parser.add_argument('-abs_3t_bed', metavar='str', type=str, help='Specify input bed file containing information of absent MEs.')
parser.add_argument('-cov', metavar='int', type=int, help='Optional. Specify coverage depth. Default: 30', default=30)
parser.add_argument('-sample_name', metavar='int', type=int, help='Optional. Specify sample name which will be labeled in the VCF output. If not specified, BAM/CRAM filename will be output.')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./genotype_out', default='./genotype_out')
parser.add_argument('-mainchr', metavar='str', type=str, help='Optional. Specify full path if you analyze non-human sample. Default: /path/to/prog/lib/human_main_chrs.txt')
parser.add_argument('-sex', metavar='str', type=str, help='Optional. Specify gender of the sample; male or female or unknown. Available only when human sample. Default: male', default='male')
parser.add_argument('-setting', metavar='str', type=str, help='Optional. Specify full path to the parameter setting file. Default: /path/to/prog/lib/parameter_settings.txt')
parser.add_argument('-overwrite', help='Optional. Specify if you overwrite previous results.', action='store_true')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Optional. Number of threads. 3 or more is recommended. Default: 1', default=1)
parser.add_argument('-v', '--version', help='Print version.', action='store_true')
args=parser.parse_args()
args.version=version
#args.readlen=150  # dummy, any read length is fine


# start
import init
init.init_geno(args, version)


# logging
import log
args.logfilename='for_debug_genotyping.log'
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


# set up
import setup
setup.setup_geno(args, init.base)
params=setup.params
args.main_chrs=setup.main_chrs
args.main_chrs_set=set(args.main_chrs)
args.fai=setup.fai_path
params.chrY=setup.chrY
params.female=setup.female
args.search_geno=False


# output file names
import utils
filenames=utils.empclass()

filenames.limited_b       =os.path.join(args.outdir, 'only_necessary.bam')
filenames.limited_c       =os.path.join(args.outdir, 'only_necessary.cram')
filenames.tmp_bam         =os.path.join(args.outdir, 'tmp.bam')

base=os.path.splitext(os.path.basename(args.ins_bed))[0]
filenames.depth_ins       =os.path.join(args.outdir, '%s_depth_ins.txt' % base)
filenames.out_spanning    =os.path.join(args.outdir, '%s_spanning_read_summary.txt.gz' % base)
filenames.disc_read_pdf   =os.path.join(args.outdir, '%s_discordant_read_num.pdf' % base)

filenames.debug_pdf1      =os.path.join(args.outdir, 'plot_out_%s_genotype_ins_for_debug.pdf' % base)
filenames.merged_pdf      =os.path.join(args.outdir, 'plot_out_%s_genotyping_insertions.pdf' % base)
filenames.ins_out_bed     =os.path.join(args.outdir, '%s_genotyped.bed' % base)
filenames.ins_out_vcf     =os.path.join(args.outdir, '%s_genotyped.vcf' % base)

filenames.depth_abs       =os.path.join(args.outdir, 'depth_abs.txt')
filenames.merged_pdf_abs  =os.path.join(args.outdir, 'plot_out_genotyping_absents.pdf')
if not args.abs_bed is None:
    base=os.path.splitext(os.path.basename(args.abs_bed))[0]
    filenames.abs_out_bed     =os.path.join(args.outdir, '%s_genotyped.bed' % base)
    filenames.abs_out_vcf     =os.path.join(args.outdir, '%s_genotyped.vcf' % base)


#  0. limit BAM/CRAM
import output_genotyped_vcf
import allele_count_ins
log.logger.info('Limit BAM/CRAM started.')
allele_count_ins.limit(args, params, filenames)
data=utils.empclass()


# 1. insertions
if not args.ins_bed is None:
    log.logger.info('Evidence search started, insertion.')

    allele_count_ins.evaluate_tsd_depth(args, params, filenames)
    data.cn_est_tsd_depth=allele_count_ins.cn_est_tsd_depth
    data.tsd_thresholds=allele_count_ins.tsd_thresholds
    data.del_thresholds=allele_count_ins.del_thresholds

    allele_count_ins.evaluate_spanning_read(args, params, filenames)
    data.cn_est_spanning=allele_count_ins.cn_est_spanning
    data.spanning_thresholds=allele_count_ins.spanning_thresholds

    allele_count_ins.evaluate_discordant(args, params, filenames)
    data.cn_est_disc=allele_count_ins.cn_est_disc
    data.disc_thresholds=allele_count_ins.disc_thresholds  # could be False

    # merge evidences; insertion
    import merge_allele_evidence_ins
    log.logger.info('Evidence merge started, insertion.')
    merge_allele_evidence_ins.plot_orig(args, params, filenames, data)
    merge_allele_evidence_ins.merge(args, params, filenames, data)
    data.merged_res=merge_allele_evidence_ins.merged_res
    merge_allele_evidence_ins.plot_merged(args, params, filenames, data)
    data.mei_filter=merge_allele_evidence_ins.mei_filter
#    merge_allele_evidence_ins.plot_custom(args, params, filenames, data)   # always commentout unless debug

    # output; insertion
    output_genotyped_vcf.output_ins_bed_vcf(args, params, filenames, data)
    log.logger.info('Did output VCF, insertion.')

    # delete unnecessary files
    if os.path.exists(filenames.tmp_bam) is True:
        os.remove(filenames.tmp_bam)
    

# 2. absent MEs
if not args.abs_bed is None:
    import allele_count_abs
    log.logger.info('Evidence search started, absent MEs.')
    data=utils.empclass()

    allele_count_abs.evaluate_spanning_read(args, params, filenames)
    data.cn_est_spanning=allele_count_abs.cn_est_spanning
    data.spanning_thresholds=allele_count_abs.spanning_thresholds

    allele_count_abs.evaluate_bp_depth(args, params, filenames)
    data.cn_est_depth=allele_count_abs.cn_est_depth
    data.abs_thresholds=allele_count_abs.abs_thresholds

    # merge evidences; absent
    import merge_allele_evidence_abs
    log.logger.info('Evidence merge started, absent MEs.')
    merge_allele_evidence_abs.merge(args, params, filenames, data)
    data.merged_res=merge_allele_evidence_abs.merged_res
    merge_allele_evidence_abs.plot_merged(args, params, filenames, data)

    # output; absent
    output_genotyped_vcf.output_abs_bed_vcf(args, params, filenames, data)
    log.logger.info('Did output VCF, absent MEs.')

    # delete unnecessary files
    if args.keep is False:
        if os.path.exists(filenames.depth) is True:
            os.remove(filenames.depth)
    if os.path.exists(filenames.tmp_bam) is True:
        os.remove(filenames.tmp_bam)


# finish!
log.logger.info('Genotyping finished!')
