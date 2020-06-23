#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,logging

'''
python /home/kooojiii/results/2020/prog_develop/koji_mei_pipeline/genotyper.py -c NA12878.final.cram -fa /home/kooojiii/Documents/genomes/hg38/1kGP/GRCh38_full_analysis_set_plus_decoy_hla.fa -fai /home/kooojiii/Documents/genomes/hg38/1kGP/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai -ins_bed result_out/MEI_final_percentile.bed -p 4 -overwrite
'''


# version
version='2020/06/21'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-b', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end BAM file.')
parser.add_argument('-c', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end CRAM file.')
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: hg38.fa')
parser.add_argument('-fai', metavar='str', type=str, help='Required. Specify fasta index of the reference genome. Example: hg38.fa.fai')
parser.add_argument('-ins_bed', metavar='str', type=str, help='Specify input bed file containing information of MEIs.')
parser.add_argument('-abs_bed', metavar='str', type=str, help='Specify input bed file containing information of absent MEs.')
parser.add_argument('-cov', metavar='int', type=int, help='Optional. Specify coverage depth. Default: 30', default=30)
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./genotype_out', default='./genotype_out')
parser.add_argument('-gender', metavar='str', type=str, help='Optional. Specify gender of the sample; male or female or unknown. Available only when human sample. Default: unknown', default='unknown')
parser.add_argument('-setting', metavar='str', type=str, help='Optional. Specify full path to the parameter setting file. Default: /path/to/prog/lib/parameter_settings.txt')
parser.add_argument('-overwrite', help='Optional. Specify if you overwrite previous results.', action='store_true')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Optional. Number of threads. 3 or more is recommended. Default: 1', default=1)
parser.add_argument('-v', '--version', help='Print version.', action='store_true')
args=parser.parse_args()
args.version=version
args.readlen=150  # any read length is fine


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
#log.logger.debug('This is version %s' % version)
#print()
#log.logger.info('Initial check started.')
#initial_check.check(args, sys.argv)


# set up
import setup
setup.setup_geno(args, init.base)
params=setup.params


# output file names
import utils
filenames=utils.empclass()

filenames.limited_b       =os.path.join(args.outdir, 'only_necessary.bam')
filenames.limited_c       =os.path.join(args.outdir, 'only_necessary.cram')
filenames.limited_f2       =os.path.join(args.outdir, 'only_necessary_f')
filenames.depth_ins       =os.path.join(args.outdir, 'depth_ins.txt')

filenames.tmp_bam         =os.path.join(args.outdir, 'tmp.bam')
filenames.out_spanning    =os.path.join(args.outdir, 'spanning_read_summary.txt.gz')

filenames.disc_read_pdf   =os.path.join(args.outdir, 'discordant_read_num.pdf')


# 0. limit BAM/CRAM
import allele_count_ins
log.logger.info('Limit BAM/CRAM started.')
#allele_count_ins.limit(args, params, filenames)
allele_count_ins.evaluate_tsd_depth(args, params, filenames)
cn_est_tsd_depth=allele_count_ins.cn_est_tsd_depth
allele_count_ins.evaluate_spanning_read(args, params, filenames)
cn_est_spanning=allele_count_ins.cn_est_spanning
allele_count_ins.evaluate_discordant(args, params, filenames)
cn_est_disc=allele_count_ins.cn_est_disc


# output comments
log.logger.info('Output vcf finished!')
