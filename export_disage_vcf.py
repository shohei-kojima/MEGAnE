#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,logging

'''
time python main.py -overwrite -b test_data/NA12878.chr22.bam -fa /home/kooojiii/Documents/genomes/hg38/hg38.fa -fadb /home/kooojiii/Documents/genomes/hg38/hg38 -rep test_data/humrepsub.fa -repout /home/kooojiii/Documents/genomes/hg38/ucsc/masked_using_RepBase24.01_humrep_humsub/hg38.fa.out -cov 35 -p 4
time python main.py -overwrite -b ../191216_2/NA12878.final.bam -fa /home/kooojiii/Documents/genomes/hg38/hg38.fa -fadb /home/kooojiii/Documents/genomes/hg38/hg38 -rep test_data/humrepsub.fa -repout /home/kooojiii/Documents/genomes/hg38/ucsc/masked_using_RepBase24.01_humrep_humsub/hg38.fa.out -cov 35 -p 4 -only_abs

python /home/kooojiii/results/2020/prog_develop/koji_mei_pipeline/export_disage_vcf.py -b NA12878.final.cram -fa /home/kooojiii/Documents/genomes/hg38/ucsc/hg38.fa -fadb /home/kooojiii/Documents/genomes/hg38/ucsc/hg38 -rep /home/kooojiii/results/2020/prog_develop/koji_mei_pipeline/lib/humrepsub.fa -mainchr /home/kooojiii/results/2020/prog_develop/koji_mei_pipeline/lib/GRCh38DH_primary_plus_alt_ucsc_style.txt -p 4
'''


# version
version='2020/05/05'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-b', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end BAM file.')  # , required=True
parser.add_argument('-c', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end CRAM file.')  # , required=True
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: hg38.fa')
parser.add_argument('-fai', metavar='str', type=str, help='Required. Specify fasta index of the reference genome. Example: hg38.fa.fai')
parser.add_argument('-rep', metavar='str', type=str, help='Required. Specify RepBase file used for repeatmasking. Example: humrep.ref')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./result_out', default='./result_out')
parser.add_argument('-mainchr', metavar='str', type=str, help='Optional. Specify full path if you analyze non-human sample. Default: /path/to/prog/lib/human_main_chrs.txt')
parser.add_argument('-monoallelic', help='Optional. Specify if you use monoalellic sample, such as mouse strains or HAP1 cells.', action='store_true')
parser.add_argument('-gender', metavar='str', type=str, help='Optional. Specify gender of the sample; male or female or unknown. Available only when human sample. Default: unknown', default='unknown')
parser.add_argument('-setting', metavar='str', type=str, help='Optional. Specify full path to the parameter setting file. Default: /path/to/prog/lib/parameter_settings.txt')
parser.add_argument('-repremove', metavar='str', type=str, help='Optional. Specify full path to a file containing the names of non-ME repeat class. Default: /path/to/prog/lib/non_ME_rep_headers.txt')
parser.add_argument('-pA_ME', metavar='str', type=str, help='Optional. Specify full path to a file containing repat class with polyA tail. Default: /path/to/prog/lib/ME_with_polyA_tail.txt')
parser.add_argument('-overwrite', help='Optional. Specify if you overwrite previous results.', action='store_true')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Optional. Number of threads. 3 or more is recommended. Default: 1', default=1)
parser.add_argument('-v', '--version', help='Print version.', action='store_true')
args=parser.parse_args()
args.version=version


# start
import init
init.init(args, version)


# logging
import log
args.logfilename='export_dosage_vcf.log'
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
setup.setup(args, init.base)
params=setup.params
args.main_chrs=setup.main_chrs
args.main_chrs_set=set(args.main_chrs)


# output file names
import utils
filenames=utils.empclass()

filenames.bp_final_g      =os.path.join(args.outdir, 'MEI_final_gaussian.bed')
filenames.bp_final_p      =os.path.join(args.outdir, 'MEI_final_percentile.bed')
filenames.bp_final_u      =os.path.join(args.outdir, 'MEI_final_user.bed')

filenames.vcf_g           =os.path.join(args.outdir, 'MEI_final_gaussian.vcf')
filenames.vcf_p           =os.path.join(args.outdir, 'MEI_final_percentile.vcf')
filenames.vcf_u           =os.path.join(args.outdir, 'MEI_final_user.vcf')

filenames.disage_pdf_g    =os.path.join(args.outdir, 'plot_gene_dosage_gaussian.pdf')
filenames.disage_pdf_p    =os.path.join(args.outdir, 'plot_gene_dosage_percentile.pdf')
filenames.disage_pdf_u    =os.path.join(args.outdir, 'plot_gene_dosage_user.pdf')


# output vcf
import infer_dosage
if os.path.exists(filenames.bp_final_g) is True:
    infer_dosage.calc_dosage(args, params, filenames, filenames.bp_final_g, filenames.disage_pdf_g, filenames.vcf_g)
    

# output comments
log.logger.info('Output vcf finished!')

