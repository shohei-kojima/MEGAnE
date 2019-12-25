#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse

'''
time python main.py -overwrite -b test_data/test.bam -rep test_data/humrepsub.fa -p 2
'''

# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-b', metavar='str', type=str, help='Required. Specify input mapped paired-end BAM file.')  # , required=True
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: hg38.fa')
parser.add_argument('-fai', metavar='str', type=str, help='Required. Specify fasta index of the reference genome. Example: hg38.fa.fai')
parser.add_argument('-rep', metavar='str', type=str, help='Required. Specify RepBase file used for repeatmasking. Example: humrep.ref')
parser.add_argument('-repout', metavar='str', type=str, help='Required. Specify RepeatMasker output. Example: hg38.fa.out')
parser.add_argument('-cov', metavar='int', type=int, help='Optional. Specify coverage depth. Default: 30', default=30)
parser.add_argument('-readlen', metavar='int', type=int, help='Optional. Specify read length. Default: 150', default=150)
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./result_out', default='./result_out')
parser.add_argument('-mainchr', metavar='str', type=str, help='Optional. Specify full path if you analyze non-human sample. Default: /path/to/prog/lib/human_main_chrs.txt')
parser.add_argument('-gender', metavar='str', type=str, help='Optional. Specify gender of the sample; male or female or unknown. Available only when human sample. Default: unknown', default='unknown')
parser.add_argument('-setting', metavar='str', type=str, help='Optional. Specify full path to the parameter setting file. Default: /path/to/prog/lib/parameter_settings.txt')
parser.add_argument('-only_ins', help='Optional. Specify if you only analyze non-reference MEI insertions.', action='store_true')
parser.add_argument('-only_abs', help='Optional. Specify if you only analyze absence of reference MEI.', action='store_true')
parser.add_argument('-overwrite', help='Optional. Specify to overwrite when outdir and result files already exist.', action='store_true')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Number of thread to use for blastn. Default: 1', default=1)
args=parser.parse_args()


# start
import init
init.init()
base=init.base


# initial check
import initial_check
initial_check.check(args)


# set up working dir
import setup
setup.setup(args, base)
main_chrs=setup.main_chrs
main_chrs_set=set(main_chrs)
params=setup.params


# output file names
reshaped_rep    =os.path.join(args.outdir, 'reshaped_repbase.fa')
rep_slide_file  =os.path.join(args.outdir, 'reshaped_repbase_slide.fa')
blast0_res      =os.path.join(args.outdir, 'blastn_rep_slide.txt')
similar_rep_list=os.path.join(args.outdir, 'similar_rep_list.txt')

overhang_fa     =os.path.join(args.outdir, 'overhang.fa')
overhang_pA     =os.path.join(args.outdir, 'overhang_pA.txt')
distant_txt     =os.path.join(args.outdir, 'distant.txt')
unmapped_fa     =os.path.join(args.outdir, 'unmapped.fa')
blast1_res      =os.path.join(args.outdir, 'blastn_result_overhang_to_rep.txt')
overhang_MEI    =os.path.join(args.outdir, 'overhang_to_MEI_list.txt')
additional_pA   =os.path.join(args.outdir, 'overhang_additional_pA.txt')
blast2_res      =os.path.join(args.outdir, 'blastn_result_unmapped_to_rep.txt')
unmapped_hit_fa =os.path.join(args.outdir, 'unmapped_ref_side.fa')
blast3_res      =os.path.join(args.outdir, 'blastn_result_unmapped_mei_to_ref.txt')
unmapped_MEI    =os.path.join(args.outdir, 'unmapped_to_MEI_list.txt')
hybrid          =os.path.join(args.outdir, 'hybrid_reads.txt')


# preprocess repbase file
import reshape_rep, blastn, utils
reshape_rep.reshape(args, reshaped_rep)

blastn.makeblastdb(reshaped_rep)
args.repdb=blastn.dbpath

reshape_rep.slide_rep_file(args, params, reshaped_rep, rep_slide_file)
#blastn.blastn(args, params, rep_slide_file, args.repdb, blast0_res)  # params, q_path, db_path, outfpath
reshape_rep.parse_slide_rep_blastn_res(args, blast0_res, similar_rep_list)

#os.remove(blast0_res)
#os.remove(rep_slide_file)
exit()


### insertion finder ###
import parse_blastn_result, find_additional_pA, extract_discordant, extract_discordant_c


# 1. process unmapped overhangs
# extract discordant reads from bam
outfiles=[overhang_fa, overhang_pA, distant_txt, unmapped_fa]
if args.p >= 2:
    from multiprocessing import Pool
    count,interval=extract_discordant.flagstat(args)
    def extract_discordant_exe(n):
        extract_discordant_c.main(args, params, main_chrs_set, outfiles, n, count, interval)
    with Pool(args.p) as p:
        p.map(extract_discordant_exe, range(args.p))
    extract_discordant.concat(args, outfiles)
else:
    extract_discordant_c.main(args, params, main_chrs_set, outfiles, None, None, None)
exit()
# blastn overhangs to repeat seqs
blastn.blastn(args, params, overhang_fa, args.repdb, blast1_res)
# find overhangs originating from MEIs
parse_blastn_result.parse(params, blast1_res, overhang_MEI)
# find additional pA overhang from the rest of the overhangs
find_additional_pA.find(params, blast1_res, overhang_fa, additional_pA)
# gzip or delete unnecessary files
utils.gzip_or_del(args, params, blast1_res)

# 2. process unmapped reads
# blastn unmapped reads to repeat seqs
blastn.blastn(args, params, unmapped_fa, args.repdb, blast2_res)
# find unmapped reads mapped to MEIs and output fasta
parse_blastn_result.unmapped_to_fa(params, unmapped_fa, blast2_res, unmapped_hit_fa)
# gzip or delete unnecessary files
utils.gzip_or_del(args, params, blast2_res)
# blastn unmapped reads to repeat seqs
blastn.blastn(args, params, unmapped_hit_fa, args.fadb, blast3_res)
# find unmapped reads mapped both to MEIs and ref genome
parse_blastn_result.find_chimeric_unmapped(params, main_chrs_set, blast3_res, unmapped_MEI)
# gzip or delete unnecessary files
utils.gzip_or_del(args, params, blast3_res)
utils.gzip_or_del(args, params, unmapped_hit_fa)

# 3. process hybrid reads
import process_distant_read
process_distant_read.process_reads(args, params, distant_txt, hybrid)  # need change, repbase

# 4. merge all results, identify MEI

