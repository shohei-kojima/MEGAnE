#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse

'''
time python main.py -overwrite -b test_data/test.bam -rep test_data/humrepsub.fa -repout /home/kooojiii/Documents/genomes/hg38/ucsc/masked_using_RepBase24.01_humrep_humsub/hg38.fa.out -p 2
'''

# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-b', metavar='str', type=str, help='Required. Specify input mapped paired-end BAM file.')  # , required=True
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: hg38.fa')
parser.add_argument('-fai', metavar='str', type=str, help='Required. Specify fasta index of the reference genome. Example: hg38.fa.fai')
parser.add_argument('-rep', metavar='str', type=str, help='Required. Specify RepBase file used for repeatmasking. Example: humrep.ref')
parser.add_argument('-repout', metavar='str', type=str, help='Required. Specify RepeatMasker output. Must be masked using the input RepBase file. Example: hg38.fa.out')
parser.add_argument('-cov', metavar='int', type=int, help='Optional. Specify coverage depth. Default: 30', default=30)
parser.add_argument('-readlen', metavar='int', type=int, help='Optional. Specify read length. Default: 150', default=150)
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./result_out', default='./result_out')
parser.add_argument('-mainchr', metavar='str', type=str, help='Optional. Specify full path if you analyze non-human sample. Default: /path/to/prog/lib/human_main_chrs.txt')
parser.add_argument('-gender', metavar='str', type=str, help='Optional. Specify gender of the sample; male or female or unknown. Available only when human sample. Default: unknown', default='unknown')
parser.add_argument('-setting', metavar='str', type=str, help='Optional. Specify full path to the parameter setting file. Default: /path/to/prog/lib/parameter_settings.txt')
parser.add_argument('-repremove', metavar='str', type=str, help='Optional. Specify full path to a file containing the names of non-ME repeat class. Default: /path/to/prog/lib/non_ME_rep_headers.txt')
parser.add_argument('-pA_ME', metavar='str', type=str, help='Optional. Specify full path to a file containing repat class with polyA tail. Default: /path/to/prog/lib/ME_with_polyA_tail.txt')
parser.add_argument('-only_ins', help='Optional. Specify if you only analyze non-reference MEI insertions.', action='store_true')
parser.add_argument('-only_abs', help='Optional. Specify if you only analyze absence of reference MEI.', action='store_true')
parser.add_argument('-overwrite', help='Optional. Specify if you overwrite previous results.', action='store_true')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Optional. Number of threads. 3 or more is recommended. Default: 1', default=1)
args=parser.parse_args()


# start
import init
init.init()
base=init.base


# initial check
import initial_check
initial_check.check(args)


# set up
import setup
setup.setup(args, base)
params=setup.params
args.main_chrs=setup.main_chrs
args.main_chrs_set=set(args.main_chrs)
args.rep_headers_to_be_removed=setup.rep_headers_to_be_removed
args.rep_with_pA=setup.rep_with_pA


# output file names
import utils
filenames=utils.empclass()

filenames.repdb           =os.path.join(args.outdir, 'repdb')
filenames.repout_bed      =os.path.join(args.outdir, 'repout.bed')
filenames.reshaped_rep    =os.path.join(args.outdir, 'reshaped_repbase.fa')
filenames.rep_slide_file  =os.path.join(args.outdir, 'reshaped_repbase_slide.fa')
filenames.blast0_res      =os.path.join(args.outdir, 'blastn_rep_slide.txt')
filenames.similar_rep_list=os.path.join(args.outdir, 'similar_rep_list.txt')

filenames.overhang_fa     =os.path.join(args.outdir, 'overhang.fa')
filenames.overhang_pA     =os.path.join(args.outdir, 'overhang_pA.txt')
filenames.distant_txt     =os.path.join(args.outdir, 'distant.txt')
filenames.unmapped_fa     =os.path.join(args.outdir, 'unmapped.fa')
filenames.mapped_fa       =os.path.join(args.outdir, 'mapped.fa')
filenames.abs_txt         =os.path.join(args.outdir, 'absent.txt')

filenames.blast1_res      =os.path.join(args.outdir, 'blastn_result_overhang_to_rep.txt')
filenames.overhang_MEI    =os.path.join(args.outdir, 'overhang_to_MEI_list.txt')
filenames.additional_pA   =os.path.join(args.outdir, 'overhang_additional_pA.txt')
filenames.blast2_res      =os.path.join(args.outdir, 'blastn_result_unmapped_to_rep.txt')
filenames.unmapped_hit_fa =os.path.join(args.outdir, 'unmapped_ref_side.fa')
filenames.blast3_res      =os.path.join(args.outdir, 'blastn_result_unmapped_mei_to_ref.txt')
filenames.unmapped_MEI    =os.path.join(args.outdir, 'unmapped_to_MEI_list.txt')
filenames.hybrid          =os.path.join(args.outdir, 'hybrid_reads.txt')

filenames.breakpoint_pairs=os.path.join(args.outdir, 'breakpoint_pairs.txt')
filenames.breakpoint_info =os.path.join(args.outdir, 'breakpoint_pairs_info.txt')
filenames.breakpoint_clean=os.path.join(args.outdir, 'breakpoint_pairs_clean.txt')

filenames.mapped_fa_select=os.path.join(args.outdir, 'mapped_selected.fa')
filenames.blast4_res      =os.path.join(args.outdir, 'blastn_result_mapped_ref.txt')
filenames.bp_pair_single  =os.path.join(args.outdir, 'breakpoint_pairs_in_TE_singletons.txt')
filenames.bp_info_single  =os.path.join(args.outdir, 'breakpoint_pairs_in_TE_singletons_info.txt')
filenames.bp_clean_single =os.path.join(args.outdir, 'breakpoint_pairs_in_TE_singletons_clean.txt')

filenames.bp_merged       =os.path.join(args.outdir, 'breakpoint_pairs_pooled_merged.txt')
filenames.hybrid_master   =os.path.join(args.outdir, 'hybrid_reads_master.txt')  # one of the final outputs
filenames.bp_merged_all   =os.path.join(args.outdir, 'breakpoint_pairs_pooled_all.txt')
filenames.bp_merged_filt  =os.path.join(args.outdir, 'breakpoint_pairs_pooled_filtered.txt')
filenames.bp_merged_group =os.path.join(args.outdir, 'breakpoint_pairs_pooled_grouped.txt')


# preprocess repbase file
import reshape_rep, blastn
reshape_rep.reshape(args, filenames)
#blastn.makeblastdb(filenames.reshaped_rep, filenames.repdb)
#reshape_rep.slide_rep_file(args, params, filenames)
#blastn.blastn(args, params, filenames.rep_slide_file, filenames.repdb, filenames.blast0_res)  # params, q_path, db_path, outfpath
reshape_rep.parse_slide_rep_blastn_res(args, filenames)
reshape_rep.reshape_repout_to_bed(args, filenames)

#os.remove(blast0_res)
#os.remove(rep_slide_file)


### insertion finder ###
import parse_blastn_result, find_additional_pA, extract_discordant, extract_discordant_c


# 1. process unmapped overhangs
if args.p >= 2:
    from multiprocessing import Pool
    count,interval=extract_discordant.flagstat(args)
    def extract_discordant_exe(n):
        extract_discordant.main(args, params, filenames, n, count, interval)  ### c
    with Pool(args.p) as p:
        p.map(extract_discordant_exe, range(args.p))
    extract_discordant.concat(args, filenames)
else:
    extract_discordant_c.main(args, params, filenames, None, None, None)
exit()
blastn.blastn(args, params, filenames.overhang_fa, filenames.repdb, filenames.blast1_res)
parse_blastn_result.parse(params, filenames.blast1_res, filenames.overhang_MEI)
find_additional_pA.find(params, filenames.blast1_res, filenames.overhang_fa, filenames.additional_pA)

utils.gzip_or_del(args, params, filenames.blast1_res)

# 2. process unmapped reads
blastn.blastn(args, params, filenames.unmapped_fa, filenames.repdb, filenames.blast2_res)
parse_blastn_result.unmapped_to_fa(params, filenames.unmapped_fa, filenames.blast2_res, filenames.unmapped_hit_fa)
utils.gzip_or_del(args, params, filenames.blast2_res)
blastn.blastn(args, params, filenames.unmapped_hit_fa, args.fadb, filenames.blast3_res)
parse_blastn_result.find_chimeric_unmapped(args, params, filenames.blast3_res, filenames.unmapped_MEI)

utils.gzip_or_del(args, params, filenames.blast3_res)
utils.gzip_or_del(args, params, filenames.unmapped_hit_fa)

# 3. process hybrid reads
import process_distant_read
process_distant_read.process_reads(args, params, filenames, filenames.distant_txt, filenames.hybrid)  # need change, repbase

# 4. merge all results, identify MEI outside of similar MEs
import pair_breakpoints
pair_breakpoints.pairing(args, params, filenames)
pair_breakpoints.add_TE_subclass(filenames, filenames.breakpoint_pairs, filenames.breakpoint_info)
pair_breakpoints.remove_cand_inside_TE(args, params, filenames)

# 5. identify MEI nested in similar MEs
import process_mapped_seq
process_mapped_seq.retrieve_mapped_seq(filenames)
process_mapped_seq.blastn_for_mapped(args, params, filenames.mapped_fa_select, args.fadb, filenames.blast4_res)
process_mapped_seq.pairing(params, filenames)
pair_breakpoints.add_TE_subclass(filenames, filenames.bp_pair_single, filenames.bp_info_single)
pair_breakpoints.remove_redundant_pairs(filenames)

# 6. pool results and filter candidates
import pool_results
pool_results.merge_breakpoints(filenames)
pool_results.add_hybrid(params, filenames)
import filter_candidates
filter_candidates.filter(args, params, filenames)
filter_candidates.grouping(filenames)
