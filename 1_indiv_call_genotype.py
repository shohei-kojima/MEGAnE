#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,shutil,logging


# version
version='v0.1.1 2021/02/15'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-b', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end BAM file.')  # , required=True
parser.add_argument('-c', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end CRAM file.')  # , required=True
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. E.g. GRCh38DH.fa')
parser.add_argument('-fadb', metavar='str', type=str, help='Required. Specify blastdb of reference genome. E.g. GRCh38DH.fa.db')
parser.add_argument('-rep', metavar='str', type=str, help='Required. Specify Repeat consensus file. E.g. Dfam_custom.ref')
parser.add_argument('-repout', metavar='str', type=str, help='Required. Specify RepeatMasker output. Must be masked using the input RepBase file. E.g. GRCh38DH.fa.out')
parser.add_argument('-repremove', metavar='str', type=str, help='Required. Specify full path to a file containing the names of non-ME repeat class. E.g. non_ME_rep.txt')
parser.add_argument('-pA_ME', metavar='str', type=str, help='Required. Specify full path to a file containing repat class with polyA tail. E.g. ME_with_pA.txt')
parser.add_argument('-mainchr', metavar='str', type=str, help='Required. Specify full path if you analyze non-human sample. E.g. main_chr_list.txt')
parser.add_argument('-cov', metavar='int or auto', type=str, help='Optional, but specify whenever possible. Specify mapping depth. If "auto" was specified, it estimate autosome depth. "auto" is only available for human samples. Default: 30', default='30')
parser.add_argument('-readlen', metavar='int or auto', type=str, help='Optional. Specify read length. If "auto" was specified, it estimate read length. Default: auto', default='auto')
parser.add_argument('-sex', metavar='str', type=str, help='Optional, but specify whenever possible. Specify "female" or "male" or "auto" or "unknown". "auto" is only available for human samples. If "auto" is specified, it automatically estimate the sex. Default: auto for human, unknown for non-human.', default='unknown')
parser.add_argument('-male_sex_chr', metavar='str', type=str, help='Optional. Specify name(s) of the male-specific chromosome(s). Default: chrY,Y', default='chrY,Y')
parser.add_argument('-female_sex_chr', metavar='str', type=str, help='Optional. Specify name(s) of the chromosome(s) that is diploid in female. Default: chrX,X', default='chrX,X')
parser.add_argument('-sample_name', metavar='str', type=str, help='Optional. Specify sample name which will be labeled in the VCF output. If not specified, BAM/CRAM filename will be output.')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./result_out', default='./result_out')
parser.add_argument('-monoallelic', help='Optional. Specify if you use monoalellic sample, such as inbread mouse strains and HAP1 cells.', action='store_true')
parser.add_argument('-unsorted', help='Optional. Specify if an input BAM/CRAM is not position sorted.', action='store_true')
parser.add_argument('-verylowdep', help='Optional. Specify if you use parameter settings for low depth (generally less than 10x).', action='store_true')
parser.add_argument('-lowdep', help='Optional. Specify if you use parameter settings for low depth (generally less than 15x).', action='store_true')
parser.add_argument('-setting', metavar='str', type=str, help='Optional. Specify full path to the parameter setting file. Default: /path/to/prog/lib/parameter_settings.txt')
parser.add_argument('-only_ins', help='Optional. Specify if you only analyze non-reference MEI insertions.', action='store_true')
parser.add_argument('-only_abs', help='Optional. Specify if you only analyze absence of reference MEI.', action='store_true')
parser.add_argument('-pybedtools_tmp', metavar='str', type=str, help='Optional. Specify directory for temporary bedtools files, e.g. /dev/shm')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-do_not_overwrite', help='Optional. Specify if you do NOT overwrite previous results.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Optional. Number of threads. 3 or more is recommended. Default: 2', default=2)
parser.add_argument('-v', '--version', action='version', version='MEGAnE %s %s' % (os.path.basename(__file__), version))
parser.add_argument('-only_geno', help='Optional. Specify if you only genotype using previously analyzed data.', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('-threshold', metavar='int', type=int, help='Optional. Specify user-defined threshold.', help=argparse.SUPPRESS)
parser.add_argument('-only_geno_precall', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('-skip_unmapped', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('-make_sex_auto', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('-no_pdf', help='Optional. Specify if you do not want to output pdf summary files.', action='store_true', help=argparse.SUPPRESS)
args=parser.parse_args()
args.version=version


# start
import init
init.init(args, version)


# logging
import log
if args.only_geno_precall is True:
    args.only_geno=True
if args.only_geno is False:
    args.logfilename='for_debug.log'
else:
    args.logfilename='for_debug_geno.log'
if os.path.exists(os.path.join(args.outdir, args.logfilename)) is True:
    os.remove(os.path.join(args.outdir, args.logfilename))
log.start_log(args)
log.logger.debug('Logging started.')


# initial check
import initial_check
log.logger.debug('This is %s version %s' % (__file__, version))
print()
log.logger.info('Initial check started.')
initial_check.check(args, sys.argv)


# set up
import setup,auto_setting
args.auto=auto_setting.init(args)
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

filenames.repdb           =os.path.join(args.outdir, 'repdb')
filenames.repout_bed      =os.path.join(args.outdir, 'repout.bed')
filenames.rep_unknown_fa  =os.path.join(args.outdir, 'rep_unknown.fa')
filenames.blast_tmp_res   =os.path.join(args.outdir, 'blastn_tmp.txt')
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

filenames.blast1_res      =os.path.join(args.outdir, 'blastn_out_overhang_to_rep.txt')
filenames.overhang_MEI    =os.path.join(args.outdir, 'overhang_to_MEI_list.txt')
filenames.additional_pA   =os.path.join(args.outdir, 'overhang_additional_pA.txt')
filenames.blast2_res      =os.path.join(args.outdir, 'blastn_out_unmapped_to_rep.txt')
filenames.unmapped_hit_fa =os.path.join(args.outdir, 'unmapped_ref_side.fa')
filenames.blast3_res      =os.path.join(args.outdir, 'blastn_out_unmapped_mei_to_ref.txt')
filenames.unmapped_MEI    =os.path.join(args.outdir, 'unmapped_to_MEI_list.txt')
filenames.hybrid          =os.path.join(args.outdir, 'hybrid_reads.txt')

filenames.breakpoint_pairs=os.path.join(args.outdir, 'breakpoint_pairs.txt')
filenames.breakpoint_info =os.path.join(args.outdir, 'breakpoint_pairs_info.txt')
filenames.breakpoint_clean=os.path.join(args.outdir, 'breakpoint_pairs_clean.txt')

filenames.mapped_fa_select=os.path.join(args.outdir, 'mapped_selected.fa')
filenames.blast4_res      =os.path.join(args.outdir, 'blastn_out_mapped_to_ref.txt')
filenames.bp_pair_single  =os.path.join(args.outdir, 'breakpoint_pairs_in_TE_singletons.txt')
filenames.bp_info_single  =os.path.join(args.outdir, 'breakpoint_pairs_in_TE_singletons_info.txt')

filenames.bp_merged       =os.path.join(args.outdir, 'breakpoint_pairs_pooled_merged.txt')
filenames.hybrid_master   =os.path.join(args.outdir, 'hybrid_reads_master.txt')  # one of the final outputs
filenames.bp_merged_all   =os.path.join(args.outdir, 'breakpoint_pairs_pooled_all.txt')
filenames.bp_merged_filt_g=os.path.join(args.outdir, 'breakpoint_pairs_pooled_filtered_gaussian.txt')
filenames.bp_merged_filt_p=os.path.join(args.outdir, 'breakpoint_pairs_pooled_filtered_percentile.txt')
filenames.bp_merged_filt_f=os.path.join(args.outdir, 'breakpoint_pairs_pooled_filtered_failed.txt')
filenames.bp_merged_filt_u=os.path.join(args.outdir, 'breakpoint_pairs_pooled_filtered_user.txt')
filenames.bp_merged_groupg=os.path.join(args.outdir, 'breakpoint_pairs_pooled_grouped_gaussian.txt')
filenames.bp_merged_groupp=os.path.join(args.outdir, 'breakpoint_pairs_pooled_grouped_percentile.txt')
filenames.bp_merged_groupf=os.path.join(args.outdir, 'breakpoint_pairs_pooled_grouped_failed.txt')
filenames.bp_merged_groupu=os.path.join(args.outdir, 'breakpoint_pairs_pooled_grouped_user.txt')
filenames.gaussian_plot   =os.path.join(args.outdir, 'plot_gaussian_fitting.pdf')

filenames.bp_final_g      =os.path.join(args.outdir, 'MEI_final_gaussian.bed')
filenames.bp_final_p      =os.path.join(args.outdir, 'MEI_final_percentile.bed')
filenames.bp_final_f      =os.path.join(args.outdir, 'MEI_final_failed.bed')
filenames.bp_final_u      =os.path.join(args.outdir, 'MEI_final_user.bed')
filenames.bp_final_d      =os.path.join(args.outdir, 'MEI_final_dummy.bed')
filenames.transd_master   =os.path.join(args.outdir, '3transduction_check_master.txt')

filenames.abs_res         =os.path.join(args.outdir, 'absent_MEs.bed')
filenames.transd_res      =os.path.join(args.outdir, 'absent_MEs_transduction.bed')
filenames.abs_dummmy      =os.path.join(args.outdir, 'absent_MEs_dummy.bed')


if args.only_geno is False:
    # 0. preprocess repbase file
    import reshape_rep, blastn
    print()
    log.logger.info('Preprocess started.')
    reshape_rep.reshape_repout_to_bed(args, filenames)
    reshape_rep.reshape(args, params, filenames)
    if do_ins is True:
        blastn.makeblastdb(filenames.reshaped_rep, filenames.repdb)
        reshape_rep.slide_rep_file(args, params, filenames)
        blastn.blastn(args, params, filenames.rep_slide_file, filenames.repdb, filenames.blast0_res)  # params, q_path, db_path, outfpath
        reshape_rep.parse_slide_rep_blastn_res(args, filenames)
        # del files
        os.remove(filenames.blast0_res)
        os.remove(filenames.rep_slide_file)
    
    # 1. process unmapped overhangs
    import parse_blastn_result, find_additional_pA, extract_discordant, extract_discordant_c
    log.logger.info('Discordant read search started.')
    if args.p >= 2:
        from multiprocessing import Pool
        def extract_discordant_exe(n):
            extract_discordant_c.main(args, params, filenames, n)  ### c
        with Pool(args.p) as p:
            p.map(extract_discordant_exe, range(args.p))
        if do_abs is True:
            extract_discordant.concat_for_abs(args, filenames)
    else:
        extract_discordant_c.main(args, params, filenames, None)  ### c
    if do_ins is True:
        log.logger.info('Clipped read processing started.')
        if args.p >= 2:
            def blast_exe(n):
                infpath =filenames.overhang_fa + str(n) + '.txt'
                outfpath=filenames.blast1_res  + str(n) + '.txt'
                blastn.blastn_single_thread(args, params, infpath, filenames.repdb, outfpath)
            with Pool(args.p) as p:
                p.map(blast_exe, range(args.p))
            extract_discordant.concat_for_ins(args, filenames)
        else:
            blastn.blastn(args, params, filenames.overhang_fa, filenames.repdb, filenames.blast1_res)
        parse_blastn_result.parse(params, filenames.blast1_res, filenames.overhang_MEI)
        find_additional_pA.find(params, filenames.blast1_res, filenames.overhang_fa, filenames.additional_pA)
        # del files
        utils.gzip_or_del(args, params, filenames.blast1_res)

        # 2. process unmapped reads
        log.logger.info('Unmapped read processing started.')
        if args.skip_unmapped is True:
            log.logger.debug('Skips unmapped read processing.')
            with open(filenames.unmapped_MEI, 'w') as outfile:
                pass  # dummy file
        else:
            blastn.blastn(args, params, filenames.unmapped_fa, filenames.repdb, filenames.blast2_res)
            parse_blastn_result.unmapped_to_fa(params, filenames.unmapped_fa, filenames.blast2_res, filenames.unmapped_hit_fa)
            utils.gzip_or_del(args, params, filenames.blast2_res)
            blastn.blastn(args, params, filenames.unmapped_hit_fa, args.fadb, filenames.blast3_res)
            parse_blastn_result.find_chimeric_unmapped(args, params, filenames.blast3_res, filenames.unmapped_MEI)
            # del files
            utils.gzip_or_del(args, params, filenames.blast3_res)
            utils.gzip_or_del(args, params, filenames.unmapped_hit_fa)
        utils.gzip_or_del(args, params, filenames.unmapped_fa)
        
        # 3. process hybrid reads
        import process_distant_read
        log.logger.info('Hybrid read processing started.')
        process_distant_read.process_reads(args, params, filenames)
        # del files
        utils.gzip_file(params, filenames.distant_txt)

        # 4. merge all results, identify MEI outside of similar MEs
        import pair_breakpoints
        log.logger.info('Integration junction search (outside of TEs) started.')
        pair_breakpoints.pairing(args, params, filenames)
        pair_breakpoints.add_TE_subclass(args, filenames, filenames.breakpoint_pairs, filenames.breakpoint_info)
        pair_breakpoints.remove_cand_inside_TE(args, params, filenames)

        # 5. identify MEI nested in similar MEs
        import process_mapped_seq
        log.logger.info('Integration junction search (nested in TEs) started.')
        process_mapped_seq.retrieve_mapped_seq(params, filenames)
        utils.gzip_or_del(args, params, filenames.mapped_fa)  # del file
        process_mapped_seq.blastn_for_mapped(args, params, filenames.mapped_fa_select, args.fadb, filenames.blast4_res)
        process_mapped_seq.pairing(params, filenames)
        pair_breakpoints.add_TE_subclass(args, filenames, filenames.bp_pair_single, filenames.bp_info_single)
        # del files
        utils.gzip_or_del(args, params, filenames.mapped_fa_select)
        utils.gzip_or_del(args, params, filenames.blast4_res)

        # 6. pool results and filter candidates
        import pool_results
        log.logger.info('Filtering started.')
        pool_results.merge_breakpoints(filenames)
        pool_results.add_hybrid(args, params, filenames)
        import filter_candidates
        filter_candidates.filter(args, params, filenames)
        args.gaussian_executed=filter_candidates.gaussian_executed
        filter_candidates.grouping(args, filenames)
        import after_processing
        after_processing.grouped_mei_to_bed(args, params, filenames)
        after_processing.retrieve_3transd_reads(args, params, filenames)
        # del files
        utils.gzip_or_del(args, params, filenames.overhang_fa)
        utils.gzip_or_del(args, params, filenames.similar_rep_list)
        utils.gzip_or_del(args, params, filenames.breakpoint_clean)
        utils.gzip_or_del(args, params, filenames.breakpoint_info)
        utils.gzip_or_del(args, params, filenames.bp_info_single)
        utils.gzip_or_del(args, params, filenames.bp_pair_single)
        utils.gzip_or_del(args, params, filenames.bp_merged)
        if os.path.exists(filenames.bp_merged_filt_g) is True:
            utils.gzip_or_del(args, params, filenames.bp_merged_filt_g)
        if os.path.exists(filenames.bp_merged_filt_p) is True:
            utils.gzip_or_del(args, params, filenames.bp_merged_filt_p)
        if os.path.exists(filenames.bp_merged_filt_f) is True:
            utils.gzip_or_del(args, params, filenames.bp_merged_filt_f)
        if os.path.exists(filenames.bp_merged_groupg) is True:
            utils.gzip_file(params, filenames.bp_merged_groupg)
        if os.path.exists(filenames.bp_merged_groupp) is True:
            utils.gzip_file(params, filenames.bp_merged_groupp)
        if os.path.exists(filenames.bp_merged_groupf) is True:
            utils.gzip_file(params, filenames.bp_merged_groupf)
        utils.gzip_file(params, filenames.breakpoint_pairs)
        utils.gzip_file(params, filenames.bp_merged_all)
        utils.gzip_file(params, filenames.overhang_MEI)
        utils.gzip_file(params, filenames.unmapped_MEI)
        utils.gzip_file(params, filenames.overhang_pA)
        utils.gzip_file(params, filenames.additional_pA)
        utils.gzip_file(params, filenames.hybrid)
        utils.gzip_file(params, filenames.hybrid_master)
        if args.keep is not True:
            os.remove(filenames.distant_txt +'.gz')
        # output comments
        log.logger.info('%s ME insertion candidates found.' % filter_candidates.ins_ns)
        log.logger.info('ME insertion search finished!')


    # 7. search for absent MEs
    if do_abs is True:
        import find_absent
        print()
        log.logger.info('Absent ME search started.')
        find_absent.find_abs(args, params, filenames)
        # gzip files
        utils.gzip_file(params, filenames.abs_txt)
        # output comments
        log.logger.info('%d absent ME candidates found.' % find_absent.abs_n)
        log.logger.info('Absent ME search finished!')


    # remove unnecessary files
    os.remove(filenames.repout_bed)
    os.remove(filenames.reshaped_rep)
    for ext in ['nhr', 'nin', 'nog', 'nsd', 'nsi', 'nsq']:
        os.remove('%s.%s' % (filenames.repdb, ext))
    log.logger.info('pME search finished!')
    print()

else:
    log.logger.info('pME search skipped.')


# 8. start genotyping
setup.setup_geno_only_load_params(args, init.base)
params=setup.params
params.chrX=setup.chrX
params.chrY=setup.chrY
params.female=setup.female
params.male=setup.male


if do_abs is True:
    if os.path.exists(filenames.abs_dummmy) is True:
        args.abs_bed=filenames.abs_dummmy
    else:
        args.abs_bed=filenames.abs_res
    if os.path.exists(filenames.transd_res) is True:
        args.abs_3t_bed=filenames.transd_res
    else:
        args.abs_3t_bed=None

filenames.limited_b       =os.path.join(args.outdir, 'only_necessary.bam')
filenames.limited_c       =os.path.join(args.outdir, 'only_necessary.cram')
filenames.limited_tb      =os.path.join(args.outdir, 'only_necessary_tmp.bam')
filenames.limited_tc      =os.path.join(args.outdir, 'only_necessary_tmp.cram')
filenames.depth           =os.path.join(args.outdir, 'read_depth.txt')
filenames.tmp_bam         =os.path.join(args.outdir, 'tmp.bam')

filenames.merged_pdf_abs  =os.path.join(args.outdir, 'plot_out_genotyping_absents.pdf')
if do_abs is True:
    base=os.path.splitext(os.path.basename(args.abs_bed))[0]
    filenames.abs_out_bed     =os.path.join(args.outdir, '%s_genotyped.bed' % base)
    filenames.abs_out_vcf     =os.path.join(args.outdir, '%s_genotyped.vcf' % base)


# monoallelic, output VCF
if args.monoallelic is True:
    import output_genotyped_vcf_mono
    for f in [filenames.bp_final_g, filenames.bp_final_p, filenames.bp_final_f, filenames.bp_final_u, filenames.bp_final_d]:
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
        output_genotyped_vcf_mono.output_abs_bed_vcf(args, params, filenames)
        log.logger.info('Did output VCF, absent MEs.')
    # cleanup
    if not args.pybedtools_tmp == args.outdir:
        shutil.rmtree(args.pybedtools_tmp)
    # all finish!
    utils.output_finish_comment(args, do_ins, do_abs, filenames)
    exit(0)


#  9. limit BAM/CRAM
import output_genotyped_vcf
import allele_count_ins
log.logger.info('Limit BAM/CRAM started.')
data=utils.empclass()
allele_count_ins.limit(args, params, filenames, data)


# 10. genotype insertions
def genotype_ins(args, params, filenames, data):
    # filenames
    base=os.path.splitext(os.path.basename(args.ins_bed))[0]
    filenames.ins_out_bed   =os.path.join(args.outdir, '%s_genotyped.bed' % base)
    filenames.ins_out_vcf   =os.path.join(args.outdir, '%s_genotyped.vcf' % base)
    
    filenames.out_spanning  =os.path.join(args.outdir, '%s_spanning_read_summary.txt.gz' % base)
    filenames.disc_read_pdf =os.path.join(args.outdir, '%s_discordant_read_num.pdf' % base)
    filenames.debug_pdf1    =os.path.join(args.outdir, 'plot_out_%s_genotype_ins_for_debug.pdf' % base)
    filenames.merged_pdf    =os.path.join(args.outdir, 'plot_out_%s_genotyping_insertions.pdf' % base)
    
    log.logger.info('Evidence search started, insertion: %s' % args.ins_bed)
    allele_count_ins.evaluate_tsd_depth(args, params, filenames)
    data.cn_est_tsd_depth=allele_count_ins.cn_est_tsd_depth
    data.tsd_thresholds=allele_count_ins.tsd_thresholds
    data.del_thresholds=allele_count_ins.del_thresholds

    allele_count_ins.evaluate_spanning_read(args, params, filenames, data)
    data.cn_est_spanning=allele_count_ins.cn_est_spanning
    data.spanning_thresholds=allele_count_ins.spanning_thresholds
    
    if args.only_geno_precall is False:
        allele_count_ins.evaluate_discordant(args, params, filenames)
        data.cn_est_disc=allele_count_ins.cn_est_disc
        data.disc_thresholds=allele_count_ins.disc_thresholds  # could be False
    else:
        log.logger.debug('-only_geno_precall specified, skip discordant reads evaluation.')
    
    # merge evidences; insertion
    import merge_allele_evidence_ins
    log.logger.info('Evidence merge started, insertion: %s' % args.ins_bed)
#    merge_allele_evidence_ins.plot_orig(args, params, filenames, data)
    if args.only_geno_precall is False:
        merge_allele_evidence_ins.merge(args, params, filenames, data)
        data.merged_res=merge_allele_evidence_ins.merged_res
        merge_allele_evidence_ins.plot_merged(args, params, filenames, data)
    else:
        merge_allele_evidence_ins.merge_wo_discordant(args, params, filenames, data)
        data.merged_res=merge_allele_evidence_ins.merged_res
        merge_allele_evidence_ins.plot_merged_wo_disc(args, params, filenames, data)
    data.mei_filter=merge_allele_evidence_ins.mei_filter
    #merge_allele_evidence_ins.plot_gt(args, params, filenames, data)   # always commentout unless debug

    # output; insertion
    output_genotyped_vcf.output_ins_bed_vcf(args, params, filenames, data)
    log.logger.info('Did output VCF, insertion: %s' % args.ins_bed)

    # delete unnecessary files
    if os.path.exists(filenames.tmp_bam) is True:
        os.remove(filenames.tmp_bam)

for f in [filenames.bp_final_g, filenames.bp_final_p, filenames.bp_final_f, filenames.bp_final_u, filenames.bp_final_d]:
    if os.path.exists(f) is True:
        args.ins_bed=f
        genotype_ins(args, params, filenames, data)


# 11. genotype absent MEs
if do_abs is True:
    # filenames
    filenames.abs_spanning  =os.path.join(args.outdir, 'abs_spanning_read_summary.txt.gz')
    
    import allele_count_abs
    log.logger.info('Evidence search started, absent MEs.')
    if args.only_geno_precall is True:
        tmp_disc_ids=data.disc_ids
        data=utils.empclass()
        data.disc_ids=tmp_disc_ids
    else:
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
    if args.only_geno_precall is False:
        merge_allele_evidence_abs.merge(args, params, filenames, data)
    else:
        merge_allele_evidence_abs.merge_wo_discordant(args, params, filenames, data)
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

# genotyping finish
log.logger.info('Genotyping finished!')

# cleanup
if not args.pybedtools_tmp == args.outdir:
    shutil.rmtree(args.pybedtools_tmp)

# all finish!
utils.output_finish_comment(args, do_ins, do_abs, filenames)
