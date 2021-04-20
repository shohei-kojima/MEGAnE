#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import log,traceback

class load:
    def __init__(self, args, f):
        log.logger.debug('started')
        try:
            # default; polymorphic ME search
            self.read_count_for_readlen_estimation=200
            self.chr1_start_depth_est=112500000
            self.chr1_end_depth_est=112600000
            self.chrX_start_depth_est=20000000
            self.chrX_end_depth_est=20100000
            self.chrY_start_depth_est=6900000
            self.chrY_end_depth_est=7000000
            self.sex_est_XY_ratio_threshold=0.3
            self.sex_est_XY_ratio_threshold_for_nochrY=0.75
            self.discordant_reads_clip_len=20
            self.read_pair_gap_len=2000
            self.max_TSD_len=50
            self.polyA_overhang_threshold=0.7
            self.mapped_region_low_complex_threshold=0.7
            self.abs_min_dist=50
            self.abs_max_dist=20000
            self.blastn_evalue=float('1e-05')
            self.blastn_ident=80
            self.blastn_word_size=11
            self.overhang_evalue_threshold=float('1e-05')
            self.gzip_compresslevel=1
            self.pA_scan_bin=12
            self.max_non_pA_count=2
            self.scan_loop_from_edge=5
            self.max_ref_genome_hits_for_unmapped=20
            self.repbase_seq_slide_bin=5
            self.min_read_num_per_breakpoint_edge=1
            self.max_breakpoint_gap=50
            self.max_TSD_len=50
            self.ref_TE_slop_len=0
            self.retrieve_mapped_seq_threshold=3
            self.blastn_evalue_for_mapped=float('1e-05')
            self.blastn_ident_for_mapped=95
            self.blastn_word_size_for_mapped=30
            self.mapped_abs_single_ident_threshold=98
            self.hybrid_read_range_from_breakpint=500
            self.hybrid_read_coeff_for_gaussian_fitting=0.1
            self.chimeric_read_coeff_for_gaussian_fitting=0.01
            self.eval_threshold_for_gaussian_fitting=float('1e-25')
            if 100 <= args.readlen <= 149:
                self.eval_threshold_for_gaussian_fitting=float('1e-15')
            elif args.readlen <= 99:
                self.eval_threshold_for_gaussian_fitting=float('1e-10')
            self.fit_gaussian_init_a_coeff=0.5
            self.fit_gaussian_init_mu_coeff=1
            self.fit_gaussian_init_sigma_coeff=0.33
            self.fit_gaussian_CI_alpha=0.99
            self.actual_cutoff_rank=0.001
            self.first_filter_eval_threshold=float('1e-15')
            self.first_filter_total_hybrid_read_num=1
            self.second_filter_hybrid_read_num=1
            self.second_filter_eval_threshold_for_few_hybrid=float('1e-25')
            self.L1_filter_min_TSD_len=5
            self.L1_filter_A_or_T_perc=50
            self.L1_filter_A_plus_T_perc=90
            self.L1_filter_eval_threshold=float('1e-25')
            self.abs_min_chimeric_num_coeff=0.03
            self.breakpoint_annotation_gap=25
            self.abs_len_to_te_ratio=0.9
            self.len_te_for_abs_ratio=0.9
            self.non_ME_len_ratio=0.5
            self.transduction_pA_len=12
            self.transduction_pA_ratio=0.75
            self.length_for_3transduction_search=1000
            # read parameter setting file
            with open(f) as infile:
                for line in infile:
                    ls=line.strip().split(' ')[0].split('=')
                    if ls[0] == 'read_count_for_readlen_estimation':
                        self.read_count_for_readlen_estimation=int(ls[1])
                    elif ls[0] == 'chr1_start_depth_est':
                        self.chr1_start_depth_est=int(ls[1])
                    elif ls[0] == 'chr1_end_depth_est':
                        self.chr1_end_depth_est=int(ls[1])
                    elif ls[0] == 'chrX_start_depth_est':
                        self.chrX_start_depth_est=int(ls[1])
                    elif ls[0] == 'chrX_end_depth_est':
                        self.chrX_end_depth_est=int(ls[1])
                    elif ls[0] == 'chrY_start_depth_est':
                        self.chrY_start_depth_est=int(ls[1])
                    elif ls[0] == 'chrY_end_depth_est':
                        self.chrY_end_depth_est=int(ls[1])
                    elif ls[0] == 'sex_est_XY_ratio_threshold':
                        self.sex_est_XY_ratio_threshold=float(ls[1])
                    elif ls[0] == 'discordant_reads_clip_len':
                        self.discordant_reads_clip_len=int(ls[1])
                    elif ls[0] == 'read_pair_gap_len':
                        self.read_pair_gap_len=int(ls[1])
                    elif ls[0] == 'max_TSD_len':
                        self.max_TSD_len=int(ls[1])
                    elif ls[0] == 'polyA_overhang_threshold':
                        self.polyA_overhang_threshold=float(ls[1])
                    elif ls[0] == 'mapped_region_low_complex_threshold':
                        self.mapped_region_low_complex_threshold=float(ls[1])
                    elif ls[0] == 'abs_min_dist':
                        self.abs_min_dist=int(ls[1])
                    elif ls[0] == 'abs_max_dist':
                        self.abs_max_dist=int(ls[1])
                    elif ls[0] == 'blastn_evalue':
                        self.blastn_evalue=float(ls[1])
                    elif ls[0] == 'blastn_ident':
                        self.blastn_ident=int(ls[1])
                    elif ls[0] == 'blastn_word_size':
                        self.blastn_word_size=int(ls[1])
                    elif ls[0] == 'overhang_evalue_threshold':
                        self.overhang_evalue_threshold=float(ls[1])
                    elif ls[0] == 'gzip_compresslevel':
                        self.gzip_compresslevel=int(ls[1])
                    elif ls[0] == 'pA_scan_bin':
                        self.pA_scan_bin=int(ls[1])
                    elif ls[0] == 'max_non_pA_count':
                        self.max_non_pA_count=int(ls[1])
                    elif ls[0] == 'scan_loop_from_edge':
                        self.scan_loop_from_edge=int(ls[1])
                    elif ls[0] == 'max_ref_genome_hits_for_unmapped':
                        self.max_ref_genome_hits_for_unmapped=int(ls[1])
                    elif ls[0] == 'repbase_seq_slide_bin':
                        self.repbase_seq_slide_bin=int(ls[1])
                    elif ls[0] == 'min_read_num_per_breakpoint_edge':
                        self.min_read_num_per_breakpoint_edge=int(ls[1])
                    elif ls[0] == 'max_breakpoint_gap':
                        self.max_breakpoint_gap=int(ls[1])
                    elif ls[0] == 'max_TSD_len':
                        self.max_TSD_len=int(ls[1])
                    elif ls[0] == 'ref_TE_slop_len':
                        if not int(ls[1]) < 0:
                            self.ref_TE_slop_len=int(ls[1])
                        else:
                            print('Warning: ref_TE_slop_len must be positive value. Changed to 0 and continue this analysis.')
                    elif ls[0] == 'retrieve_mapped_seq_threshold':
                        self.retrieve_mapped_seq_threshold=int(ls[1])
                    elif ls[0] == 'blastn_evalue_for_mapped':
                        self.blastn_evalue_for_mapped=float(ls[1])
                    elif ls[0] == 'blastn_ident_for_mapped':
                        self.blastn_ident_for_mapped=int(ls[1])
                    elif ls[0] == 'blastn_word_size_for_mapped':
                        self.blastn_word_size_for_mapped=int(ls[1])
                    elif ls[0] == 'mapped_abs_single_ident_threshold':
                        self.mapped_abs_single_ident_threshold=int(ls[1])
                    elif ls[0] == 'hybrid_read_range_from_breakpint':
                        if not int(ls[1]) < self.max_TSD_len:
                            self.hybrid_read_range_from_breakpint=int(ls[1])
                        else:
                            print('Warning: hybrid_read_range_from_breakpint must be larger than max_TSD_len. Changed to max_TSD_len and continue this analysis.')
                            self.hybrid_read_range_from_breakpint=self.max_TSD_len
                    elif ls[0] == 'hybrid_read_coeff_for_gaussian_fitting':
                        self.hybrid_read_coeff_for_gaussian_fitting=float(ls[1])
                    elif ls[0] == 'chimeric_read_coeff_for_gaussian_fitting':
                        self.chimeric_read_coeff_for_gaussian_fitting=float(ls[1])
#                    elif ls[0] == 'eval_threshold_for_gaussian_fitting':
#                        self.eval_threshold_for_gaussian_fitting=float(ls[1])
                    elif ls[0] == 'fit_gaussian_init_a_coeff':
                        self.fit_gaussian_init_a_coeff=float(ls[1])
                    elif ls[0] == 'fit_gaussian_init_mu_coeff':
                        self.fit_gaussian_init_mu_coeff=float(ls[1])
                    elif ls[0] == 'fit_gaussian_init_sigma_coeff':
                        self.fit_gaussian_init_sigma_coeff=float(ls[1])
                    elif ls[0] == 'fit_gaussian_CI_alpha':
                        self.fit_gaussian_CI_alpha=float(ls[1])
                    elif ls[0] == 'actual_cutoff_rank':
                        self.actual_cutoff_rank=float(ls[1])
                    elif ls[0] == 'first_filter_eval_threshold':
                        self.first_filter_eval_threshold=float(ls[1])
                    elif ls[0] == 'first_filter_total_hybrid_read_num':
                        self.first_filter_total_hybrid_read_num=int(ls[1])
                    elif ls[0] == 'second_filter_hybrid_read_num':
                        self.second_filter_hybrid_read_num=int(ls[1])
                    elif ls[0] == 'second_filter_eval_threshold_for_few_hybrid':
                        self.second_filter_eval_threshold_for_few_hybrid=float(ls[1])
                    elif ls[0] == 'L1_filter_min_TSD_len':
                        self.L1_filter_min_TSD_len=int(ls[1])
                    elif ls[0] == 'L1_filter_A_or_T_perc':
                        self.L1_filter_A_or_T_perc=int(ls[1])
                    elif ls[0] == 'L1_filter_A_plus_T_perc':
                        self.L1_filter_A_plus_T_perc=int(ls[1])
                    elif ls[0] == 'L1_filter_eval_threshold':
                        self.L1_filter_eval_threshold=float(ls[1])
                    elif ls[0] == 'abs_min_chimeric_num_coeff':
                        self.abs_min_chimeric_num_coeff=float(ls[1])
                    elif ls[0] == 'breakpoint_annotation_gap':
                        self.breakpoint_annotation_gap=int(ls[1])
                    elif ls[0] == 'transduction_pA_len':
                        self.transduction_pA_len=int(ls[1])
                    elif ls[0] == 'abs_len_to_te_ratio':
                        self.abs_len_to_te_ratio=float(ls[1])
                    elif ls[0] == 'len_te_for_abs_ratio':
                        self.len_te_for_abs_ratio=float(ls[1])
                    elif ls[0] == 'non_ME_len_ratio':
                        self.non_ME_len_ratio=float(ls[1])
                    elif ls[0] == 'transduction_pA_ratio':
                        self.transduction_pA_ratio=float(ls[1])
                    elif ls[0] == 'length_for_3transduction_search':
                        self.length_for_3transduction_search=int(ls[1])
            params_for_debug=[]
            for k,v in self.__dict__.items():
                params_for_debug.append('%s=%s' % (k, str(v)))
            log.logger.debug('parameters:\n'+ '\n'.join(params_for_debug))
        except:
            log.logger.error('\n'+ traceback.format_exc())
            exit(1)


class load_geno:
    def __init__(self, args, f):
        log.logger.debug('started')
        try:
            # default; genotyping
            self.ins_disc_ids_threshold_coeff=0.05
            self.abs_disc_ids_threshold=1
            self.ins_slop_len=300
            self.abs_slop_len=300
            self.ins_slop_len_for_disc_detection=10
            self.ins_discard_flank_len=5
            self.abs_slop_len_for_disc_detection=10
            self.min_tsd_len_to_remove_1nt=4
            self.tsd_only_highest_len_precall=10
            self.tsd_flank_len=4
            self.tsd_flank_len_for_precall=5
            self.ins_remove_multi_nt=1
            self.ins_remove_adj=1
            self.del_mono_high_threshold_coeff=0.666
            self.del_bi_high_threshold_coeff=1.333
            self.tsd_outlier_low_coeff=0.25
            self.tsd_outlier_low_for_precall=1.1
            self.tsd_threshold_correction=-0.07
            self.tsd_outlier=3
            self.del_outlier=2
            self.del_outlier_low_for_precall=0.8
            self.min_len_overhang_for_spanning=31
            self.min_overhang_match_for_spanning=31
            self.max_overhang_mismatch_ratio_for_spanning=0.1
            self.spanning_high_threshold_coeff=0.666
            self.spanning_zero_threshold_coeff=0.333
            self.spanning_outlier_coeff=1.5
            self.spanning_outlier_coeff_for_precall=0.66
            self.fit_gaussian_init_a_coeff=0.5
            self.fit_gaussian_init_mu_coeff=1
            self.fit_gaussian_init_sigma_coeff=0.33
            self.discordant_outlier_coeff=3.5
            self.spanning_threshold_coeff_for_merge=0.01
            self.custom_tsd_mono_high_conf_threshold=False
            self.custom_tsd_threshold=False
            self.custom_tsd_bi_high_conf_threshold=False
            self.custom_tsd_outlier=False
            self.custom_tsd_outlier_low_threshold=False
            self.custom_spanning_zero_threshold=False
            self.custom_spanning_high_threshold=False
            self.custom_spanning_outlier=False
            self.custom_disc_mono_high_conf_threshold=-1
            self.custom_disc_threshold=-1
            self.custom_disc_di_high_conf_threshold=-1
            self.custom_disc_outlier_threshold=False
            
            self.abs_slop_len=300
            self.min_len_overhang_for_spanning_abs=31
            self.min_overhang_match_for_spanning_abs=31
            self.max_overhang_mismatch_ratio_for_spanning_abs=0.1
            self.spanning_high_threshold_coeff_abs=0.666
            self.spanning_zero_threshold_coeff_abs=0.333
            self.spanning_outlier_coeff_abs=1.5
            self.spanning_outlier_coeff_for_precall_abs=0.5
            self.abs_flank_len=5
            self.abs_flank_len_for_precall=15
            self.abs_remove_multi_nt=2
            self.abs_remove_adj=2
            self.mono_peak_notfound=0.7
            self.abs_depth_outlier=0.8
            self.abs_depth_outlier_precall=0.75
            self.spanning_threshold_coeff_for_merge_abs=0.01
            self.custom_spanning_zero_threshold_abs=False
            self.custom_spanning_high_threshold_abs=False
            self.custom_spanning_outlier_abs=False
            self.custom_abs_mono_high_conf_threshold=False
            self.custom_abs_depth_outlier=False
            
            self.bam_sort_maxmem='2G'
            
            params_for_debug=[]
            for k,v in self.__dict__.items():
                params_for_debug.append('%s=%s' % (k, str(v)))
            log.logger.debug('parameters:\n'+ '\n'.join(params_for_debug))
        except:
            log.logger.error('\n'+ traceback.format_exc())
            exit(1)


class load_merge_vcf:
    def __init__(self, args):
        log.logger.debug('started')
        try:
            # default; merge vcf
            self.blastn_evalue=float('1e-05')
            self.blastn_ident=80
            self.blastn_word_size=11
            self.pass_perc_threshold=0.25
            self.overhang_evalue_threshold=float('1e-05')
            self.slop_len_for_abs=10
            self.min_support_reads_ins=1
            self.processed_interval=10
            
            params_for_debug=[]
            for k,v in self.__dict__.items():
                params_for_debug.append('%s=%s' % (k, str(v)))
            log.logger.debug('parameters:\n'+ '\n'.join(params_for_debug))
        except:
            log.logger.error('\n'+ traceback.format_exc())
            exit(1)
