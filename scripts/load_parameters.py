#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

class load:
    def __init__(self, f):
        # default
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
        self.eval_threshold_for_gaussian_fitting=float('1e-25')
        self.fit_gaussian_init_a_coeff=0.5
        self.fit_gaussian_init_mu_coeff=1
        self.fit_gaussian_init_sigma_coeff=0.33
        self.fit_gaussian_CI_alpha=0.99
        self.first_filter_eval_threshold=float('1e-15')
        self.first_filter_total_hybrid_read_num=1
        self.second_filter_hybrid_read_num=1
        self.second_filter_eval_threshold_for_few_hybrid=float('1e-25')
        self.L1_filter_min_TSD_len=5
        self.L1_filter_A_or_T_perc=50
        self.L1_filter_A_plus_T_perc=90
        self.L1_filter_eval_threshold=float('1e-25')
        # read parameter setting file
        with open(f) as infile:
            for line in infile:
                ls=line.strip().split('=')
                if ls[0] == 'discordant_reads_clip_len':
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
                elif ls[0] == 'eval_threshold_for_gaussian_fitting':
                    self.eval_threshold_for_gaussian_fitting=float(ls[1])
                elif ls[0] == 'fit_gaussian_init_a_coeff':
                    self.fit_gaussian_init_a_coeff=float(ls[1])
                elif ls[0] == 'fit_gaussian_init_mu_coeff':
                    self.fit_gaussian_init_mu_coeff=float(ls[1])
                elif ls[0] == 'fit_gaussian_init_sigma_coeff':
                    self.fit_gaussian_init_sigma_coeff=float(ls[1])
                elif ls[0] == 'fit_gaussian_CI_alpha':
                    self.fit_gaussian_CI_alpha=float(ls[1])
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
