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
        self.blastn_evalue=float('1e-05')
        self.blastn_ident=80
        self.blastn_word_size=11
        self.overhang_evalue_threshold=1e-05
        self.gzip_compresslevel=1
        self.pA_scan_bin=12
        self.max_non_pA_count=2
        self.scan_loop_from_edge=5
        self.max_ref_genome_hits_for_unmapped=20
        self.repbase_seq_slide_bin=5
        self.min_read_num_per_breakpoint_edge=1
        self.max_breakpoint_gap=50
        self.max_TSD_len=50
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
