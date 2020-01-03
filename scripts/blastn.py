#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import time
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline


def blastn(args, params, q_path, db_path, outfpath):
    NcbiblastnCommandline(db=db_path, query=q_path, evalue=params.blastn_evalue, perc_identity=params.blastn_ident, word_size=params.blastn_word_size, num_threads=args.p, outfmt=6, out=outfpath)()
    time.sleep(20)


def blastn_single_thread(args, params, q_path, db_path, outfpath):
    NcbiblastnCommandline(db=db_path, query=q_path, evalue=params.blastn_evalue, perc_identity=params.blastn_ident, word_size=params.blastn_word_size, num_threads=1, outfmt=6, out=outfpath)()
    time.sleep(20)


def blastn_for_unknown_rep_ident(args, params, q_path, db_path, outfpath):
    NcbiblastnCommandline(db=db_path, query=q_path, evalue=params.blastn_evalue, perc_identity=params.blastn_ident, word_size=7, num_threads=args.p, outfmt=6, out=outfpath)()
    time.sleep(20)


def makeblastdb(fasta_file, dbpath):
    NcbimakeblastdbCommandline(input_file=fasta_file, dbtype='nucl', out=dbpath, parse_seqids=True)()
    time.sleep(20)
