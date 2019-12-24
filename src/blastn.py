#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

from Bio.Blast.Applications import NcbiblastnCommandline

def blastn(params, q_path, db_path, outfpath):
    NcbiblastnCommandline(db=db_path, query=q_path, evalue=params.blastn_evalue, perc_identity=params.blastn_ident, word_size=params.blastn_word_size, num_threads=params.p, outfmt=6, out=outfpath)()
    pass
