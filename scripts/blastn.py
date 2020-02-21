#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import log,traceback


def blastn(args, params, q_path, db_path, outfpath):
    log.logger.debug('started')
    try:
        NcbiblastnCommandline(db=db_path, query=q_path, evalue=params.blastn_evalue, perc_identity=params.blastn_ident, word_size=params.blastn_word_size, num_threads=args.p, outfmt=6, out=outfpath)()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def blastn_single_thread(args, params, q_path, db_path, outfpath):
    log.logger.debug('started')
    try:
        NcbiblastnCommandline(db=db_path, query=q_path, evalue=params.blastn_evalue, perc_identity=params.blastn_ident, word_size=params.blastn_word_size, num_threads=1, outfmt=6, out=outfpath)()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def blastn_for_unknown_rep_ident(args, params, q_path, db_path, outfpath):
    log.logger.debug('started')
    try:
        NcbiblastnCommandline(db=db_path, query=q_path, evalue=params.blastn_evalue, perc_identity=params.blastn_ident, word_size=7, num_threads=args.p, outfmt=6, out=outfpath)()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def makeblastdb(fasta_file, dbpath):
    log.logger.debug('started')
    try:
        NcbimakeblastdbCommandline(input_file=fasta_file, dbtype='nucl', out=dbpath, parse_seqids=True)()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
