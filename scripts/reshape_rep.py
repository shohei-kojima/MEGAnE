#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os
from utils import parse_fasta
import blastn
import log,traceback


def reshape(args, params, filenames):
    log.logger.debug('started')
    try:
        fa=parse_fasta(args.rep)
        if len(fa) == 0:
            log.logger.error('No sequence found in your input rep file (-rep flag).')
            exit(1)
        fa_keep={}
        fa_unknown={}
        known_name_to_clas={}
        rep_kept_n=0
        rep_removed_n=0
        for header in fa:
            hs=header.split('\t')
            if len(hs) == 3:
                if not hs[1] in args.rep_headers_to_be_removed:
                    fa_keep[header]=fa[header]
                    known_name_to_clas[hs[0].replace('>', '')]=hs[1]
                    rep_kept_n += 1
                else:
                    rep_removed_n += 1
            elif len(hs) == 1:
                fa_unknown[header]=fa[header]
                rep_kept_n += 1
            else:
                log.logger.error('Fasta header of your input rep file (-rep flag) does not look like the appropriate fasta file conveted from Dfam.')
                exit(1)
        if rep_kept_n == 0:
            log.logger.error('All sequences in your input rep file (-rep flag) match with non-MEs listed in a file specified with -repremove flag.')
            exit(1)
        else:
            log.logger.info('N=%d repeats found in %s. N=%d will be analyzed. N=%d will be excluded due to non-ME repeats.' % (len(fa), args.rep, rep_kept_n, rep_removed_n))
        with open(filenames.reshaped_rep, 'w') as outfile:
            for header in fa_keep:
                outfile.write('%s\n%s\n' % (header, fa_keep[header]))
            outfile.flush()
            os.fdatasync(outfile.fileno())
        dels=set()
        for header in fa_unknown:
            h=header.replace('>', '')
            if h in known_name_to_clas:
                dels.add(header)
        if len(dels) >= 1:
            for d in dels:
                del(fa_unknown[d])
        if len(fa_unknown) >= 1:
            log.logger.debug('%d MEs with ambiguous subclass found.' % len(fa_unknown))
            with open(filenames.rep_unknown_fa, 'w') as outfile:
                for header in fa_unknown:
                    outfile.write('%s\n%s\n' % (header, fa_unknown[header]))
                outfile.flush()
                os.fdatasync(outfile.fileno())
            blastn.makeblastdb(filenames.reshaped_rep, filenames.repdb)
            blastn.blastn_for_unknown_rep_ident(args, params, filenames.rep_unknown_fa, filenames.repdb, filenames.blast_tmp_res)  # determine TE class of unknown rep
            hits={}
            with open(filenames.blast_tmp_res) as infile:
                for line in infile:
                    ls=line.split()
                    if not ls[0] == ls[1]:
                        if not ls[0] in hits:
                            hits[ls[0]]=ls[1]
            for header in hits:
                hit=hits[header]
                if hit in known_name_to_clas:
                    fa_keep['>%s\t%s' % (header, known_name_to_clas[hit])]=fa_unknown['>%s' % header]
            with open(filenames.reshaped_rep, 'w') as outfile:
                for header in fa_keep:
                    outfile.write(header +'\n'+ fa_keep[header] +'\n')
                outfile.flush()
                os.fdatasync(outfile.fileno())
            os.remove(filenames.rep_unknown_fa)
            os.remove(filenames.blast_tmp_res)
        else:
            log.logger.debug('No MEs with ambiguous subclass found.')
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def slide_rep_file(args, params, filenames):
    log.logger.debug('started')
    try:
        def bin_slide(string, bin, interval):
            l=[]
            if len(string) >= bin:
                for i in range(0, len(string)-bin+1, interval):
                    l.append(string[i:i+bin])
            else:
                l.append(string)
            return l

        fa=parse_fasta(filenames.reshaped_rep)
        with open(filenames.rep_slide_file, 'w') as outfile:
            for h in fa:
                seqs=bin_slide(fa[h], args.readlen, params.repbase_seq_slide_bin)
                tmp=[]
                n=0
                for seq in seqs:
                    tmp.append(h.strip() +'/'+ str(n))
                    tmp.append(seq)
                    n += 1
                outfile.write('\n'.join(tmp) +'\n')
            outfile.flush()
            os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def parse_slide_rep_blastn_res(args, filenames):
    log.logger.debug('started')
    try:
        d={}
        with open(filenames.blast0_res) as infile:
            for line in infile:
                ls=line.split()
                id=ls[0].split('/')[0]
                if not id in d:
                    d[id]=set()
                d[id].add(ls[1])
        for id in d:
            d[id]=sorted(list(d[id]))
        with open(filenames.similar_rep_list, 'w') as outfile:
            for id in d:
                outfile.write('%s\t%s\n' % (id, ';'.join(d[id])))
            outfile.flush()
            os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def reshape_repout_to_bed(args, filenames):
    log.logger.debug('started')
    try:
        non_rep_line={'SW', 'bit', 'score'}
        with open(filenames.repout_bed, 'w') as outfile:
            with open(args.repout) as infile:
                for line in infile:
                    ls=line.split()
                    if len(ls) < 15:
                        continue
                    if ls[0] in non_rep_line:
                        continue
                    start= int(ls[5]) - 1  # 0-based
                    strand='+' if ls[8] == '+' else '-'
                    outfile.write('%s\t%d\t%s\t%s:%s\t%s\n' % (ls[4], start, ls[6], ls[9], ls[10], strand))
            outfile.flush()
            os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

