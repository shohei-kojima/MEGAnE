#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os
import log,traceback


def parse(params, f_blast_res, outfpath):  
    log.logger.debug('started')
    try:
        evalue_threshold=params.overhang_evalue_threshold

        def determine_best_hit(list):
            evals=[]
            for l in list:
                evals.append(float(l[10]))
            eval=min(evals)
            if eval <= evalue_threshold:
                hits=[]
                for l in list:
                    if float(l[10]) == eval:
                        hits.append(l[1])
                return l[0], ';'.join(hits), str(eval)
            else:
                return 'NA', 'NA', 'NA'

        # main
        with open(outfpath, 'w') as outifle:
            currentq='any'
            tmp=[]
            with open(f_blast_res) as infile:
                for line in infile:
                    ls=line.split()
                    if ls[0] == currentq:
                        tmp.append(ls)
                    elif len(tmp) >= 1:
                        id,hits,eval=determine_best_hit(tmp)
                        if not id == 'NA':
                            outifle.write(id +'\t'+ hits +'\t'+ eval +'\n')
                        tmp=[ls]
                        currentq=ls[0]
                    else:
                        tmp=[ls]
                        currentq=ls[0]
                    pass
            id,hits,eval=determine_best_hit(tmp)
            if not id == 'NA':
                outifle.write(id +'\t'+ hits +'\t'+ eval +'\n')
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def unmapped_to_fa(params, unmapped_fa, blast_res, outfpath):
    log.logger.debug('started')
    try:
        fai={}
        res={}
        with open(blast_res) as infile:
            for line in infile:
                ls=line.split()
                if not ls[0] in fai:
                    fai[ls[0]]=[]
                    res[ls[0]]=[]
                res[ls[0]].append(line)

        # load query length
        with open(unmapped_fa) as infile:
            for line in infile:
                header=line.strip().replace('>', '')
                seq=next(infile).strip()
                if header in fai:
                    fai[header].append(len(seq))
                    fai[header].append(seq)

        # output remaining seqs where should be reference seq
        with open(outfpath, 'w') as outfile:
            for id in res:
                tmp={}
                seqlen,seq=fai[id]
                for line in res[id]:
                    ls=line.split()
                    qstart,qend= int(ls[6]) - 1, int(ls[7])  # 1-base to 0-base
                    sstart,send= int(ls[8]) - 1, int(ls[9])  # 1-base to 0-base
                    strand='+'
                    if sstart >= send:
                        strand='-'
                    rc= seqlen - qend
                    if (qstart >= params.discordant_reads_clip_len) or (rc >= params.discordant_reads_clip_len):
                        clipseq=''
                        if qstart > rc:
                            clipseq += seq[:qstart]
                            dir='L'
                        elif rc > qstart:
                            clipseq += seq[qend:]
                            dir='R'
                        if len(clipseq) >= 1:
                            if not clipseq in tmp:
                                tmp[clipseq]=[]
                            tmp[clipseq].append(id +'/'+ dir +'/'+ ls[1] +'/'+ strand)
                if len(tmp) >= 1:
                    for cs in tmp:
                        hs='>' + ';'.join(tmp[cs])
                        if len(hs) >= 2:
                            outfile.write(hs +'\n'+ cs +'\n')
            outfile.flush()
            os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def find_chimeric_unmapped(args, params, blast_res, outfpath):
    log.logger.debug('started')
    try:
        # load blast results
        res={}
        with open(blast_res) as infile:
            for line in infile:
                ls=line.split()
                if not ls[0] in res:
                    res[ls[0]]=[]
                res[ls[0]].append(line)

        # select hits with less than params.max_ref_genome_hits_for_unmapped
        res_few={}
        for id in res:
            if len(res[id]) <= params.max_ref_genome_hits_for_unmapped:
                res_few[id]=[]
                evals=[]
                for line in res[id]:
                    evals.append(float(ls[10]))
                eval=min(evals)
                for line in res[id]:
                    if float(ls[10]) == eval:
                        res_few[id].append(line)
        del(res)

        # reshape results
        tmp_d={}
        for id in res_few:
            refs=[]
            meis=[]
            for i in id.split(';'):
                tmp,dir=i.rsplit('/', 1)
                tmp,te=tmp.rsplit('/', 1)
                readname,breakpoint=tmp.rsplit('/', 1)
                if breakpoint == 'L':
                    breakpoint='R'
                elif breakpoint == 'R':
                    breakpoint='L'
                meis.append([breakpoint +'/'+ readname, dir, te])
            for line in res_few[id]:
                ls=line.split()
                if ls[1] in args.main_chrs_set:
                    sstart,send= int(ls[8]) - 1, int(ls[9])  # 1-base to 0-base
                    dir='+'
                    if sstart > send:
                        dir='-'
                        sstart,send=send,sstart
                    refs.append([ls[1] +':'+ str(sstart) +'-'+ str(send), dir, ls[10]])
            for mid,mdir,te in meis:
                for rid,rdir,eval in refs:
                    dir='-'
                    if mdir == rdir:
                        dir='+'
                    finalid= rid +'/'+ mid +'/'+ dir +';'
                    if not finalid in tmp_d:
                        tmp_d[finalid]={}
                    if not eval in tmp_d[finalid]:
                        tmp_d[finalid][eval]=[]
                    tmp_d[finalid][eval].append(te)

        # group same eval hits
        with open(outfpath, 'w') as outfile:
            for id in tmp_d:
                evals=list(tmp_d[id].keys())
                eval=min(evals)
                te=';'.join(tmp_d[id][eval])
                outfile.write(id +'\t'+ te +'\t'+ eval +'\n')
            outfile.flush()
            os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

