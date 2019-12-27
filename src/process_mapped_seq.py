#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


def retrieve_mapped_seq(filenames):
    # load read names
    R_readnames,L_readnames={},{}
    R_read_lens,L_read_lens={},{}
    ids=[]
    with open(filenames.breakpoint_pairs) as infile:
        for line in infile:
            ls=line.split()
            id= ls[0] +':'+ ls[1] +'-'+ ls[2] +'/'+ ls[7]
            ids.append(id)
            R_read_lens[id]=[]
            L_read_lens[id]=[]
            if not ls[8] == 'NA':
                for readname in ls[8].split(';'):
                    if not readname in R_readnames:
                        R_readnames[readname]=set()
                    R_readnames[readname].add(id)
            if not ls[9] == 'NA':
                for readname in ls[9].split(';'):
                    if not readname in L_readnames:
                        L_readnames[readname]=set()
                    L_readnames[readname].add(id)

    # load longest reads
    with open(filenames.mapped_fa) as infile:
        for line in infile:
            seq=next(infile).strip()
            ls=line.strip().replace('>', '').split(';')[:-1]
            for readname in ls:
                if readname in R_readnames:
                    for id in R_readnames[readname]:
                        if len(R_read_lens[id]) == 0:
                            R_read_lens[id]=[readname, seq]
                        elif len(seq) > len(R_read_lens[id][1]):
                            R_read_lens[id]=[readname, seq]
                if readname in L_readnames:
                    for id in L_readnames[readname]:
                        if len(L_read_lens[id]) == 0:
                            L_read_lens[id]=[readname, seq]
                        elif len(seq) > len(L_read_lens[id][1]):
                            L_read_lens[id]=[readname, seq]

    max_lens={}
    for id in ids:
        max_lens[id]=[[],[]]  # R,L
        if len(R_read_lens[id]) >= 1:
            max_lens[id][0]=R_read_lens[id]
        if len(L_read_lens[id]) >= 1:
            max_lens[id][1]=L_read_lens[id]

    # group duplicated seqs
    seqs={}
    for id in max_lens:
        if len(max_lens[id][0]) == 2:
            h=id +'//'+ max_lens[id][0][0]
            if not max_lens[id][0][1] in seqs:
                seqs[max_lens[id][0][1]]='>'
            seqs[max_lens[id][0][1]] += h +';'
        if len(max_lens[id][1]) == 2:
            h=id +'//'+ max_lens[id][1][0]
            if not max_lens[id][1][1] in seqs:
                seqs[max_lens[id][1][1]]='>'
            seqs[max_lens[id][1][1]] += h +';'

    # output fasta
    with open(filenames.mapped_fa_select, 'w') as outfile:
        for seq in seqs:
            tmp=seqs[seq] +'\n'+ seq +'\n'
            outfile.write(tmp)
            pass


def blastn_for_mapped(args, params, q_path, db_path, outfpath):
    NcbiblastnCommandline(db=db_path, query=q_path, evalue=params.blastn_evalue_for_mapped, perc_identity=params.blastn_ident_for_mapped, word_size=params.blastn_word_size_for_mapped, num_threads=args.p, culling_limit=2, outfmt=6, out=outfpath)()


def pairing(params, filenames):
    def check_true_pos(list):
        retain=False
        if len(list) == 1:
            retain=True
        elif (float(list[0][2]) > float(list[1][2])) and (float(list[0][2]) >= params.mapped_abs_single_ident_threshold) and (float(list[1][2]) < params.mapped_abs_single_ident_threshold):
            retain=True
        if retain is True:
            ls=list[0]
            headers=ls[0].split(';')[:-1]
            for header in headers:
                id,h=header.split('//')
                chr,tmp=h.split(':', 1)
                start,tmp=tmp.split('-', 1)
                end,tmp=tmp.split('/', 1)
                if (chr == ls[1]) and (start == str(int(ls[8]) - 1)) and (end == ls[9]):  # 0-based
                    return id, h
        return None, None

    # load blastn result, find singleton hits
    currentq='any'
    header_singleton={}
    tmp=[]
    with open(filenames.blast4_res) as infile:
        for line in infile:
            ls=line.split()
            if ls[0] == currentq:
                tmp.append(ls)
            elif len(tmp) >= 1:
                id,header=check_true_pos(tmp)
                if not header == None:
                    if not id in header_singleton:
                        header_singleton[id]=[]
                    header_singleton[id].append(header)
                tmp=[ls]
                currentq=ls[0]
            else:
                tmp=[ls]
                currentq=ls[0]
    id,header=check_true_pos(tmp)
    if not header == None:
        if not id in header_singleton:
            header_singleton[id]=[]
        header_singleton[id].append(header)

    with open(filenames.bp_pair_single, 'w') as outfile:
        with open(filenames.breakpoint_pairs) as infile:
            for line in infile:
                ls=line.split()
                id=ls[0] +':'+ ls[1] +'-'+ ls[2] +'/'+ ls[7]
                if id in header_singleton:
                    outfile.write(line)

