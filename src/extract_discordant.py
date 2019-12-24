#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,pysam,itertools,math,shutil,cython
from Bio.Seq import Seq


nt=['A', 'T', 'G', 'C']

cigar_op=['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
cigar_ref_retain=['M', 'D', 'N', '=', 'X']
cigar_read_retain=['M', 'I', '=', 'X']


# flagstat
def flagstat(args):
    flag=pysam.flagstat(args.b, '-@ %d' % args.p)  # multi process
    count= int(flag.split()[0])
    interval= math.ceil(count / args.p)
    return count,interval


# concatenate result files
def concat(args, outfiles):
    for file_base in outfiles:
        shutil.move(file_base +'0.txt', file_base)
        with open(file_base, 'a') as outfile:
            for n in range(1, args.p):
                with open(file_base + str(n) +'.txt') as infile:
                    shutil.copyfileobj(infile, outfile)
                os.remove(file_base + str(n) +'.txt')


# main
def main(args, params, main_chrs, outfiles, n, count, interval):
    if not n is None:  # multi process
        start= interval * n
        end= min(interval * (n + 1), count)
        f_overhang  =open(outfiles[0] + str(n) + '.txt', 'w')
        f_pA        =open(outfiles[1] + str(n) + '.txt', 'w')
        f_distant   =open(outfiles[2] + str(n) + '.txt', 'w')
        f_unmapped  =open(outfiles[3] + str(n) + '.txt', 'w')
    else:  # single process
        f_overhang  =open(outfiles[0], 'w')
        f_pA        =open(outfiles[1], 'w')
        f_distant   =open(outfiles[2], 'w')
        f_unmapped  =open(outfiles[3], 'w')

    def complement(string):
        c=Seq(string).reverse_complement()
        return str(c)

    def count_clip(cigar):
        length=0
        tmp=''
        for c in cigar:
            if not c in cigar_op:
                tmp += c
            elif c == 'S':
                length += int(tmp)
                tmp=''
            else:
                tmp=''
        return length

    def determine_breakpoint_from_cigar(string):
        left,right=0,0
        tmp=''
        for c in string:
            if not c in cigar_op:
                tmp += c
            elif c == 'S':
                left += int(tmp)
                break
            else:
                break
        tmp=''
        if string[-1] == 'S':
            for c in string[-2::-1]:
                if not c in cigar_op:
                    tmp = c + tmp
                else:
                    right += int(tmp)
                    break
        if left > right:
            breakpoint='L'
        elif left == right:
            breakpoint='NA'
        else:
            breakpoint='R'
        return breakpoint,left,right

    def calc_ref_len(cigar):
        length=0
        tmp=''
        for c in cigar:
            if not c in cigar_op:
                tmp += c
            elif c in cigar_ref_retain:
                length += int(tmp)
                tmp=''
            else:
                tmp=''
        return length

    def calc_read_len(cigar):
        length=0
        tmp=''
        for c in cigar:
            if not c in cigar_op:
                tmp += c
            elif c in cigar_read_retain:
                length += int(tmp)
                tmp=''
            else:
                tmp=''
        return length

    def check_repeat(seq):
        judge=False
        seqlen=len(seq)
        for n in nt:
            c=seq.count(n)
            if (c / seqlen) >= params.mapped_region_low_complex_threshold:
                judge=True
                break
        return judge

    
    with pysam.AlignmentFile(args.b, 'rb') as infile:
        if not n is None:  # multi process
            infile=itertools.islice(infile, start, end, 1)
        for line in infile:
            line=line.tostring()
            ls=line.strip().split('\t')
            if int(ls[1]) < 2048:  # remove supplementary alignment
                b=bin(int(ls[1]))
                if b[-1] == '1':   # paired-end
                    if b[-3] == '1':   # retrieving unmapped
                        strand='/1' if b[-7] == '1' else '/2'
                        f_unmapped.write('>'+ ls[0] + strand +'\n'+ ls[9] +'\n')
                    else:
                        if ('S' in ls[5]) or ('SA:Z:' in line) or ('XA:Z:' in line):   # start retrieving overhangs
                            deletion=False
                            fseq=ls[9]
                            strand='+'
                            if len(b) >= 7:
                                if b[-5] == '1':
                                    strand='-'
                                    fseq=complement(fseq)
                            if len(b) >= 9:
                                first_or_second='/1' if b[-7] == '1' else '/2'
                            elems=[]
                            if ls[2] in main_chrs:
                                if 'S' in ls[5]:
                                    elems.append([ls[2], ls[3], ls[5], strand])  # chr, pos, cigar, strand
                            if 'SA:Z:' in line:
                                for l in ls[11:]:
                                    if 'SA:Z:' in l:
                                        lsp=l.replace('SA:Z:', '').split(';')
                                        for i in lsp[:-1]:
                                            isp=i.split(',')
                                            if isp[0] in main_chrs:
                                                elems.append([isp[0], isp[1], isp[3], isp[2]])
                                        break
                                for chr,pos,_,_ in elems[1:]:
                                    if chr == ls[2]:
                                        if abs(int(ls[3]) - int(pos)) < params.max_TSD_len:
                                            deletion=True
                                            break
                            if deletion is False:
                                if 'XA:Z:' in line:
                                    for l in ls[11:]:
                                        if 'XA:Z:' in l:
                                            lsp=l.replace('XA:Z:', '').split(';')
                                            for i in lsp[:-1]:
                                                isp=i.split(',')
                                                if isp[0] in main_chrs:
                                                    elems.append([isp[0], isp[1][1:], isp[2], isp[1][0]])
                                            break
                                d={}
                                rseq=complement(fseq)
                                for chr,pos,cigar,strand in elems:
                                    breakpoint,L_clip_len,R_clip_len=determine_breakpoint_from_cigar(cigar)
                                    if not breakpoint == 'NA':
                                        if max(L_clip_len, R_clip_len) >= params.discordant_reads_clip_len:
                                            seq=fseq
                                            if strand == '-':
                                                seq=rseq
                                            seqlen=len(seq)
                                            read_seq= seq[L_clip_len:seqlen - R_clip_len]
                                            rep=check_repeat(read_seq)
                                            if rep is False:
                                                ref_len=calc_ref_len(cigar)
                                                start= int(pos) - 1  # 0-based
                                                end= start + ref_len
                                                if breakpoint == 'L':
                                                    clip_seq= seq[0:L_clip_len]
                                                    Acount=clip_seq.count('A')
                                                else:
                                                    clip_seq= seq[seqlen - R_clip_len:]
                                                    Acount=clip_seq.count('T')
                                                if (Acount / len(clip_seq)) >= params.polyA_overhang_threshold:
                                                    f_pA.write(chr +':'+ str(start) +'-'+ str(end) +'/'+ breakpoint +'/'+ ls[0] + first_or_second +'/'+ strand +'\t'+ str(len(clip_seq)) +'\n')  # reads with pA
                                                else:
                                                    if not clip_seq in d:
                                                        d[clip_seq]=[]
                                                    h= chr +':'+ str(start) +'-'+ str(end) +'/'+ breakpoint +'/'+ ls[0] + first_or_second +'/'+ strand
                                                    d[clip_seq].append(h)
                                if len(d) >= 1:
                                    out_overhang=''
                                    for seq in d:
                                        out_overhang += '>'
                                        for i in d[seq]:
                                            out_overhang += i +';'
                                    out_overhang += '\n'+ seq +'\n'
                                    f_overhang.write(out_overhang)   # end retrieving overhangs
                        ins=int(ls[8])
                        if (ins == 0) or (ins <= -params.read_pair_gap_len) or (params.read_pair_gap_len <= ins):    # start retrieving distant reads
                            if b[-4:-2] == '00':
                                if not 'S' in ls[5]:
                                    tmp=[]
                                    for l in ls[::-1]:
                                        if 'MC:Z:' in l:
                                            mcigar=l.replace('MC:Z:', '')
                                            break
                                    if len(b) >= 9:
                                        first_or_second='/1' if b[-7] == '1' else '/2'
                                    readname= ls[0] + first_or_second
                                    if not 'S' in mcigar:
                                        if len(b) >= 7:
                                            dir='+' if b[-5] == '0' else '-'
                                        start= int(ls[3]) - 1  # 0-based
                                        end= start + calc_ref_len(ls[5])
                                        tmp.append(ls[2] +':'+ str(start) +'-'+ str(end) +'/'+ dir)
                                    if 'XA:Z:' in line:
                                        for l in ls[::-1]:
                                            if 'XA:Z:' in l:
                                                xs=l.replace('XA:Z:', '').split(';')[:-1]
                                                break
                                        for x in xs:
                                            chr,pos,xcigar,_=x.split(',')
                                            if not 'S' in xcigar:
                                                dir=pos[0]
                                                start=pos[1:]
                                                end= int(start) + calc_ref_len(xcigar)
                                                tmp.append(chr +':'+ start +'-'+ str(end) +'/'+ dir)
                                    if len(tmp) >= 1:
                                        f_distant.write(readname +'\t'+ ';'.join(tmp) +'\n')   # end retrieving distant reads
    f_overhang.close()
    f_pA.close()
    f_distant.close()
    f_unmapped.close()


