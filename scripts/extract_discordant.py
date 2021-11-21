#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,pysam,itertools,math,shutil,string
import log,traceback


nt=('A', 'T', 'G', 'C')
strands=('+', '-')

cigar_op={'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}
cigar_ref_retain={'M', 'D', 'N', '=', 'X'}
cigar_read_retain={'M', 'I', '=', 'X'}


# flagstat
def flagstat(args):
    log.logger.debug('started')
    try:
        flag=pysam.flagstat(args.b, '-@ %d' % (args.p - 1))  # multi process
        count= int(flag.split()[0])
        interval= math.ceil(count / args.p)
        return count,interval
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


# concatenate result files
def concat_for_ins(args, params, filenames):
    log.logger.debug('started')
    try:
        outfiles=[filenames.overhang_fa, filenames.overhang_pA, filenames.distant_txt, filenames.mapped_fa, filenames.blast1_res]
        for file_base in outfiles:
            shutil.move(file_base +'0.txt', file_base)
            with open(file_base, 'a') as outfile:
                for n in range(1, args.p):
                    with open(file_base + str(n) +'.txt') as infile:
                        shutil.copyfileobj(infile, outfile)
                    os.remove(file_base + str(n) +'.txt')
                outfile.flush()
                os.fdatasync(outfile.fileno())
        # unmapped read, remove short seq
        file_base=filenames.unmapped_fa
        unmapped_min_len= params.blastn_word_size + 2
        discarded=0
        with open(file_base, 'w') as outfile:
            for n in range(args.p):
                with open(file_base + str(n) +'.txt') as infile:
                    for line in infile:
                        if line[0] == '>':
                            tmp=line
                        else:
                            nonN_len=len(line.replace('N', ''))
                            if nonN_len >= unmapped_min_len:
                                tmp += line
                                outfile.write(tmp)
                            else:
                                discarded += 1
                os.remove(file_base + str(n) +'.txt')
            outfile.flush()
            os.fdatasync(outfile.fileno())
        log.logger.debug('%d unmapped reads discarded due to shorter than word size' % discarded)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def concat_for_abs(args, filenames):
    log.logger.debug('started')
    try:
        outfiles=[filenames.abs_txt]
        for file_base in outfiles:
            shutil.move(file_base +'0.txt', file_base)
            with open(file_base, 'a') as outfile:
                for n in range(1, args.p):
                    with open(file_base + str(n) +'.txt') as infile:
                        shutil.copyfileobj(infile, outfile)
                    os.remove(file_base + str(n) +'.txt')
                outfile.flush()
                os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def summarize(res):
    log.logger.debug('started')
    try:
        nonXY_n,X_n,Y_n,chimeric_n,hybrid_n,pA_n,unmapped_n,absent_n=0,0,0,0,0,0,0,0
        for nonXY,X,Y,chimeric,hybrid,pA,unmapped,absent in res:
            nonXY_n    += nonXY
            X_n        += X
            Y_n        += Y
            chimeric_n += chimeric
            hybrid_n   += hybrid
            pA_n       += pA
            unmapped_n += unmapped
            absent_n   += absent
        log.logger.info('Screening results:nonXY_reads=%d,X_reads=%d,Y_reads=%s,chimeric_reads=%d,hybrid_reads=%d,pA_reads=%d,unmapped_reads=%d,absent_reads=%d' % (nonXY_n, X_n, Y_n, chimeric_n, hybrid_n, pA_n, unmapped_n, absent_n))
        return [nonXY_n, X_n, Y_n]
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


# main
def main(args, params, filenames, n):
    log.logger.debug('started')
    try:
        do_ins=False if args.only_abs is True else True
        do_abs=False if args.only_ins is True else True
        
        if not n is None:  # multi process
            if do_ins is True:
                f_overhang  =open(filenames.overhang_fa + str(n) + '.txt', 'w')
                f_pA        =open(filenames.overhang_pA + str(n) + '.txt', 'w')
                f_distant   =open(filenames.distant_txt + str(n) + '.txt', 'w')
                f_unmapped  =open(filenames.unmapped_fa + str(n) + '.txt', 'w')
                f_mapped    =open(filenames.mapped_fa   + str(n) + '.txt', 'w')
            if do_abs is True:
                f_abs       =open(filenames.abs_txt     + str(n) + '.txt', 'w')
        else:  # single process
            if do_ins is True:
                f_overhang  =open(filenames.overhang_fa, 'w')
                f_pA        =open(filenames.overhang_pA, 'w')
                f_distant   =open(filenames.distant_txt, 'w')
                f_unmapped  =open(filenames.unmapped_fa, 'w')
                f_mapped    =open(filenames.mapped_fa  , 'w')
            if do_abs is True:
                f_abs       =open(filenames.abs_txt    , 'w')

        def complement(string):
            seq_c=string.translate(str.maketrans('ATGC', 'TACG'))[::-1]
            return seq_c

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

        if not args.b is None:
            infile=pysam.AlignmentFile(args.b, 'rb')
        elif not args.c is None:
            infile=pysam.AlignmentFile(args.c, 'rc', reference_filename=args.fa)
        if not n is None:  # multi process
            infile=itertools.islice(infile, n, None, args.p)
        chrX_set=params.chrX
        chrY_set=params.chrY
        nonXY_n=0
        X_n=0
        Y_n=0
        chimeric_n=0
        hybrid_n=0
        pA_n=0
        unmapped_n=0
        absent_n=0
        for line in infile:
            line=line.tostring()
            ls=line.strip().split('\t')
            if int(ls[1]) < 2048:  # remove supplementary alignment
                b=bin(int(ls[1]))
                if b[-1] == '1':   # paired-end
                    if b[-3] == '1':   # output unmapped
                        if do_ins is True:
                            if len(b) >= 12:
                                if b[-10] == '1':  # read fails platform/vendor quality checks
                                    continue
                            strand='/1' if b[-7] == '1' else '/2'
                            f_unmapped.write('>%s%s\n%s\n' % (ls[0], strand, ls[9]))
                            unmapped_n += 1
                    else:
                        if ls[2] in chrX_set:
                            X_n += 1
                        elif ls[2] in chrY_set:
                            Y_n += 1
                        else:
                            nonXY_n += 1
                        if do_ins is True:
                            if (('S' in ls[5]) or ('SA:Z:' in line) or ('XA:Z:' in line)) and not ('H' in ls[5]):   # start retrieving overhangs
                                deletion=False
                                fseq=ls[9].upper()
                                strand='+'
                                if b[-5] == '1':
                                    strand='-'
                                    fseq=complement(fseq)
                                first_or_second='/1' if b[-7] == '1' else '/2'
                                elems=[]
                                if ls[2] in args.main_chrs_set:
                                    if 'S' in ls[5]:
                                        elems.append([ls[2], ls[3], ls[5], strand])  # chr, pos, cigar, strand
                                if 'SA:Z:' in line:
                                    for l in ls[11:]:
                                        if 'SA:Z:' in l:
                                            lsp=l.replace('SA:Z:', '').split(';')
                                            for i in lsp[:-1]:
                                                isp=i.split(',')
                                                if isp[0] in args.main_chrs_set:
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
                                                    if isp[0] in args.main_chrs_set:
                                                        elems.append([isp[0], isp[1][1:], isp[2], isp[1][0]])
                                                break
                                    d={}
                                    d_mapped={}
                                    rseq=complement(fseq)
                                    for chr,pos,cigar,strand in elems:
                                        breakpoint,L_clip_len,R_clip_len=determine_breakpoint_from_cigar(cigar)
                                        if not breakpoint == 'NA':
                                            if max(L_clip_len, R_clip_len) >= params.discordant_reads_clip_len:
                                                seq=fseq if strand == '+' else rseq
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
                                                    mapped_seq= seq[L_clip_len:seqlen - R_clip_len]
                                                    if (Acount / len(clip_seq)) >= params.polyA_overhang_threshold:
                                                        f_pA.write('%s:%d-%d/%s/%s%s/%s\t%d\n' % (chr, start, end, breakpoint, ls[0], first_or_second, strand, len(clip_seq)))  # reads with pA
                                                        pA_n += 1
                                                        if not mapped_seq in d_mapped:
                                                            d_mapped[mapped_seq]=[]
                                                        d_mapped[mapped_seq].append('%s:%d-%d/%s/%s%s/%s' % (chr, start, end, breakpoint, ls[0], first_or_second, strand))
                                                    else:
                                                        if not clip_seq in d:
                                                            d[clip_seq]=[]
                                                        if not mapped_seq in d_mapped:
                                                            d_mapped[mapped_seq]=[]
                                                        h='%s:%d-%d/%s/%s%s/%s' % (chr, start, end, breakpoint, ls[0], first_or_second, strand)
                                                        d[clip_seq].append(h)
                                                        d_mapped[mapped_seq].append(h)
                                    if len(d) >= 1:
                                        out_overhang=[ '>%s;\n%s\n' % (';'.join(d[seq]), seq) for seq in d ]
                                        f_overhang.write(''.join(out_overhang))   # end retrieving overhang seqs
                                        chimeric_n += 1
                                    if len(d_mapped) >= 1:
                                        out_mapped=[ '>%s;\n%s\n' % (';'.join(d_mapped[seq]), seq) for seq in d_mapped ]
                                        f_mapped.write(''.join(out_mapped))   # end retrieving mapped seqs
                            
                            ins=int(ls[8])
                            if (ins == 0) or (ins <= -params.read_pair_gap_len) or (params.read_pair_gap_len <= ins):    # start retrieving distant reads
                                if b[-4:-2] == '00':   # both mapped
                                    retain=False
                                    if not 'S' in ls[5]:
                                        retain=True
                                    else:
                                        breakpoint,left,right=determine_breakpoint_from_cigar(ls[5])
                                        if breakpoint == 'NA':
                                            retain=True
                                        elif max(left, right) < params.discordant_reads_clip_len:
                                            retain=True
                                    if retain is True:
                                        tmp=[]
                                        if ls[2] in args.main_chrs_set:
                                            dir='+' if b[-5] == '0' else '-'
                                            start= int(ls[3]) - 1   # 0-based
                                            end= start + calc_ref_len(ls[5])
                                            tmp.append('%s:%d-%d/%s' % (ls[2], start, end, dir))
                                        if 'XA:Z:' in line:
                                            for l in ls[11:]:
                                                if 'XA:Z:' in l:
                                                    xs=l.replace('XA:Z:', '').split(';')[:-1]
                                                    break
                                            for x in xs:
                                                chr,pos,xcigar,_=x.split(',')
                                                retain=False
                                                if not 'S' in xcigar:
                                                    retain=True
                                                else:
                                                    breakpoint,left,right=determine_breakpoint_from_cigar(ls[5])
                                                    if breakpoint == 'NA':
                                                        retain=True
                                                    elif max(left, right) < params.discordant_reads_clip_len:
                                                        retain=True
                                                if retain is True:
                                                    if chr in args.main_chrs_set:
                                                        start= int(pos[1:]) - 1   # 0-based
                                                        end= start + calc_ref_len(xcigar)
                                                        tmp.append('%s:%d-%d/%s' % (chr, start, end, pos[0]))
                                        if len(tmp) >= 1:
                                            first_or_second='/1' if b[-7] == '1' else '/2'
                                            readname= ls[0] + first_or_second
                                            f_distant.write(readname +'\t'+ ';'.join(tmp) +'\n')   # end retrieving distant reads
                                            hybrid_n += 1
                        
                        if do_abs is True:   # retrieve reads with absent ME
                            if 'SA:Z:' in line:
                                first_or_second='/1' if b[-7] == '1' else '/2'
                                saz,xaz={},{}
                                saz['+'],xaz['+'],saz['-'],xaz['-']={},{},{},{}
                                if ('S' in ls[5]) and not ('H' in ls[5]):
                                    if ls[2] in args.main_chrs_set:
                                        breakpoint,l_len,r_len=determine_breakpoint_from_cigar(ls[5])
                                        if not breakpoint == 'NA':
                                            clip_len=l_len if breakpoint == 'L' else r_len
                                            if clip_len > params.discordant_reads_clip_len:
                                                strand='+' if b[-5] == 0 else '-'
                                                start= int(ls[3]) - 1  # 0-based
                                                end= start + calc_ref_len(ls[5])
                                                cigar_reshape= '%d-%d' % (l_len, l_len + calc_read_len(ls[5]))
                                                if not ls[2] in xaz[strand]:
                                                    xaz[strand][ls[2]]={}
                                                if breakpoint == 'L':
                                                    if not 'L' in xaz[strand][ls[2]]:
                                                        xaz[strand][ls[2]]['L']=[]
                                                    xaz[strand][ls[2]]['L'].append((start, end, cigar_reshape))
                                                else:
                                                    if not 'R' in xaz[strand][ls[2]]:
                                                        xaz[strand][ls[2]]['R']=[]
                                                    xaz[strand][ls[2]]['R'].append((start, end, cigar_reshape))
                                for l in ls[11:]:
                                    if 'SA:Z:' in l:
                                        ss=l.replace('SA:Z:', '').split(';')[:-1]
                                        for i in ss:
                                            ispl=i.split(',')
                                            if 'S' in ispl[3]:
                                                if ispl[0] in args.main_chrs_set:
                                                    breakpoint,l_len,r_len=determine_breakpoint_from_cigar(ispl[3])
                                                    if not breakpoint == 'NA':
                                                        clip_len=l_len if breakpoint == 'L' else r_len
                                                        if clip_len > params.discordant_reads_clip_len:
                                                            start= int(ispl[1]) - 1  # 0-based
                                                            end= start + calc_ref_len(ispl[3])
                                                            cigar_reshape= '%d-%d' % (l_len, l_len + calc_read_len(ispl[3]))
                                                            if not ispl[0] in saz[ispl[2]]:
                                                                saz[ispl[2]][ispl[0]]={}
                                                            if breakpoint == 'L':
                                                                if not 'L' in saz[ispl[2]][ispl[0]]:
                                                                    saz[ispl[2]][ispl[0]]['L']=[]
                                                                saz[ispl[2]][ispl[0]]['L'].append((start, end, cigar_reshape))
                                                            else:
                                                                if not 'R' in saz[ispl[2]][ispl[0]]:
                                                                    saz[ispl[2]][ispl[0]]['R']=[]
                                                                saz[ispl[2]][ispl[0]]['R'].append((start, end, cigar_reshape))
                                        break
                                for l in ls[11:]:
                                    if 'XA:Z:' in l:
                                        ss=l.replace('XA:Z:', '').split(';')[:-1]
                                        for i in ss:
                                            ispl=i.split(',')
                                            if 'S' in ispl[2]:
                                                if ispl[0] in args.main_chrs_set:
                                                    breakpoint,l_len,r_len=determine_breakpoint_from_cigar(ispl[2])
                                                    if not breakpoint == 'NA':
                                                        clip_len=l_len if breakpoint == 'L' else r_len
                                                        if clip_len > params.discordant_reads_clip_len:
                                                            start= int(ispl[1][1:]) - 1  # 0-based
                                                            end= start + calc_ref_len(ispl[2])
                                                            cigar_reshape= '%d-%d' % (l_len, l_len + calc_read_len(ispl[2]))
                                                            if not ispl[0] in xaz[ispl[1][0]]:
                                                                xaz[ispl[1][0]][ispl[0]]={}
                                                            if breakpoint == 'L':
                                                                if not 'L' in xaz[ispl[1][0]][ispl[0]]:
                                                                    xaz[ispl[1][0]][ispl[0]]['L']=[]
                                                                xaz[ispl[1][0]][ispl[0]]['L'].append((start, end, cigar_reshape))
                                                            else:
                                                                if not 'R' in xaz[ispl[1][0]][ispl[0]]:
                                                                    xaz[ispl[1][0]][ispl[0]]['R']=[]
                                                                xaz[ispl[1][0]][ispl[0]]['R'].append((start, end, cigar_reshape))
                                        break
                                for strand in strands:
                                    for chr in saz[strand]:
                                        if chr in xaz[strand]:
                                            if ('R' in saz[strand][chr]) and ('L' in xaz[strand][chr]):
                                                for r_s,r_e,r_pos in saz[strand][chr]['R']:
                                                    for l_s,l_e,l_pos in xaz[strand][chr]['L']:
                                                        if params.abs_min_dist <= (l_s - r_e) <= params.abs_max_dist:
                                                            f_abs.write('%s%s\t%s\t%d\t%d\t%s:%d-%d\t%s:%d-%d\t%s\t%s\n' %(ls[0], first_or_second, chr, r_e, l_s, chr, r_s, r_e, chr, l_s, l_e, r_pos, l_pos))
                                                            absent_n += 1
                                            if ('R' in xaz[strand][chr]) and ('L' in saz[strand][chr]):
                                                for r_s,r_e,r_pos in xaz[strand][chr]['R']:
                                                    for l_s,l_e,l_pos in saz[strand][chr]['L']:
                                                        if params.abs_min_dist <= (l_s - r_e) <= params.abs_max_dist:
                                                            f_abs.write('%s%s\t%s\t%d\t%d\t%s:%d-%d\t%s:%d-%d\t%s\t%s\n' %(ls[0], first_or_second, chr, r_e, l_s, chr, r_s, r_e, chr, l_s, l_e, r_pos, l_pos))   # end retrieving reads with absent ME
                                                            absent_n += 1
        if do_ins is True:
            f_overhang.flush()
            f_pA.flush()
            f_distant.flush()
            f_unmapped.flush()
            f_mapped.flush()
            os.fdatasync(f_overhang.fileno())
            os.fdatasync(f_pA.fileno())
            os.fdatasync(f_distant.fileno())
            os.fdatasync(f_unmapped.fileno())
            os.fdatasync(f_mapped.fileno())
            f_overhang.close()
            f_pA.close()
            f_distant.close()
            f_unmapped.close()
            f_mapped.close()
        if do_abs is True:
            f_abs.flush()
            os.fdatasync(f_abs.fileno())
            f_abs.close()
        return [nonXY_n, X_n, Y_n, chimeric_n, hybrid_n, pA_n, unmapped_n, absent_n]
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

