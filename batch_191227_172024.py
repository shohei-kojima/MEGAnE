#!/usr/bin/env python

"""
2019/12/27 17:20:24 +0900

# batch_191227_172024.py
# generated at /home/kooojiii/results/misc/1000_hgp/aimup
# usage: python %prog
# python3.7

"""




def cigar_retain(string):
    length=0
    tmp=''
    for c in string:
        if not c in cigar_op:
            tmp += c
        elif c in cigar_op_retain:
            length += int(tmp)
            tmp=''
        else:
            tmp=''
    return length

def count_clip(cigar):
    length=0
    tmp=''
    for c in cigar:
        if not c in cigar_op:
            tmp += c
        elif c in cigar_op_clip:
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
        elif c in cigar_op:
            if c in cigar_op_clip:
                left += int(tmp)
            break
    tmp=''
    if string[-1] in cigar_op_clip:
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
    return breakpoint


strands=('+', '-')
for line in fileinput.input():
    ls=line.strip().split('\t')
    if 'SA:Z:' in line:
        saz,xaz={},{}
        saz['+'],xaz['+'],saz['-'],xaz['-']={},{},{},{}
        if 'S' in ls[5]:
            breakpoint,l_len,r_len=determine_breakpoint_from_cigar(ls[5])
            if not breakpoint == 'NA':
                clip_len=l_len if breakpoint == 'L' else r_len
                if clip_len > params.clip_len:
                    strand='+' if b[-5] == 0 else '-'
                    start= int(ls[3]) - 1  # 0-based
                    end= start + calc_ref_len(ls[5])
                    if not ls[2] in xaz[strand]:
                        xaz[strand][ls[2]]={}
                    if breakpoint == 'L':
                        if not 'L' in xaz[strand][ls[2]]:
                            xaz[strand][ls[2]]['L']=[]
                        xaz[strand][ls[2]]['L'].append((start, end))
                    else:
                        if not 'R' in xaz[strand][ls[2]]:
                            xaz[strand][ls[2]]['R']=[]
                        xaz[strand][ls[2]]['R'].append((start, end))
        for l in ls[::-1]:
            if 'SA:Z:' in l:
                ss=l.replace('SA:Z:', '').split(';')[:-1]
                for i in ss:
                    ispl=i.split(',')
                    if 'S' in ispl[3]:
                        breakpoint,l_len,r_len=determine_breakpoint_from_cigar(ispl[3])
                        if not breakpoint == 'NA':
                            clip_len=l_len if breakpoint == 'L' else r_len
                            if clip_len > params.clip_len:
                                start= int(ispl[1]) - 1  # 0-based
                                end= start + calc_ref_len(ispl[3])
                                if not ispl[0] in saz[ispl[2]]:
                                    saz[ispl[2]][ispl[0]]={}
                                if breakpoint == 'L':
                                    if not 'L' in saz[ispl[2]][ispl[0]]:
                                        saz[ispl[2]][ispl[0]]['L']=[]
                                    saz[ispl[2]][ispl[0]]['L'].append((start, end))
                                else:
                                    if not 'R' in saz[ispl[2]][ispl[0]]:
                                        saz[ispl[2]][ispl[0]]['R']=[]
                                    saz[ispl[2]][ispl[0]]['R'].append((start, end))
                break
        for l in ls[::-1]:
            if 'XA:Z:' in l:
                ss=l.replace('XA:Z:', '').split(';')[:-1]
                for i in ss:
                    ispl=i.split(',')
                    if 'S' in ispl[2]:
                        breakpoint,l_len,r_len=determine_breakpoint_from_cigar(ispl[2])
                        if not breakpoint == 'NA':
                            clip_len=l_len if breakpoint == 'L' else r_len
                            if clip_len > params.clip_len:
                                start= int(ispl[1][1:]) - 1  # 0-based
                                end= start + calc_ref_len(ispl[2])
                                if not ispl[0] in xaz[ispl[1][0]]:
                                    xaz[ispl[1][0]][ispl[0]]={}
                                if breakpoint == 'L':
                                    if not 'L' in xaz[ispl[1][0]][ispl[0]]:
                                        xaz[ispl[1][0]][ispl[0]]['L']=[]
                                    xaz[ispl[1][0]][ispl[0]]['L'].append((start, end))
                                else:
                                    if not 'R' in xaz[ispl[1][0]][ispl[0]]:
                                        xaz[ispl[1][0]][ispl[0]]['R']=[]
                                    xaz[ispl[1][0]][ispl[0]]['R'].append((start, end))
                break
        for strand in strands:
            for chr in saz[strand]:
                if chr in xaz[strand]:
                    if (len(saz[strand][chr]['R']) >= 1) and (len(xaz[strand][chr]['L']) >= 1):
                        for r_s,r_e in saz[strand][chr]['R']:
                            for l_s,l_e in xaz[strand][chr]['L']:
                                if params.abs_min_dist <= (l_s - r_e) <= params.abs_max_dist:
                                    f_abs.write('%s\t%s\t%d\t%d\t%s:%d-%d\t%s:%d-%d\n' %(ls[0], chr, r_e, l_s, chr, r_s, r_e, chr, l_s, l_e))
                    if (len(saz[strand][chr]['L']) >= 1) and (len(xaz[strand][chr]['R']) >= 1):
                        for r_s,r_e in xaz[strand][chr]['R']:
                            for l_s,l_e in saz[strand][chr]['L']:
                                if params.abs_min_dist <= (l_s - r_e) <= params.abs_max_dist:
                                    f_abs.write('%s\t%s\t%d\t%d\t%s:%d-%d\t%s:%d-%d\n' %(ls[0], chr, r_e, l_s, chr, r_s, r_e, chr, l_s, l_e))

