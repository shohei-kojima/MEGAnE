#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os
from pybedtools import BedTool
import log,traceback


def merge_breakpoints(filenames):
    def merge(list):
        chr=list[0][0]
        te=list[0][7]
        # ls[3]
        ss=set()
        for ls in list:
            ss.add(ls[3])
        n,pos=[],[]
        for s in ss:
            n.append(int(s.split(':')[1]))
            pos.append(int(s.split(':')[0]))
        max_n=0
        total=0
        true_pos=[]
        for num,p in zip(n,pos):
            total += num
            if num > max_n:
                max_n=num
        for num,p in zip(n,pos):
            if num == max_n:
                true_pos.append(p)
        start=max(true_pos)
        r_pos_num= str(max(true_pos)) +':'+ str(total)
        # 4
        ss=set()
        for ls in list:
            ss.add(ls[4])
        n,pos=[],[]
        for s in ss:
            n.append(int(s.split(':')[1]))
            pos.append(int(s.split(':')[0]))
        max_n=0
        total=0
        true_pos=[]
        for num,p in zip(n,pos):
            total += num
            if num > max_n:
                max_n=num
        for num,p in zip(n,pos):
            if num == max_n:
                true_pos.append(p)
        end=min(true_pos)
        l_pos_num= str(min(true_pos)) +':'+ str(total)
        # 5
        ss=set()
        for ls in list:
            ss.add(ls[3].split(':')[0] +':'+ ls[5])
        n,not_NA=0,False
        for s in ss:
            if not 'NA' in s:
                n += int(s.split(':')[1])
                not_NA=True
        if not_NA is False:
            n='NA'
        else:
            n=str(n)
        r_pA= n
        # 6
        ss=set()
        for ls in list:
            ss.add(ls[4].split(':')[0] +':'+ ls[6])
        n,not_NA=0,False
        for s in ss:
            if not 'NA' in s:
                n += int(s.split(':')[1])
                not_NA=True
        if not_NA is False:
            n='NA'
        else:
            n=str(n)
        l_pA= n
        tmp=[]
        for n in [8, 9]:
            s=''
            for ls in list:
                s += ls[n] +';'
            tmp.append(s)
        for n in range(10, len(list[0])):
            s=''
            for ls in list:
                if ls[n] == 'NA':
                    s += ls[n] +';;'
                else:
                    s += ls[n] +';'
            tmp.append(s)
        # output
        if start > end:
            start,end=str(end),str(start)
        else:
            start,end=str(start),str(end)
        l=[chr, start, end, r_pos_num, l_pos_num, r_pA, l_pA, te]
        l.extend(tmp)
        ret_str='\t'.join(l)
        return ret_str +'\n'

    # merge results, delete redundant
    bed_d={}
    first_clean=set()
    with open(filenames.breakpoint_clean) as infile:
        for line in infile:
            first_clean.add(line)
            ls=line.split()
            if not ls[7] in bed_d:
                bed_d[ls[7]]=[]
            bed_d[ls[7]].append(ls)
    with open(filenames.bp_info_single) as infile:
        for line in infile:
            if not line in first_clean:
                ls=line.split()
                if not ls[7] in bed_d:
                    bed_d[ls[7]]=[]
                bed_d[ls[7]].append(ls)
    lines={}
    for m in bed_d:
        lines[m]={}
        bed_d[m]=sorted(bed_d[m], key=lambda x:(x[0], int(x[1]), int(x[2])))
        for ls in bed_d[m]:
            if not ls[0] in lines[m]:
                lines[m][ls[0]]=[]
            lines[m][ls[0]].append(ls)
    del(bed_d)

    with open(filenames.bp_merged, 'w') as outfile:
        for m in lines:
            for chr in lines[m]:
                prev_end=-1
                prev_line=''
                tmp=[]
                for ls in lines[m][chr]:
                    start=int(ls[1])
                    if start <= prev_end:
                        tmp.append(prev_ls)
                    else:
                        if len(tmp) >= 1:
                            tmp.append(prev_ls)
                            outfile.write(merge(tmp))
                            tmp=[]
                        else:
                            outfile.write(prev_line)
                    prev_end=int(ls[2])
                    prev_line='\t'.join(ls) +'\n'
                    prev_ls=ls
                if len(tmp) >= 1:
                    tmp.append(prev_ls)
                    outfile.write(merge(tmp))
                    tmp=[]
                else:
                    outfile.write(prev_line)
        outfile.flush()
        os.fdatasync(outfile.fileno())


def add_hybrid(params, filenames):
    log.logger.debug('started')
    try:
        bed_up,bed_down='',''
        with open(filenames.bp_merged) as infile:
            for line in infile:
                ls=line.split()
                r=int(ls[3].split(':')[0])
                l=int(ls[4].split(':')[0])
                s= r - params.hybrid_read_range_from_breakpint
                e= l + params.hybrid_read_range_from_breakpint
                s=0 if s < 0 else s
                te=set()
                for t in ls[10:12]:
                    if not t == 'NA':
                        for i in t.split(';'):
                            if not i == '':
                                te.add(i.split(',')[0])
                te=';'.join(list(te))
                id=ls[0] +':'+ ls[1] +'-'+ ls[2] +'/'+ ls[7] +'/ID='+ ls[12]
                bed_up += ls[0] +'\t'+ str(s) +'\t'+ str(r) +'\t'+ id +'\t'+ te +'\n'
                bed_down += ls[0] +'\t'+ str(l) +'\t'+ str(e) +'\t'+ id +'\t'+ te +'\n'
        bed_up=BedTool(bed_up, from_string=True)
        bed_down=BedTool(bed_down, from_string=True)

        # load hybrid reads
        bed_hybrid=BedTool(filenames.hybrid)

        # list up hybrid reads near from integration junction candidates
        hybrid={}
        bed_up=bed_up.intersect(bed_hybrid, wa=True, wb=True)
        for line in bed_up:
            line=str(line)
            ls=line.split()
            if ls[10] == '+':
                chimeric_te=set(ls[4].split(';'))
                hybrid_te=set(ls[8].split(';'))
                shared_te= chimeric_te & hybrid_te
                if len(shared_te) >= 1:
                    if not ls[3] in hybrid:
                        hybrid[ls[3]]={}
                    if not 'up' in hybrid[ls[3]]:
                        hybrid[ls[3]]['up']=set()
                    hybrid_read_info= ls[5] +':'+ ls[6] +'-'+ ls[7] +'(+)'
                    hybrid[ls[3]]['up'].add(ls[9] +'\t'+ ls[8] +'\t'+ hybrid_read_info)

        bed_down=bed_down.intersect(bed_hybrid, wa=True, wb=True)
        for line in bed_down:
            line=str(line)
            ls=line.split()
            if ls[10] == '-':
                chimeric_te=set(ls[4].split(';'))
                hybrid_te=set(ls[8].split(';'))
                shared_te= chimeric_te & hybrid_te
                if len(shared_te) >= 1:
                    if not ls[3] in hybrid:
                        hybrid[ls[3]]={}
                    if not 'down' in hybrid[ls[3]]:
                        hybrid[ls[3]]['down']=set()
                    hybrid_read_info= ls[5] +':'+ ls[6] +'-'+ ls[7] +'(-)'
                    hybrid[ls[3]]['down'].add(ls[9] +'\t'+ ls[8] +'\t'+ hybrid_read_info)

        # write hybrid read info to a file, master file
        with open(filenames.hybrid_master, 'w') as outfile:
            for pos in hybrid:
                if 'up' in hybrid[pos]:
                    for line in hybrid[pos]['up']:
                        rname,te,info=line.split()
                        outfile.write(pos +'\t'+ rname +'\t'+ info +'\t'+ te +'\n')
                if 'down' in hybrid[pos]:
                    for line in hybrid[pos]['down']:
                        rname,te,info=line.split()
                        outfile.write(pos +'\t'+ rname +'\t'+ info +'\t'+ te +'\n')
            outfile.flush()
            os.fdatasync(outfile.fileno())

        # output
        with open(filenames.bp_merged_all, 'w') as outfile:
            with open(filenames.bp_merged) as infile:
                for line in infile:
                    ls=line.split()
                    id=ls[0] +':'+ ls[1] +'-'+ ls[2] +'/'+ ls[7] +'/ID='+ ls[12]
                    if id in hybrid:
                        r,r_num=ls[3].split(':')
                        l,l_num=ls[4].split(':')
                        r_num,l_num=int(r_num),int(l_num)
                        r_hub_num=len(hybrid[id]['up']) if 'up' in hybrid[id] else 0
                        l_hub_num=len(hybrid[id]['down']) if 'down' in hybrid[id] else 0
                        r_num += r_hub_num
                        l_num += l_hub_num
                        s=[ls[0], ls[1], ls[2], r +':'+ str(r_num), l +':'+ str(l_num)]
                        s.extend(ls[5:12])
                        s.extend([str(r_hub_num), str(l_hub_num)])
                        s.append(ls[12])
                    else:
                        s=[]
                        s.extend(ls[:12])
                        s.extend(['0', '0'])
                        s.append(ls[12])
                    outline= '\t'.join(s) +'\n'
                    outfile.write(outline)
            outfile.flush()
            os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
