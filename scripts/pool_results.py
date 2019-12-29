#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


from pybedtools import BedTool


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
        for n in range(len(list[0]) - 8):
            s=''
            n=n+8
            for ls in list:
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
                bed_d[ls[7]]=''
            bed_d[ls[7]] += line
    with open(filenames.bp_info_single) as infile:
        for line in infile:
            if not line in first_clean:
                ls=line.split()
                if not ls[7] in bed_d:
                    bed_d[ls[7]]=''
                bed_d[ls[7]] += line
    lines={}
    for m in bed_d:
        bed=BedTool(bed_d[m], from_string=True).sort()
        for line in bed:
            line=str(line)
            ls=line.split()
            if not ls[7] in lines:
                lines[ls[7]]={}
            if not ls[0] in lines[ls[7]]:
                lines[ls[7]][ls[0]]=''
            lines[ls[7]][ls[0]] += line
    del(bed_d)

    with open(filenames.bp_merged, 'w') as outfile:
        for clas in lines:
            for chr in lines[clas]:
                bed=BedTool(lines[clas][chr], from_string=True)
                bed=bed.sort()
                prev_end=-1
                prev_line=''
                tmp=[]
                for line in bed:
                    line=str(line)
                    ls=line.split()
                    start=int(ls[1])
                    end=int(ls[2])
                    if start <= prev_end:
                        tmp.append(prev_ls)
                    else:
                        if len(tmp) >= 1:
                            tmp.append(prev_ls)
                            outfile.write(merge(tmp))
                            tmp=[]
                        else:
                            outfile.write(prev_line)
                    prev_end=end
                    prev_line=line
                    prev_ls=ls
                if len(tmp) >= 1:
                    tmp.append(prev_ls)
                    outfile.write(merge(tmp))
                    tmp=[]
                else:
                    outfile.write(prev_line)


def add_hybrid(params, filenames):
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
                            te.add(i)
            te=';'.join(list(te))
            id=ls[0] +':'+ ls[1] +'-'+ ls[2] +'/'+ ls[7]
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

    # output
    with open(filenames.bp_merged_all, 'w') as outfile:
        with open(filenames.bp_merged) as infile:
            for line in infile:
                ls=line.split()
                id=ls[0] +':'+ ls[1] +'-'+ ls[2] +'/'+ ls[7]
                if id in hybrid:
                    r,r_num=ls[3].split(':')
                    l,l_num=ls[4].split(':')
                    r_num,l_num=int(r_num),int(l_num)
                    r_hub_num=len(hybrid[id]['up']) if 'up' in hybrid[id] else 0
                    l_hub_num=len(hybrid[id]['down']) if 'down' in hybrid[id] else 0
                    r_num += r_hub_num
                    l_num += l_hub_num
                    s=[ls[0], ls[1], ls[2], r +':'+ str(r_num), l +':'+ str(l_num)]
                    s.extend(ls[5:])
                    s.extend([str(r_hub_num), str(l_hub_num)])
                else:
                    s=[]
                    s.extend(ls)
                    s.extend(['0', '0'])
                outline= '\t'.join(s) +'\n'
                outfile.write(outline)
