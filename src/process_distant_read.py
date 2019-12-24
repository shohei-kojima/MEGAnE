#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

def process_reads(args, params, distant_txt, outfpath):
    
    def process_bed(bed):
        ids={}
        bed=BedTool(bed, from_string=True)
        bed=bed.intersect(simple, v=True)
        bed=bed.intersect(tes, f=0.9, wa=True, wb=True)
        if len(bed) >= 1:
            for line in bed:
                line=str(line)
                ls=line.split()
                te=ls[7].split(':')[0]
                if not ls[3] in ids:
                    ids[ls[3]]=set()
                ids[ls[3]].add(te)
        return ids

    # read repbase file, pick up simple repeats
    name_to_clas={}
    with open(args.rep) as infile:
        for line in infile:
            if '>' in line:
                ls=line.strip().replace('>', '').split('\t')
                if len(ls) >= 2:
                    clas=ls[1]
                    if ls[0] == 'SVA2':
                        clas='SINE'
                elif ('Alu' in ls[0]) or ('FLA' in ls[0]):
                    clas='SINE1/7SL'
                elif 'LTR26' in ls[0]:
                    clas='ERV1'
                name_to_clas[ls[0]]=clas

    # list up simple repeat in the ref genome
    simple=''
    tes=''
    with open(args.repout) as infile:
        for line in infile:
            ls=line.split()
            name,clas=ls[3].split(':')
            if clas == 'Simple_repeat':
                simple += line
            elif name_to_clas[name] in simple_rep:
                simple += line
            else:
                tes += line
            pass
    simple=BedTool(simple, from_string=True)
    tes=BedTool(tes, from_string=True)

    # pairing of distant reads
    d={}
    with open(distant_txt) as infile:
        for line in infile:
            ls=line.split()
            id,dir=ls[0].split('/')
            if not id in d:
                d[id]=set()
            d[id].add(dir)

    retain=set()
    for id in d:
        if len(d[id]) == 2:
            retain.add(id)

    # extract reads intersect with ref TEs
    d={}
    bed=[]
    with open(distant_txt) as infile:
        for line in infile:
            ls=line.split()
            rn,dir=ls[0].split('/')
            if rn in retain:
                id,poss=line.split()
                for pos in poss.split(';'):
                    tmp,_=pos.rsplit('/', 1)
                    tmp,end=tmp.rsplit('-', 1)
                    chr,start=tmp.rsplit(':', 1)
                    bed.append('\t'.join([chr, start, end, id]))
                if len(bed) >= 100000:  # chunk
                    bed='\n'.join(bed)
                    tmp=process_bed(bed)
                    d.update(tmp)
                    bed=[]
    del(simple)
    del(tes)

    # pairing, 2nd
    pair_d={}
    for fid in d:
        id,first_or_second=fid.split('/')
        if not id in pair_d:
            pair_d[id]={}
        if not first_or_second in pair_d[id]:
            pair_d[id][first_or_second]=d[fid]

    retain={}
    for id in pair_d:
        if len(pair_d[id]) == 1:
            key=list(pair_d[id].keys())[0]
            first_or_second='1' if key == '2' else '2'
            retain[id +'/'+ first_or_second]=sorted(list(pair_d[id][key]))
        elif not len(pair_d[id]['1'] & pair_d[id]['2']) >= 1:
            keys=list(pair_d[id].keys())
            for key in keys:
                first_or_second='1' if key == '2' else '2'
                retain[id +'/'+ first_or_second]=sorted(list(pair_d[id][key]))

    # output bed
    with open(outfpath, 'w') as outfile:
        with open(distant_txt) as infile:
            for line in infile:
                ls=line.split()
                if ls[0] in retain:
                    te=';'.join(retain[ls[0]])
                    for pos in ls[1].split(';'):
                        tmp,strand=pos.rsplit('/', 1)
                        tmp,end=tmp.rsplit('-', 1)
                        chr,start=tmp.rsplit(':', 1)
                        outfile.write('\t'.join([chr, start, end, te, ls[0], strand]) +'\n')
