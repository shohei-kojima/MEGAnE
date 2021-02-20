#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,datetime,collections,gzip
import log,traceback


def fill_in_ins(args, params, filenames):
    log.logger.debug('started')
    try:
        # load scaffold
        geno_orig={}
        bps={}
        with gzip.open(filenames.merged_vcf_ins) as infile:
            for line in infile:
                line=line.decode()
                if line[0] == '#':
                    if '#CHROM' in line:
                        hs=line.strip().split('\t')
                else:
                    ls=line.strip().split('\t')
                    for info in ls[7].split(';'):
                        if 'MEI=' in info:
                            mes=info.replace('MEI=', '')
                        elif '1END=' in info:
                            end=info.replace('1END=', '')
                            break
                    for me in mes.split('|'):
                        if not me in bps:
                            bps[me]={}
                        if not ls[0] in bps[me]:
                            bps[me][ls[0]]={}
                        for p in range(int(ls[1]), int(end) + 1):
                            bps[me][ls[0]][p]=ls[2]
                    geno_orig[ls[2]]={}
                    for h,v in zip(hs[9:], ls[9:]):
                        if not v == '0/0':
                            geno_orig[ls[2]][h]=v
        sample_ids=hs[9:]
        sample_ids_n=len(sample_ids)
        # fill-in
        processed_n=0
        for sample_id in sample_ids:
            dir=args.sample_id_to_dir[sample_id]
            f=args.dirs[dir][2]
            reads=collections.Counter()
            with gzip.open(f) as infile:
                for line in infile:
                    ls=line.decode().strip().split('\t')
                    if float(ls[2]) < params.overhang_evalue_threshold:
                        chr,tmp=ls[0].split(':', 1)
                        start,tmp=tmp.split('-', 1)
                        end,lr,_=tmp.split('/', 2)
                        bp=int(start) if lr == 'L' else int(end)
                        me,_=ls[1].split(',', 1)
                        if me in bps:
                            if chr in bps[me]:
                                if bp in bps[me][chr]:
                                    id=bps[me][chr][bp]
                                    if not sample_id in geno_orig[id]:
                                        reads[id] += 1
            for id in reads:
                if reads[id] >= params.min_support_reads_ins:
                    geno_orig[id][sample_id]='0/.'
            processed_n += 1
            if (processed_n % params.processed_interval) == 0:
                log.logger.info('%d samples processed...' % processed_n)
        
        # output
        missing_line_added=False
        with gzip.open(filenames.filled_vcf_ins, 'wt') as outfile:
            with gzip.open(filenames.merged_vcf_ins) as infile:
                for line in infile:
                    line=line.decode()
                    if line[0] == '#':
                        if '##FILTER=<ID=' in line and missing_line_added is False:
                            missing_line='##FILTER=<ID=M,Description="Too many missing genotypes">\n'
                            outfile.write(missing_line)
                            missing_line_added=True
                        outfile.write(line)
                    else:
                        ls=line.split('\t', 10)
                        tmp=ls[:9]
                        zero=0
                        missing=0
                        ac=0
                        for sample_id in sample_ids:
                            if sample_id in geno_orig[ls[2]]:
                                v=geno_orig[ls[2]][sample_id]
                                tmp.append(v)
                                if v == '0/.':
                                    missing += 1
                                elif v == '1/1':
                                    ac += 2
                                else:
                                    ac += 1
                            else:
                                tmp.append('0/0')
                                zero += 1
                        change_to_nonpass=False
                        if (ac / (sample_ids_n * 2)) < 0.05:
                            if missing >= (ac * 2):
                                change_to_nonpass=True
                        elif missing >= ac or missing > (zero / 2):
                            change_to_nonpass=True
                        if change_to_nonpass is True:
                            if tmp[6] == 'PASS':
                                tmp[6]='M'
                            else:
                                tmp[6]='%s;M' % tmp[6]
                        outfile.write('\t'.join(tmp) +'\n')
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def fill_in_abs(args, params, filenames):
    log.logger.debug('started')
    try:
        slop_len=params.slop_len_for_abs
        # load scaffold
        geno_orig={}
        bps={}
        with gzip.open(filenames.merged_vcf_abs) as infile:
            for line in infile:
                line=line.decode()
                if line[0] == '#':
                    if '#CHROM' in line:
                        hs=line.strip().split('\t')
                else:
                    ls=line.strip().split('\t')
                    for info in ls[7].split(';'):
                        if '1END=' in info:
                            end=info.replace('1END=', '')
                            break
                    if not ls[0] in bps:
                        bps[ls[0]]={}
                    for p in range(int(ls[1]) - slop_len, int(ls[1]) + slop_len):
                        bps[ls[0]][p]=ls[2]
                    for p in range(int(end) - slop_len, int(end) + slop_len):
                        bps[ls[0]][p]=ls[2]
                    geno_orig[ls[2]]={}
                    for h,v in zip(hs[9:], ls[9:]):
                        if not v == '0/0':
                            geno_orig[ls[2]][h]=v
        
        sample_ids=hs[9:]
        sample_ids_n=len(sample_ids)
        # fill-in
        processed_n=0
        for sample_id in sample_ids:
            dir=args.sample_id_to_dir[sample_id]
            f=args.dirs[dir][1]
            to_be_added=set()
            with gzip.open(f) as infile:
                for line in infile:
                    ls=line.decode().strip().split('\t')
                    chr=ls[1]
                    start=int(ls[2])
                    end=int(ls[3])
                    if chr in bps:
                        if start in bps[chr]:
                            if end in bps[chr]:
                                if bps[chr][start] == bps[chr][end]:
                                    to_be_added.add(bps[chr][start])
            for id in to_be_added:
                if not sample_id in geno_orig[id]:
                    geno_orig[id][sample_id]='0/.'
            processed_n += 1
            if (processed_n % params.processed_interval) == 0:
                log.logger.info('%d samples processed...' % processed_n)
        
        # output
        missing_line_added=False
        with gzip.open(filenames.filled_vcf_abs, 'wt') as outfile:
            with gzip.open(filenames.merged_vcf_abs) as infile:
                for line in infile:
                    line=line.decode()
                    if line[0] == '#':
                        if '##FILTER=<ID=' in line and missing_line_added is False:
                            missing_line='##FILTER=<ID=M,Description="Too many missing genotypes">\n'
                            missing_line_added=True
                        outfile.write(line)
                    else:
                        ls=line.split('\t', 10)
                        tmp=ls[:9]
                        zero=0
                        missing=0
                        ac=0
                        for sample_id in sample_ids:
                            if sample_id in geno_orig[ls[2]]:
                                v=geno_orig[ls[2]][sample_id]
                                tmp.append(v)
                                if v == '0/.':
                                    missing += 1
                                elif v == '1/1':
                                    ac += 2
                                else:
                                    ac += 1
                            else:
                                tmp.append('0/0')
                                zero += 1
                        change_to_nonpass=False
                        if (ac / (sample_ids_n * 2)) < 0.05:
                            if missing >= (ac * 2):
                                change_to_nonpass=True
                        elif ((ac + missing) / (sample_ids_n * 2)) < 0.75:
                            if missing >= ac or missing > zero:
                                change_to_nonpass=True
                        if change_to_nonpass is True:
                            if tmp[6] == 'PASS':
                                tmp[6]='M'
                            else:
                                tmp[6]='%s;M' % tmp[6]
                        outfile.write('\t'.join(tmp) +'\n')
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
