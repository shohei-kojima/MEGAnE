#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,datetime,collections,gzip
from multiprocessing import Pool
import mmap,io
import numpy as np
import pybedtools
from pybedtools import BedTool
from utils import parse_fasta_for_merge_vcf, load_me_classification
import log,traceback


def rename_L1(me, me_to_class):
    if ('_5end' in me) or ('_3end' in me) or ('_orf2' in me):
        if me_to_class[me] == 'LINE/L1':
            me=me.replace('_5end', '').replace('_3end', '').replace('_orf2', '')
    return me


def judge_split_L1(me, me_to_class):
    if ('_5end' in me) or ('_3end' in me) or ('_orf2' in me):
        if me_to_class[me] == 'LINE/L1':
            return True
    return False


def check_paths(args):
    log.logger.debug('started')
    try:
        if args.input_scaffold is None:
            if args.merge_mei is True:
                necessary_files=['MEI_final_gaussian_genotyped.vcf', 'breakpoint_pairs_pooled_all.txt.gz', 'overhang_to_MEI_list.txt.gz']
            elif args.merge_absent_me is True:
                necessary_files=['absent_MEs_genotyped.vcf', 'absent.txt.gz']
            len_necessary_files=len(necessary_files)
        dirs={}
        n=0
        if args.input_scaffold is None:
            with open(args.f) as infile:
                for line in infile:
                    n += 1
                    dir=line.strip()
                    if os.path.exists(dir) is False:
                        log.logger.error('%s does not exist.' % dir)
                        exit(1)
                    tmp=[]
                    for f in necessary_files:
                        f=os.path.join(dir, f)
                        if os.path.exists(f) is False:
                            log.logger.warning('%s did not found. This sample will NOT be merged.' % f)
                            continue
                        tmp.append(f)
                    if len(tmp) == len_necessary_files:
                        dirs[dir]=tmp
            log.logger.info('%d directories specified. %d directories will be analyzed.' % (n, len(dirs)))
        else:
            sample_id_to_dir={}
            with open(args.chunk_f) as infile:
                for line in infile:
                    n += 1
                    ls=line.split()
                    overhang_to_MEI_f=ls[1]
                    dir=os.path.dirname(overhang_to_MEI_f)
                    if args.merge_mei is True:
                        tmp=[0, 0]
                    elif args.merge_absent_me is True:
                        tmp=[0]
                    if os.path.exists(overhang_to_MEI_f) is False:
                        log.logger.error('%s did not found. Please check the path again.' % overhang_to_MEI_f)
                        exit(1)
                    tmp.append(overhang_to_MEI_f)
                    tmp.append(ls[0])
                    dirs[dir]=tmp
                    sample_id_to_dir[ls[0]]=dir
            log.logger.info('%d directories specified. %d directories will be analyzed.' % (n, len(dirs)))
            args.sample_id_to_dir=sample_id_to_dir
        
        if not n == len(dirs):
            failed_n= n - len(dirs)
            log.logger.warning('MEGAnE did not find necessary files in %d directory(ies). These directory(ies) will not be used for this analysis.' % failed_n)
        args.dirs=dirs
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def load_geno_ins_mmap(args, dir, te_len, chrY_set, me_to_class):
    f=args.dirs[dir][0]
    vcf_loaded=set()
    tmp_bed={}
    tmp_vlc={}
    with open(f) as infile:
        with mmap.mmap(infile.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
            for line in iter(mapped.readline, b''):
                line=line.decode()
                if line[0] == '#':
                    if line[:6] == '#CHROM':
                        ls=line.split()
                        sample_id=ls[9]
                elif not line[0] == '#':
                    ls=line.strip().split('\t')
                    if ls[0] in args.chr and not 'Y' in ls[6]:
                        te=ls[7].split(';')[0].replace('SVTYPE=', '')
                        start=str(int(ls[1]) - 1)
                        end=ls[7].split(';')[3].replace('MEI_rpos=', '')
                        pred_class=[]
                        te_lens=[]
                        strands=[]
                        infos=ls[7].split(';')
                        if infos[1] == 'MEPRED=PASS':
                            mepred='PASS'
                            if 'MEI=' in infos[4]:
                                preds=infos[4].replace('MEI=', '')
                                for pred in preds.split('|'):
                                    if len(pred.split(',')) == 3:
                                        clas,bp,strand=pred.split(',')
                                        pred_class.append('%s,%s' % (clas, clas))
                                        if 'NA' in bp:
                                            s,e=0,0
                                        else:
                                            s,e=bp.split('/')
                                        if strand == '+/+':
                                            strands.append('+')
                                            tmp_len= int(e) - int(s)
                                            if tmp_len > 0:
                                                te_lens.append(tmp_len)
                                        elif strand == '-/-':
                                            strands.append('-')
                                            tmp_len= int(s) - int(e)
                                            if tmp_len > 0:
                                                te_lens.append(tmp_len)
                                        else:
                                            strands.append(strand)
                            elif infos[4] == 'MEI_lbp=pT':
                                preds=infos[5].replace('MEI_rbp=', '')
                                for pred in preds.split('|'):
                                    clas,bp,strand=pred.split(',')
                                    _clas=rename_L1(clas, me_to_class)
                                    pred_class.append('%s,%s' % (_clas, _clas))
                                    if judge_split_L1(clas, me_to_class) is False:
                                        te_lens.append(te_len[clas] - int(bp))
                                    strands.append(strand)
                            elif infos[5] == 'MEI_rbp=pA':
                                preds=infos[4].replace('MEI_lbp=', '')
                                for pred in preds.split('|'):
                                    _clas,bp,strand=pred.split(',')
                                    clas=rename_L1(clas, me_to_class)
                                    pred_class.append('%s,%s' % (_clas, _clas))
                                    if judge_split_L1(clas, me_to_class) is False:
                                        te_lens.append(te_len[clas] - int(bp))
                                    strands.append(strand)
                        else:
                            mepred='FAILED'
                            left, right=infos[4].replace('MEI_lbp=', ''), infos[5].replace('MEI_rbp=', '')
                            for pred in left.split('|'):
                                pred_class.append(rename_L1(pred.split(',')[0], me_to_class))
                            for pred in right.split('|'):
                                pred_class.append(rename_L1(pred.split(',')[0], me_to_class))
                        if len(te_lens) > 0:
                            telen= str(int(np.round(np.mean(te_lens))))
                        else:
                            telen= 'NA'
                        if len(strands) > 0:
                            strand=collections.Counter(strands).most_common()[0][0]
                        else:
                            strand= 'NA'
                        if not te in tmp_bed:
                            tmp_bed[te]=[]
                        if len(pred_class) > 0:
                            tmp_bed[te].append('\t'.join([ls[0], start, end, '%s,%s,%s,%s,%s,%s' % (sample_id, ls[9], mepred, telen, ls[6], strand), ','.join(pred_class)]))
                            vcf_loaded.add(ls[2])
    # load very low confidence
    f=args.dirs[dir][1]
    with open(f) as infile0:
        with mmap.mmap(infile0.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
            with io.TextIOWrapper(gzip.GzipFile(fileobj=mapped)) as infile:
                for line in infile:
                    ls=line.strip().split('\t')
                    if ls[0] in args.chr:
                        if not ls[-1] in vcf_loaded and not ls[0] in chrY_set:
                            if not ls[7] in tmp_vlc:
                                tmp_vlc[ls[7]]=[]
                            tmp_vlc[ls[7]].append('%s\t%s\t%s\t%s\n' % (ls[0], ls[1], ls[2], sample_id))
    return [tmp_bed, tmp_vlc, sample_id, dir]


def load_geno_ins(args, dir, te_len, chrY_set, me_to_class):
    f=args.dirs[dir][0]
    vcf_loaded=set()
    tmp_bed={}
    tmp_vlc={}
    with open(f) as infile:
        for line in infile:
            if line[0] == '#':
                if line[:6] == '#CHROM':
                    ls=line.split()
                    sample_id=ls[9]
            elif not line[0] == '#':
                ls=line.strip().split('\t')
                if ls[0] in args.chr and not 'Y' in ls[6]:
                    te=ls[7].split(';')[0].replace('SVTYPE=', '')
                    start=str(int(ls[1]) - 1)
                    end=ls[7].split(';')[3].replace('MEI_rpos=', '')
                    pred_class=[]
                    te_lens=[]
                    strands=[]
                    infos=ls[7].split(';')
                    if infos[1] == 'MEPRED=PASS':
                        mepred='PASS'
                        if 'MEI=' in infos[4]:
                            preds=infos[4].replace('MEI=', '')
                            for pred in preds.split('|'):
                                if len(pred.split(',')) == 3:
                                    clas,bp,strand=pred.split(',')
                                    pred_class.append('%s,%s' % (clas, clas))
                                    if 'NA' in bp:
                                        s,e=0,0
                                    else:
                                        s,e=bp.split('/')
                                    if strand == '+/+':
                                        strands.append('+')
                                        tmp_len= int(e) - int(s)
                                        if tmp_len > 0:
                                            te_lens.append(tmp_len)
                                    elif strand == '-/-':
                                        strands.append('-')
                                        tmp_len= int(s) - int(e)
                                        if tmp_len > 0:
                                            te_lens.append(tmp_len)
                                    else:
                                        strands.append(strand)
                        elif infos[4] == 'MEI_lbp=pT':
                            preds=infos[5].replace('MEI_rbp=', '')
                            for pred in preds.split('|'):
                                clas,bp,strand=pred.split(',')
                                _clas=rename_L1(clas, me_to_class)
                                pred_class.append('%s,%s' % (_clas, _clas))
                                if judge_split_L1(clas, me_to_class) is False:
                                    te_lens.append(te_len[clas] - int(bp))
                                strands.append(strand)
                        elif infos[5] == 'MEI_rbp=pA':
                            preds=infos[4].replace('MEI_lbp=', '')
                            for pred in preds.split('|'):
                                clas,bp,strand=pred.split(',')
                                _clas=rename_L1(clas, me_to_class)
                                pred_class.append('%s,%s' % (_clas, _clas))
                                if judge_split_L1(clas, me_to_class) is False:
                                    te_lens.append(te_len[clas] - int(bp))
                                strands.append(strand)
                    else:
                        mepred='FAILED'
                        left, right=infos[4].replace('MEI_lbp=', ''), infos[5].replace('MEI_rbp=', '')
                        for pred in left.split('|'):
                            pred_class.append(rename_L1(pred.split(',')[0], me_to_class))
                        for pred in right.split('|'):
                            pred_class.append(rename_L1(pred.split(',')[0], me_to_class))
                    if len(te_lens) > 0:
                        telen= str(int(np.round(np.mean(te_lens))))
                    else:
                        telen= 'NA'
                    if len(strands) > 0:
                        strand=collections.Counter(strands).most_common()[0][0]
                    else:
                        strand= 'NA'
                    if not te in tmp_bed:
                        tmp_bed[te]=[]
                    if len(pred_class) > 0:
                        tmp_bed[te].append('\t'.join([ls[0], start, end, '%s,%s,%s,%s,%s,%s' % (sample_id, ls[9], mepred, telen, ls[6], strand), ','.join(pred_class)]))
                        vcf_loaded.add(ls[2])
    # load very low confidence
    f=args.dirs[dir][1]
    with gzip.open(f, 'rt') as infile:
        for line in infile:
            ls=line.strip().split('\t')
            if ls[0] in args.chr:
                if not ls[-1] in vcf_loaded and not ls[0] in chrY_set:
                    if not ls[7] in tmp_vlc:
                        tmp_vlc[ls[7]]=[]
                    tmp_vlc[ls[7]].append('%s\t%s\t%s\t%s\n' % (ls[0], ls[1], ls[2], sample_id))
    return [tmp_bed, tmp_vlc, sample_id, dir]


def merge_vcf_ins(args, params, filenames):
    log.logger.debug('started')
    try:
        pybedtools.set_tempdir(args.pybedtools_tmp)
        chrY_set=set()
        if args.no_sex_chr is False:
            chrY_set=set([ chr for chr in args.male_sex_chr.split(',') ])
        else:
            pass
        me_to_class,_=load_me_classification(args.rep)
        
        # load te length
        te_len={}
        tmp=[]
        with open(args.rep) as infile:
            for line in infile:
                if '>' in line:
                    if len(tmp) > 0:
                        te_len[te]=len(''.join(tmp))
                    te=line.split()[0].replace('>', '')
                    tmp=[]
                else:
                    tmp.append(line.strip())
        te_len[te]=len(''.join(tmp))
        
        # load file paths
        file_num=len(args.dirs)
        
        # load high confidence and unique MEIs
        all={}
        if args.chr is None:
            chrs_set=set()
            for dir in args.dirs:
                f=args.dirs[dir][0]
                with open(f) as infile:
                    for line in infile:
                        if not line[0] == '#':
                            ls=line.split('\t')
                            if ls[6] == 'PASS':
                                chrs_set.add(ls[0])
                                te=ls[7].split(';')[0].replace('SVTYPE=', '')
                                pos=[ls[0], str(int(ls[1]) - 1), ls[7].split(';')[3].replace('MEI_rpos=', '')]
                                if not te in all:
                                    all[te]=[]
                                all[te].append(pos)
            log.logger.debug('%d chrs found.' % len(chrs_set))
            args.chr=chrs_set
        else:
            found_on_chr=False
            for dir in args.dirs:
                f=args.dirs[dir][0]
                with open(f) as infile:
                    for line in infile:
                        if not line[0] == '#':
                            ls=line.split('\t')
                            if ls[0] in args.chr:  # check chr
                                found_on_chr=True
                                if ls[6] == 'PASS':
                                    te=ls[7].split(';')[0].replace('SVTYPE=', '')
                                    pos=[ls[0], str(int(ls[1]) - 1), ls[7].split(';')[3].replace('MEI_rpos=', '')]
                                    if not te in all:
                                        all[te]=[]
                                    all[te].append(pos)
            if found_on_chr == False:
                log.logger.error('No MEs found in specified chromosomes. Exit.')
                exit(1)
        if len(all) == 0:
            log.logger.error('No MEs with "FILTER PASS" found. Exit.')
            exit(1)
        # header
        header=[]
        header.append('##fileformat=VCFv4.1\n')
        header.append('##fileDate=%s\n' % str(datetime.datetime.now()).split('.')[0])
        header.append('##source=MEI merge version "%s"\n' % args.version)
        header.append('##reference=%s\n' % args.fa)
        # load header
        for dir in args.dirs:
            f=args.dirs[dir][0]
            vcf_loaded=set()
            with open(f) as infile:
                for line in infile:
                    if line[0] == '#':
                        if '##contig=' in line or '##ALT=' in line or '##FILTER=' in line:
                            header.append(line)
                        if line[:6] == '#CHROM':
                            break
            break
        
        # count how many samples have a certain MEI
        bed={}
        vlc={}
        sample_id_to_dir={}
        if args.mmap is True:
            te_lens=[ te_len.copy() for _ in range(args.p) ]
            dirs=[]
            for dir in args.dirs:
                dirs.append(dir)
            print_chunk=20
            chunk=0
            with Pool(processes=args.p) as p:
                for n in range(0, file_num, args.p):
                    end= n + args.p
                    if end > file_num:
                        end=file_num
                    jobs=[]
                    for dir,te_len_each in zip(dirs[n:end], te_lens):
                        jobs.append(p.apply_async(load_geno_ins_mmap, (args, dir, te_len_each, chrY_set, me_to_class)))
                    res=[]
                    for job in jobs:
                        res.append(job.get())
                    for tmp_bed,tmp_vlc,sample_id,dir in res:
                        for me in tmp_bed:
                            if not me in bed:
                                bed[me]=[]
                            bed[me].extend(tmp_bed[me])
                        for me in tmp_vlc:
                            if not me in vlc:
                                vlc[me]=[]
                            vlc[me].extend(tmp_vlc[me])
                        sample_id_to_dir[sample_id]=dir
                        args.dirs[dir].append(sample_id)
                    chunk += 1
                    if (chunk % print_chunk) == 0:
                        log.logger.info('Loading genotypes, %d samples processed...' % (chunk * args.p))
            log.logger.info('Loading genotypes, %d samples processed...' % end)
        else:
            processed_n=0
            for dir in args.dirs:
                tmp_bed,tmp_vlc,sample_id,dir=load_geno_ins(args, dir, te_len, chrY_set, me_to_class)
                for me in tmp_bed:
                    if not me in bed:
                        bed[me]=[]
                    bed[me].extend(tmp_bed[me])
                for me in tmp_vlc:
                    if not me in vlc:
                        vlc[me]=[]
                    vlc[me].extend(tmp_vlc[me])
                sample_id_to_dir[sample_id]=dir
                args.dirs[dir].append(sample_id)
                processed_n += 1
                if (processed_n % 50) == 0:
                    log.logger.info('Loading genotypes, %d samples processed...' % processed_n)
            log.logger.info('Loading genotypes, %d samples processed...' % processed_n)
        
        args.sample_id_to_dir=sample_id_to_dir
        count={}
        poss={}
        vlc_intersects={}
        for m in bed:
            if m in all:
                if not m in count:
                    count[m]={}
                    poss[m]={}
                all[m]=sorted(all[m], key=lambda x:(x[0], int(x[1]), int(x[2])))
                all[m]=[ '\t'.join(l) for l in all[m] ]
                all[m]=BedTool('\n'.join(all[m]), from_string=True).merge()
                for line in all[m]:
                    ls=str(line).split()
                    id='%s:%s-%s:%s' % (ls[0], ls[1], ls[2], m)
                    count[m][id]=[]
                    poss[m][id]=[[], []]
                log.logger.info('Merging %s....' % m)
                
                tmp=BedTool('\n'.join(bed[m]), from_string=True)
                tmp=tmp.intersect(all[m], wa=True, wb=True, nonamecheck=True)
                ss,es=[],[]
                for line in tmp:
                    ls=str(line).split()
                    id='%s:%s-%s:%s' % (ls[5], ls[6], ls[7], m)
                    poss[m][id][0].append(ls[6])
                    poss[m][id][1].append(ls[7])
                    count[m][id].append([ls[3], ls[4]])
                 
                if m in vlc:
                    tmp=BedTool('\n'.join(vlc[m]), from_string=True)
                    tmp=tmp.intersect(all[m], wa=True, wb=True, nonamecheck=True)
                    for line in tmp:
                        ls=str(line).split()
                        id='%s:%s-%s:%s' % (ls[4], ls[5], ls[6], m)
                        if not id in vlc_intersects:
                            vlc_intersects[id]=set()
                        vlc_intersects[id].add(ls[3])
        # determine most common breakpoint
        failed={}
        for m in poss:
            failed[m]=set()
            for id in poss[m]:
                if len(poss[m][id][0]) >= 1:
                    start=collections.Counter(poss[m][id][0]).most_common()[0][0]
                    end=  collections.Counter(poss[m][id][1]).most_common()[0][0]
                    poss[m][id]=[start, end]
                else:
                    failed[m].add(id)
        for m in failed:
            log.logger.debug('%s=%d dropped,%d kept' % (m, len(failed[m]), len(poss[m])))
            if len(failed[m]) >= 1:
                log.logger.debug('%s=%d ids failed' % (m, len(failed[m])))
        # retrieve fasta
        bed=[]
        for m in count:
            for id in count[m]:
                if id in failed[m]:
                    continue
                chr=id.split(':')[0]
                start= int(poss[m][id][0]) - 1  # add one base before bp
                if start < 0:
                    start=0
                bed.append('%s\t%d\t%s\t%s\n' % (chr, start, poss[m][id][1], id))
                log.logger.debug('%s,%d,%s,%s' % (chr, start, poss[m][id][1], id))
        log.logger.debug('len of bed=%d' % len(bed))
        bed=BedTool(''.join(bed), from_string=True)
        log.logger.debug('bed file obj was made')
        fa=bed.sequence(fi=args.fa, name=True)
        log.logger.debug('fa was made')
        tmp=parse_fasta_for_merge_vcf(fa.seqfn)
        log.logger.debug('len of tmp=%d' % len(tmp))
        fa={}
        for h in tmp:
            fa[h.split('::')[0].replace('>', '')]=tmp[h].upper()
        log.logger.debug('len of fa=%d' % len(fa))
        # output
        n=0
        d={}
        sample_ids=set()
        for m in count:
            log.logger.debug('processing=%s' % m)
            for id in count[m]:
                if id in failed[m]:
                    continue
                afs=[]    # allele freq
                quals=[]  # me prediction quality
                lens=[]   # me length
                filts=[]  # me filter
                strands=[]  # me stands
                for info,_ in count[m][id]:
                    infos=info.split(',')
                    afs.append(int(infos[1][0]))  # in case '1.5' -> '1'
                    quals.append(infos[2])
                    if not infos[3] == 'NA':
                        lens.append(int(infos[3]))
                    filts.append(infos[4])
                    strands.append(infos[5])
                ac= sum(afs)
                pass_perc= quals.count('PASS') / len(quals)
                mepred='PASS' if pass_perc >= params.pass_perc_threshold else 'LC'  # flag if MEPRED=FAILED was majority
                pass_perc= filts.count('PASS') / len(filts)
                if pass_perc > params.pass_perc_threshold:
                    filt=collections.Counter(filts).most_common()[0][0]
                else:
                    filts_no_pass=[ v for v in filts if not v == 'PASS' ]
                    filt=collections.Counter(filts_no_pass).most_common()[0][0]
                if (len(lens) / len(filts)) > params.pass_perc_threshold:
                    melen=str(int(np.round(np.mean(lens))))
                else:
                    melen='NA'
                if len(strands) > 0:
                    strand=collections.Counter(strands).most_common()[0][0]
                else:
                    strand='NA'
                preds=[]  # te class pred
                preds=','.join([ l[1] for l in count[m][id] ])
                preds=preds.split(',')
                c=collections.Counter(preds).most_common()
                max_num=c[0][1]
                pred_class=[]
                for t in c:
                    if t[1] == max_num:
                        pred_class.append(t[0])
                    elif t[1] < max_num:
                        break
                pred_class='|'.join(pred_class)
                d[id]={}
                d[id]['CHROM']=id.split(':')[0]
                d[id]['POS']=poss[m][id][0]  # 0-based to 1-based; add one base before bp
                d[id]['ID']='%s_ins_%d' % (args.cohort_name, n)
                d[id]['REF']=fa['%s' % id]
                d[id]['ALT']='<INS:ME>'
                d[id]['QUAL']='.'
                d[id]['FILTER']=filt
                d[id]['INFO']='SVTYPE=%s;MEI=%s;MEPRED=%s;0START=%s;0END=%s;SVLEN=%s;MESTRAND=%s;AC=%d' % (m, pred_class, mepred, poss[m][id][0], poss[m][id][1], melen, strand, ac)
                d[id]['FORMAT']='GT'
                for info,_ in count[m][id]:
                    sample_id,c,_,_,_,_=info.split(',')
                    if c == '0':
                        gt='0/0'
                    elif c == '1':
                        gt='0/1'
                    elif c == '1.5':
                        gt='./1'
                    else:
                        gt='1/1'
                    d[id][sample_id]=gt
                    sample_ids.add(sample_id)
                n += 1
            log.logger.debug('processed=%s,%d' % (m, n))
        log.logger.debug('total processed=%d' % n)
        sample_ids=sorted(list(sample_ids))
        for id in d:
            for sample_id in sample_ids:
                if not sample_id in d[id]:
                    d[id][sample_id]='0/0'
                    if id in vlc_intersects:
                        if sample_id in vlc_intersects[id]:
                            d[id][sample_id]='0/.'
        log.logger.debug('genotyping done')
        header.append('##INFO=<ID=SVTYPE,Number=.,Type=String,Description="Mobile element class">\n')
        header.append('##INFO=<ID=MEI,Number=.,Type=String,Description="Mobile element subclass">\n')
        header.append('##INFO=<ID=MEPRED,Number=.,Type=String,Description="MEI prediction status">\n')
        header.append('##INFO=<ID=0START,Number=1,Type=Integer,Description="Position of left breakpoint. 0-based numbering (eg BED style).">\n')
        header.append('##INFO=<ID=0END,Number=1,Type=Integer,Description="Position of right breakpoint. 0-based numbering (eg BED style).">\n')
        header.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
        header.append('##INFO=<ID=MESTRAND,Number=.,Type=String,Description="Strandness of MEI">\n')
        header.append('##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count.">\n')
        header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header.append('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
        log.logger.debug('just before outputting filenames.merged_vcf_ins')
        with gzip.open(filenames.merged_vcf_ins, 'wt') as outfile:
            outfile.write(''.join(header))
            outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % '\t'.join(sample_ids))
            col_names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + sample_ids
            for id in d:
                tmp=[]
                for c in col_names:
                    tmp.append(d[id][c])
                outfile.write('\t'.join(tmp) +'\n')
        log.logger.debug('all done')
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def merge_vcf_abs(args, params, filenames):
    log.logger.debug('started')
    try:
        pybedtools.set_tempdir(args.pybedtools_tmp)
        # load te length
        clas_to_fam={}
        tmp=[]
        with open(filenames.reshaped_rep) as infile:
            for line in infile:
                if '>' in line:
                    ls=line.split()
                    clas=ls[0].replace('>', '')
                    clas_to_fam[clas]=ls[1]
        
        # load file paths
        file_num=len(args.dirs)

        # load high confidence and unique MEIs
        all={}
        if args.chr is None:
            chrs_set=set()
            for dir in args.dirs:
                f=args.dirs[dir][0]
                with open(f) as infile:
                    for line in infile:
                        if not line[0] == '#':
                            ls=line.split()
                            if ls[6] == 'PASS':
                                chrs_set.add(ls[0])
                                tmp=ls[7].split(';SVLEN=')[0].replace('MEI=', '')
                                fam_set=set()
                                for v in tmp.split(';'):
                                    clas,_=v.split(':')
                                    if clas in clas_to_fam:
                                        fam=clas_to_fam[clas]
                                        fam_set.add(fam)
                                leng=ls[7].split(';SVLEN=')[1].split(';')[0]
                                start= int(ls[1]) - 1
                                end= start + int(leng)
                                pos=[ls[0], str(start), str(end)]
                                if len(fam_set) >= 1:
                                    te='|'.join(sorted(list(fam_set)))
                                    if not te in all:
                                        all[te]=[]
                                    all[te].append(pos)
            log.logger.debug('%d chrs found.' % len(chrs_set))
            args.chr=chrs_set
        else:
            found_on_chr=False
            for dir in args.dirs:
                f=args.dirs[dir][0]
                with open(f) as infile:
                    for line in infile:
                        if not line[0] == '#':
                            ls=line.split('\t')
                            if ls[0] in args.chr:
                                found_on_chr=True
                                if ls[6] == 'PASS':
                                    tmp=ls[7].split(';SVLEN=')[0].replace('MEI=', '')
                                    fam_set=set()
                                    for v in tmp.split(';'):
                                        clas,_=v.split(':')
                                        if clas in clas_to_fam:
                                            fam=clas_to_fam[clas]
                                            fam_set.add(fam)
                                    leng=ls[7].split(';SVLEN=')[1].split(';')[0]
                                    start= int(ls[1]) - 1
                                    end= start + int(leng)
                                    pos=[ls[0], str(start), str(end)]
                                    if len(fam_set) >= 1:
                                        te='|'.join(sorted(list(fam_set)))
                                        if not te in all:
                                            all[te]=[]
                                        all[te].append(pos)
            if found_on_chr == False:
                log.logger.error('No MEs found in specified chromosomes. Exit.')
                exit(1)
        if len(all) == 0:
            log.logger.error('No MEs with "FILTER PASS" found. Exit.')
            exit(1)
        
        # header
        header=[]
        header.append('##fileformat=VCFv4.1\n')
        header.append('##fileDate=%s\n' % str(datetime.datetime.now()).split('.')[0])
        header.append('##source=Absent ME merge version "%s"\n' % args.version)
        header.append('##reference=%s\n' % args.fa)
        
        # count how many samples have a certain MEI
        load_header=False
        bed={}
        sample_id_to_dir={}
        for dir in args.dirs:
            f=args.dirs[dir][0]
            with open(f) as infile:
                for line in infile:
                    if line[0] == '#':
                        if load_header is False:
                            if '##contig=' in line or '##ALT=' in line or '##FILTER=' in line:
                                header.append(line)
                        if line[:6] == '#CHROM':
                            ls=line.split()
                            sample_id=ls[9]
                            sample_id_to_dir[sample_id]=dir
                            args.dirs[dir].append(sample_id)
                    else:
                        ls=line.strip().split('\t')
                        if ls[0] in args.chr:
                            if not 'Y' in ls[6]:
                                tmp=ls[7].split(';SVLEN=')[0].replace('MEI=', '')
                                fam_set=set()
                                clases=set()
                                for v in tmp.split(';'):
                                    clas,_=v.split(':')
                                    if clas in clas_to_fam:
                                        fam=clas_to_fam[clas]
                                        fam_set.add(fam)
                                        clases.add(clas)
                                if len(fam_set) >= 1:
                                    te='|'.join(sorted(list(fam_set)))
                                    clases=sorted(list(clases))
                                    leng=ls[7].split(';SVLEN=')[1].split(';')[0]
                                    start= int(ls[1]) - 1
                                    end= start + int(leng)
                                    infos=ls[7].split(';')
                                    if not te in bed:
                                        bed[te]=[]
                                    bed[te].append('\t'.join([ls[0], str(start), str(end), '%s,%s,%s,%s,%s' % (sample_id, ls[9], '.', leng, ls[6]), ','.join(clases)]))
                load_header=True
        args.sample_id_to_dir=sample_id_to_dir
        count={}
        poss={}
        for m in bed:
            if m in all:
                if not m in count:
                    count[m]={}
                    poss[m]={}
                all[m]=sorted(all[m], key=lambda x:(x[0], int(x[1]), int(x[2])))
                all[m]=[ '\t'.join(l) for l in all[m] ]
                all[m]=BedTool('\n'.join(all[m]), from_string=True).merge()
                for line in all[m]:
                    ls=str(line).split()
                    id='%s:%s-%s' % (ls[0], ls[1], ls[2])
                    count[m][id]=[]
                    poss[m][id]=[[], []]
                if not '|' in m:
                    log.logger.info('Merging %s....' % m)
                
                tmp=BedTool('\n'.join(bed[m]), from_string=True)
                tmp=tmp.intersect(all[m], wa=True, wb=True, nonamecheck=True)
                ss,es=[],[]
                for line in tmp:
                    ls=str(line).split()
                    id='%s:%s-%s' % (ls[5], ls[6], ls[7])
                    poss[m][id][0].append(ls[6])
                    poss[m][id][1].append(ls[7])
                    count[m][id].append([ls[3], ls[4]])
        # determine most common breakpoint
        for m in poss:
            for id in poss[m]:
                start=collections.Counter(poss[m][id][0]).most_common()[0][0]
                end=  collections.Counter(poss[m][id][1]).most_common()[0][0]
                poss[m][id]=[start, end]
        # retrieve fasta
        bed=[]
        for m in count:
            for id in count[m]:
                chr=id.split(':')[0]
                start= int(poss[m][id][0]) - 1  # add one base before bp
                bed.append('%s\t%d\t%s\t%s%s\n' % (chr, start, poss[m][id][1], m, id))
        bed=BedTool(''.join(bed), from_string=True)
        fa=bed.sequence(fi=args.fa, name=True)
        tmp=parse_fasta_for_merge_vcf(fa.seqfn)
        fa={}
        for h in tmp:
            fa[h.split('::')[0].replace('>', '')]=tmp[h].upper()
        # output
        n=0
        d={}
        sample_ids=set()
        for m in count:
            for id in count[m]:
                afs=[]    # allele freq
                lens=[]   # me length
                filts=[]  # me filter
                clases=set()
                for info,clas in count[m][id]:
                    infos=info.split(',')
                    afs.append(int(infos[1]))
                    if not infos[3] == 'NA':
                        lens.append(int(infos[3]))
                    filts.append(infos[4])
                    for v in clas.split(','):
                        clases.add(v)
                ac= sum(afs)
                pass_perc= filts.count('PASS') / len(filts)
                if pass_perc > params.pass_perc_threshold:
                    filt=collections.Counter(filts).most_common()[0][0]
                else:
                    filts_no_pass=[ v for v in filts if not v == 'PASS' ]
                    filt=collections.Counter(filts_no_pass).most_common()[0][0]
                melen= str(int(poss[m][id][1]) - int(poss[m][id][0]))
                pred_class='|'.join(sorted(list(clases)))
                d[id]={}
                d[id]['CHROM']=id.split(':')[0]
                d[id]['POS']=poss[m][id][0]  # 0-based to 1-based; add one base before bp
                d[id]['ID']='%s_abs_%d' % (args.cohort_name, n)
                d[id]['REF']=fa['%s%s' % (m, id)]
                d[id]['ALT']=fa['%s%s' % (m, id)][0]
                d[id]['QUAL']='.'
                d[id]['FILTER']=filt
                d[id]['INFO']='SVTYPE=%s;MEI=%s;0START=%s;0END=%s;SVLEN=%s;AC=%d' % (m, pred_class, poss[m][id][0], poss[m][id][1], melen, ac)
                d[id]['FORMAT']='GT'
                for info,_ in count[m][id]:
                    sample_id,c,_,_,_=info.split(',')
                    if c == '0':
                        gt='0/0'
                    elif c == '1':
                        gt='0/1'
                    elif c == '1.5':
                        gt='./1'
                    else:
                        gt='1/1'
                    d[id][sample_id]=gt
                    sample_ids.add(sample_id)
                n += 1
        sample_ids=sorted(list(sample_ids))
        for id in d:
            for sample_id in sample_ids:
                if not sample_id in d[id]:
                    d[id][sample_id]='0/0'
        header.append('##INFO=<ID=SVTYPE,Number=.,Type=String,Description="Mobile element class">\n')
        header.append('##INFO=<ID=MEI,Number=.,Type=String,Description="Mobile element subclass">\n')
        header.append('##INFO=<ID=0START,Number=1,Type=Integer,Description="Position of left breakpoint. 0-based numbering (eg BED style).">\n')
        header.append('##INFO=<ID=0END,Number=1,Type=Integer,Description="Position of right breakpoint. 0-based numbering (eg BED style).">\n')
        header.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
        header.append('##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count.">\n')
        header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header.append('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
        with gzip.open(filenames.merged_vcf_abs, 'wt') as outfile:
            outfile.write(''.join(header))
            outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % '\t'.join(sample_ids))
            col_names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + sample_ids
            for id in d:
                tmp=[]
                for c in col_names:
                    tmp.append(d[id][c])
                outfile.write('\t'.join(tmp) +'\n')
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def output_summary_for_fill_in(args, params, filenames):
    log.logger.debug('started')
    try:
        out=[]
        for dir in args.dirs:
            out.append('%s\t%s\n' % (args.dirs[dir][-1], args.dirs[dir][-2]))
        if args.merge_mei is True:
            with open(filenames.scaffold_info_i, 'w') as outfile:
                outfile.write(''.join(out))
        elif args.merge_absent_me is True:
            with open(filenames.scaffold_info_a, 'w') as outfile:
                outfile.write(''.join(out))
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
