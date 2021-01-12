#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,datetime,collections
import numpy as np
from pybedtools import BedTool
from utils import parse_fasta_for_merge_vcf
import log,traceback


def merge_vcf_ins(args, params, filenames):
    log.logger.debug('started')
    try:
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
        files=[]
        with open(args.f) as infile:
            for line in infile:
                line=line.strip()
                if os.path.exists(line) is False:
                    exit(1)
                else:
                    files.append(line)
        file_num=len(files)
                    
        # load high confidence and unique MEIs
        all={}
        for f in files:
            with open(f) as infile:
                for line in infile:
                    if not line[0] == '#':
                        ls=line.split()
                        if ls[6] == 'PASS':
                            te=ls[7].split(';')[0].replace('SVTYPE=', '')
                            pos=[ls[0], str(int(ls[1]) - 1), ls[7].split(';')[3].replace('MEI_rpos=', '')]
                            if not te in all:
                                all[te]=[]
                            all[te].append(pos)
        # header
        header=[]
        header.append('##fileformat=VCFv4.1\n')
        header.append('##fileDate=%s\n' % str(datetime.datetime.now()).split('.')[0])
        header.append('##source=MEI merge version "%s"\n' % args.version)
        header.append('##reference=%s\n' % args.fa)
        
        # count how many samples have a certain MEI
        load_header=False
        bed={}
        for f in files:
            sample_id=f
            with open(f) as infile:
                for line in infile:
                    if line[0] == '#':
                        if load_header is False:
                            if '##contig=' in line or '##ALT=' in line or '##FILTER=' in line:
                                header.append(line)
                        if line[:6] == '#CHROM':
                            ls=line.split()
                            sample_id=ls[9]
                    else:
                        ls=line.split()
                        if not 'Y' in ls[6]:
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
                                        pred_class.append('%s,%s' % (clas, clas))
                                        te_lens.append(te_len[clas] - int(bp))
                                        strands.append(strand)
                                elif infos[5] == 'MEI_rbp=pA':
                                    preds=infos[4].replace('MEI_lbp=', '')
                                    for pred in preds.split('|'):
                                        clas,bp,strand=pred.split(',')
                                        pred_class.append('%s,%s' % (clas, clas))
                                        te_lens.append(te_len[clas] - int(bp))
                                        strands.append(strand)
                            else:
                                mepred='FAILED'
                                left, right=infos[4].replace('MEI_lbp=', ''), infos[5].replace('MEI_rbp=', '')
                                for pred in left.split('|'):
                                    pred_class.append(pred.split(',')[0])
                                for pred in right.split('|'):
                                    pred_class.append(pred.split(',')[0])
                            if len(te_lens) > 0:
                                telen= str(int(np.round(np.mean(te_lens))))
                            else:
                                telen= 'NA'
                            if len(strands) > 0:
                                strand=collections.Counter(strands).most_common()[0][0]
                            else:
                                strand= 'NA'
                            if not te in bed:
                                bed[te]=[]
                            if len(pred_class) > 0:
                                bed[te].append('\t'.join([ls[0], start, end, '%s,%s,%s,%s,%s,%s' % (sample_id, ls[9], mepred, telen, ls[6], strand), ','.join(pred_class)]))
                load_header=True
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
                log.logger.info('Merging %s....' % m)
                
                tmp=BedTool('\n'.join(bed[m]), from_string=True)
                tmp=tmp.intersect(all[m], wa=True, wb=True)
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
                quals=[]  # me prediction quality
                lens=[]   # me length
                filts=[]  # me filter
                strands=[]  # me stands
                for info,_ in count[m][id]:
                    infos=info.split(',')
                    afs.append(int(infos[1]))
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
                d[id]['REF']=fa['%s%s' % (m, id)]
                d[id]['ALT']='<INS:ME>'
                d[id]['QUAL']='.'
                d[id]['FILTER']=filt
                d[id]['INFO']='SVTYPE=%s;MEI=%s;MEPRED=%s;0START=%s;1END=%s;SVLEN=%s;MESTRAND=%s;AC=%d' % (m, pred_class, mepred, poss[m][id][0], poss[m][id][1], melen, strand, ac)
                d[id]['FORMAT']='GT'
                for info,_ in count[m][id]:
                    sample_id,c,_,_,_,_=info.split(',')
                    if c == '0':
                        gt='0/0'
                    elif c == '1':
                        gt='0/1'
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
        header.append('##INFO=<ID=MEPRED,Number=.,Type=String,Description="MEI prediction status">\n')
        header.append('##INFO=<ID=0START,Number=1,Type=Integer,Description="Position of left breakpoint. 0-based numbering (eg BED style).">\n')
        header.append('##INFO=<ID=1END,Number=1,Type=Integer,Description="Position of right breakpoint. 1-based numbering (eg BED style).">\n')
        header.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
        header.append('##INFO=<ID=MESTRAND,Number=.,Type=String,Description="Strandness of MEI">\n')
        header.append('##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count.">\n')
        header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header.append('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
        with open(filenames.merged_vcf_ins, 'w') as outfile:
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


def merge_vcf_abs(args, params, filenames):
    log.logger.debug('started')
    try:
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
        files=[]
        with open(args.f) as infile:
            for line in infile:
                line=line.strip()
                if os.path.exists(line) is False:
                    exit(1)
                else:
                    files.append(line)
        file_num=len(files)
        
        # load high confidence and unique MEIs
        all={}
        for f in files:
            with open(f) as infile:
                for line in infile:
                    if not line[0] == '#':
                        ls=line.split()
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
        # header
        header=[]
        header.append('##fileformat=VCFv4.1\n')
        header.append('##fileDate=%s\n' % str(datetime.datetime.now()).split('.')[0])
        header.append('##source=Absent ME merge version "%s"\n' % args.version)
        header.append('##reference=%s\n' % args.fa)
        
        # count how many samples have a certain MEI
        load_header=False
        bed={}
        for f in files:
            sample_id=f
            with open(f) as infile:
                for line in infile:
                    if line[0] == '#':
                        if load_header is False:
                            if '##contig=' in line or '##ALT=' in line or '##FILTER=' in line:
                                header.append(line)
                        if line[:6] == '#CHROM':
                            ls=line.split()
                            sample_id=ls[9]
                    else:
                        ls=line.split()
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
                tmp=tmp.intersect(all[m], wa=True, wb=True)
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
                d[id]['INFO']='SVTYPE=%s;MEI=%s;0START=%s;1END=%s;SVLEN=%s;AC=%d' % (m, pred_class, poss[m][id][0], poss[m][id][1], melen, ac)
                d[id]['FORMAT']='GT'
                for info,_ in count[m][id]:
                    sample_id,c,_,_,_=info.split(',')
                    if c == '0':
                        gt='0/0'
                    elif c == '1':
                        gt='0/1'
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
        header.append('##INFO=<ID=1END,Number=1,Type=Integer,Description="Position of right breakpoint. 1-based numbering (eg BED style).">\n')
        header.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
        header.append('##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count.">\n')
        header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header.append('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
        with open(filenames.merged_vcf_abs, 'w') as outfile:
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

