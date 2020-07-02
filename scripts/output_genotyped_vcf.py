#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,datetime
from utils import parse_fasta
from pybedtools import BedTool
import log,traceback


def output_ins_bed_vcf(args, params, filenames, data):
    log.logger.debug('started')
    try:
        # load orig bed
        orig={}
        bed=[]
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.strip().split('\t')
                orig[ls[-1]]=ls
                bed.append('%s\t%s\t%d\t%s\n' % (ls[0], ls[1], int(ls[1]) + 1, ls[-1]))
        # save results
        if not args.sample_name is None:
            input_sample=args.sample_name
        else:
            if args.b is not None:
                input_sample=os.path.basename(args.b)
            else:
                input_sample=os.path.basename(args.c)
        header=[]
        header.append('##fileformat=VCFv4.1\n')
        header.append('##fileDate=%s\n' % str(datetime.datetime.now()).split('.')[0])
        header.append('##source=MEI search version "%s"\n' % args.version)
        header.append('##reference=%s\n' % args.fa)
        with open(args.fai) as infile:
            for line in infile:
                ls=line.split()
                if ls[0] in args.main_chrs_set:
                    header.append('##contig=<ID=%s,length=%s>\n' % (ls[0], ls[1]))
        header.append('##ALT=<ID=INS:ME,Description="Insertion of mobile element">\n')
        header.append('##FILTER=<ID=LC,Description="Low confidence">\n')
        header.append('##FILTER=<ID=NU,Description="Not unique">\n')
        header.append('##FILTER=<ID=S,Description="Shorter than 50-bp">\n')
        header.append('##FILTER=<ID=G,Description="Outliers during genotyping">\n')
        header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        header.append('##INFO=<ID=MEPRED,Number=1,Type=String,Description="ME prediction status">\n')
        header.append('##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">\n')
        header.append('##INFO=<ID=MEI_rpos,Number=1,Type=Integer,Description="Reference poition of MEI right breakpoint">\n')
        header.append('##INFO=<ID=MEI,Number=.,Type=String,Description="Mobile element info">\n')
        header.append('##INFO=<ID=MEI_lbp,Number=.,Type=String,Description="Mobile element info of left breakpint">\n')
        header.append('##INFO=<ID=MEI_rbp,Number=.,Type=String,Description="Mobile element info of right breakpint">\n')
        header.append('##INFO=<ID=CN_conf,Number=.,Type=String,Description="Confidence of genotyping (e.g. CN). H=high_confidence,L=low_confidence">\n')
        header.append('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
        header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % input_sample)
        # retrieve ref seq
        bed=BedTool(''.join(bed), from_string=True)
        fa=bed.sequence(fi=args.fa, name=True)
        tmp=parse_fasta(fa.seqfn)
        fa={}
        for h in tmp:
            fa[h.split('::')[0]]=tmp[h].upper()
        # output
        out_vcf=[]
        out_bed=[]
        for id in orig:
            ls=orig[id]
            onebase_start= int(ls[1]) + 1  # 0-based to 1-based
            new_id=ls[-1].replace('ID=', '')
            new_id=new_id[:-1] if new_id[-1] == ';' else new_id
            seq=fa['>%s' % ls[-1]]
            filt=[]
            if ls[6] == 'confidence:low':
                filt.append('LC')
            if 'unique:no' in ls[7]:
                filt.append('NU')
            if '50bp_or_longer:no' in ls[7]:
                filt.append('S')
            if not data.merged_res[id][1] == 'PASS':
                filt.append(data.merged_res[id][1])
            if len(filt) == 0:
                filt='PASS'
            else:
                filt=';'.join(filt)
            left= int(ls[4].split(',')[0].replace('MEI_left:ref_pos=', ''))
            right= int(ls[5].split(',')[0].replace('MEI_right:ref_pos=', ''))
            homlen= (left - right) if left > right else 0
            meinfo=ls[8].split(',', 1)[1].replace('MEI_left_breakpoint', 'MEI_lbp').replace('MEI_right_breakpoint', 'MEI_rbp')
            mepred='PASS' if 'subfamily_pred:status=PASS' in ls[8] else 'FAILED'
            tmp=ls[:-1]
            tmp.append('%d;%s;%s' % tuple(data.merged_res[id]))
            tmp.append(id)
            out_bed.append('\t'.join(tmp) +'\n')
            info='SVTYPE=%s;MEPRED=%s;HOMLEN=%d;MEI_rpos=%s;%s;CN_conf=%s' % (ls[3], mepred, homlen, ls[2], meinfo, data.merged_res[id][2])
            vcfline='%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n' % (ls[0], onebase_start, id, seq, '<INS:ME>', '.', filt, info, 'CN', data.merged_res[id][0])
            out_vcf.append(vcfline)
        with open(filenames.ins_out_bed, 'w') as outfile:
            outfile.write(''.join(out_bed))
            outfile.flush()
            os.fdatasync(outfile.fileno())
        with open(filenames.ins_out_vcf, 'w') as outfile:
            outfile.write(''.join(header))
            outfile.write(''.join(out_vcf))
            outfile.flush()
            os.fdatasync(outfile.fileno())
        
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def output_abs_bed_vcf(args, params, filenames, data):
    log.logger.debug('started')
    try:
        # load orig bed
        orig={}
        bed=[]
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.strip().split('\t')
                orig[ls[-1]]=ls
                bed.append('%s\t%s\t%d\t%s\n' % (ls[0], ls[1], int(ls[1]) + 1, ls[-1]))
        # save results
        if not args.sample_name is None:
            input_sample=args.sample_name
        else:
            if args.b is not None:
                input_sample=os.path.basename(args.b)
            else:
                input_sample=os.path.basename(args.c)
        header=[]
        header.append('##fileformat=VCFv4.1\n')
        header.append('##fileDate=%s\n' % str(datetime.datetime.now()).split('.')[0])
        header.append('##source=MEI search version "%s"\n' % args.version)
        header.append('##reference=%s\n' % args.fa)
        with open(args.fai) as infile:
            for line in infile:
                ls=line.split()
                if ls[0] in args.main_chrs_set:
                    header.append('##contig=<ID=%s,length=%s>\n' % (ls[0], ls[1]))
        header.append('##ALT=<ID=INS:ME,Description="Insertion of mobile element">\n')
        header.append('##FILTER=<ID=LC,Description="Low confidence">\n')
        header.append('##FILTER=<ID=NU,Description="Not unique">\n')
        header.append('##FILTER=<ID=S,Description="Shorter than 50-bp">\n')
        header.append('##FILTER=<ID=G,Description="Outliers during genotyping">\n')
        header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        header.append('##INFO=<ID=MEPRED,Number=1,Type=String,Description="ME prediction status">\n')
        header.append('##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">\n')
        header.append('##INFO=<ID=MEI_rpos,Number=1,Type=Integer,Description="Reference poition of MEI right breakpoint">\n')
        header.append('##INFO=<ID=MEI,Number=.,Type=String,Description="Mobile element info">\n')
        header.append('##INFO=<ID=MEI_lbp,Number=.,Type=String,Description="Mobile element info of left breakpint">\n')
        header.append('##INFO=<ID=MEI_rbp,Number=.,Type=String,Description="Mobile element info of right breakpint">\n')
        header.append('##INFO=<ID=CN_conf,Number=.,Type=String,Description="Confidence of genotyping (e.g. CN). H=high_confidence,L=low_confidence">\n')
        header.append('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
        header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % input_sample)
        # retrieve ref seq
        bed=BedTool(''.join(bed), from_string=True)
        fa=bed.sequence(fi=args.fa, name=True)
        tmp=parse_fasta(fa.seqfn)
        fa={}
        for h in tmp:
            fa[h.split('::')[0]]=tmp[h].upper()
        # output
        out_vcf=[]
        out_bed=[]
        for id in orig:
            ls=orig[id]
            onebase_start= int(ls[1]) + 1  # 0-based to 1-based
            new_id=ls[-1].replace('ID=', '')
            new_id=new_id[:-1] if new_id[-1] == ';' else new_id
            seq=fa['>%s' % ls[-1]]
            filt=[]
            if ls[6] == 'confidence:low':
                filt.append('LC')
            if 'unique:no' in ls[7]:
                filt.append('NU')
            if '50bp_or_longer:no' in ls[7]:
                filt.append('S')
            if not data.merged_res[id][1] == 'PASS':
                filt.append(data.merged_res[id][1])
            if len(filt) == 0:
                filt='PASS'
            else:
                filt=';'.join(filt)
            left= int(ls[4].split(',')[0].replace('MEI_left:ref_pos=', ''))
            right= int(ls[5].split(',')[0].replace('MEI_right:ref_pos=', ''))
            homlen= (left - right) if left > right else 0
            meinfo=ls[8].split(',', 1)[1].replace('MEI_left_breakpoint', 'MEI_lbp').replace('MEI_right_breakpoint', 'MEI_rbp')
            mepred='PASS' if 'subfamily_pred:status=PASS' in ls[8] else 'FAILED'
            tmp=ls[:-1]
            tmp.append('%d;%s;%s' % tuple(data.merged_res[id]))
            tmp.append(id)
            out_bed.append('\t'.join(tmp) +'\n')
            info='SVTYPE=%s;MEPRED=%s;HOMLEN=%d;MEI_rpos=%s;%s;CN_conf=%s' % (ls[3], mepred, homlen, ls[2], meinfo, data.merged_res[id][2])
            vcfline='%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n' % (ls[0], onebase_start, id, seq, '<INS:ME>', '.', filt, info, 'CN', data.merged_res[id][0])
            out_vcf.append(vcfline)
        with open(filenames.ins_out_bed, 'w') as outfile:
            outfile.write(''.join(out_bed))
            outfile.flush()
            os.fdatasync(outfile.fileno())
        with open(filenames.ins_out_vcf, 'w') as outfile:
            outfile.write(''.join(header))
            outfile.write(''.join(out_vcf))
            outfile.flush()
            os.fdatasync(outfile.fileno())
        
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

