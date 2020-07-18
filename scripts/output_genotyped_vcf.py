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
        header.append('##FILTER=<ID=R,Description="No discordant read stat available">\n')
        header.append('##FILTER=<ID=Y,Description="Variants on chrY. This is only available when female">\n')
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
        remove_Y=True if args.sex in params.female else False
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
            if remove_Y is True and ls[0] in params.chrY:
                filt.append('Y')
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
            tmp[6]='%s;%s' % (tmp[6], filt)
            tmp.append('%d;%s;%s' % tuple(data.merged_res[id]))
            tmp.append(';'.join([ str(i) for i in data.cn_est_tsd_depth[id] ]))
            tmp.append(';'.join([ str(i) for i in data.cn_est_spanning[id] ]))
            if not data.disc_thresholds is False:
                tmp.append(';'.join([ str(i) for i in data.cn_est_disc[id] ]))
            else:
                tmp.append('NA')
            tmp.append(id)
            out_bed.append('\t'.join(tmp) +'\n')
            info='SVTYPE=%s;MEPRED=%s;HOMLEN=%d;MEI_rpos=%s;%s;CN_conf=%s' % (ls[3], mepred, homlen, ls[2], meinfo, data.merged_res[id][2])
            vcfline='%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n' % (ls[0], onebase_start, new_id, seq, '<INS:ME>', '.', filt, info, 'CN', data.merged_res[id][0])
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
        n=0
        with open(args.abs_bed) as infile:
            for line in infile:
                id='ID=%d' % n
                n += 1
                ls=line.strip().split('\t')
                orig[id]=ls
                bed.append('%s\t%s\t%d\t%s\n' % (ls[0], ls[1], int(ls[1]) + 1, id))
        if not args.abs_3t_bed is None:
            n=0
            with open(args.abs_3t_bed) as infile:
                for line in infile:
                    id='ID=3T%d' % n
                    n += 1
                    ls=line.strip().split('\t')
                    orig[id]=ls
                    bed.append('%s\t%s\t%d\t%s\n' % (ls[0], ls[1], int(ls[1]) + 1, id))
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
        header.append('##source=Absent ME search version "%s"\n' % args.version)
        header.append('##reference=%s\n' % args.fa)
        with open(args.fai) as infile:
            for line in infile:
                ls=line.split()
                if ls[0] in args.main_chrs_set:
                    header.append('##contig=<ID=%s,length=%s>\n' % (ls[0], ls[1]))
        header.append('##ALT=<ID=DEL:ME,Description="Absence of reference mobile element">\n')
        header.append('##FILTER=<ID=3,Description="Potential 3\' transduction">\n')
        header.append('##FILTER=<ID=D,Description="Relative depth of breakpoint is outlier">\n')
        header.append('##FILTER=<ID=S,Description="Spanning read num is outlier">\n')
        header.append('##FILTER=<ID=Y,Description="Variants on chrY. This is only available when female">\n')
        header.append('##INFO=<ID=MEI,Number=.,Type=String,Description="Mobile element info">\n')
        header.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
        header.append('##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">\n')
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
        remove_Y=True if args.sex in params.female else False
        out_vcf=[]
        out_bed=[]
        for id in orig:
            ls=orig[id]
            onebase_start= int(ls[1]) + 1  # 0-based to 1-based
            seq=fa['>%s' % id]
            homlen=ls[5].replace('TSD_len=', '')
            svlen= int(ls[2]) - int(ls[1])
            if remove_Y is True and ls[0] in params.chrY:
                if data.merged_res[id][1] == 'PASS':
                    data.merged_res[id][1]='Y'
                else:
                    data.merged_res[id][1]=data.merged_res[id][1] +';Y'
            tmp=[]
            tmp.extend(ls)
            tmp.append('%s;%s' % tuple(data.merged_res[id]))
            tmp.append(';'.join([ str(i) for i in data.cn_est_depth[id] ]))
            tmp.append(';'.join([ str(i) for i in data.cn_est_spanning[id] ]))
            tmp.append(id)
            out_bed.append('\t'.join(tmp) +'\n')
            info='MEI=%s;SVLEN=%d;HOMLEN=%s' % (ls[3], svlen, homlen)
            vcfline='%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ls[0], onebase_start, id, seq, '<DEL:ME>', '.', data.merged_res[id][1], info, 'CN', data.merged_res[id][0])
            out_vcf.append(vcfline)
        with open(filenames.abs_out_bed, 'w') as outfile:
            outfile.write(''.join(out_bed))
            outfile.flush()
            os.fdatasync(outfile.fileno())
        with open(filenames.abs_out_vcf, 'w') as outfile:
            outfile.write(''.join(header))
            outfile.write(''.join(out_vcf))
            outfile.flush()
            os.fdatasync(outfile.fileno())
        
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

