#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,gzip,subprocess,statistics
import pysam
from pybedtools import BedTool
from scipy import stats
import numpy as np
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import log,traceback

matplotlib.rcParams['lines.linewidth']=0.5
matplotlib.rcParams['axes.linewidth']=0.5
matplotlib.rcParams['xtick.major.width']=0.5
matplotlib.rcParams['ytick.major.width']=0.5
matplotlib.rcParams['font.size']=5


cigar_op={'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}
cigar_ref_retain={'M', 'D', 'N', '=', 'X'}
cigar_read_retain={'M', 'I', '=', 'X'}
nums={'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}

def limit(args, params, filenames):
    log.logger.debug('started')
    try:
        if args.ins_bed is not None and args.abs_bed is not None:
            insbed=BedTool(args.ins_bed).slop(b=params.ins_slop_len, g=args.fai)
            absbed=BedTool(args.abs_bed).slop(b=params.abs_slop_len, g=args.fai)
            slopbed= insbed + absbed
            slopbed=slopbed.sort().merge()
        elif args.ins_bed is not None:
            slopbed=BedTool(args.ins_bed).slop(b=params.ins_slop_len, g=args.fai).sort().merge()
        else:
            slopbed=BedTool(args.abs_bed).slop(b=params.abs_slop_len, g=args.fai).sort().merge()
        if args.b is not None:
            cmd='samtools view -@ %d %s -bh -M -L %s -o %s' % (args.p, args.b, slopbed.fn, filenames.limited_b)
        else:
            cmd='samtools view -@ %d %s -T %s -Ch -M -L %s -o %s' % (args.p, args.c, args.fa, slopbed.fn, filenames.limited_c)
        log.logger.debug('samtools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during samtools running.')
            exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


# main
def evaluate_tsd_depth(args, params, filenames):
    log.logger.debug('started')
    try:
        # convert to depth
#        insbed=BedTool(args.ins_bed).slop(b=params.tsd_flank_len + 1, g=args.fai).sort().merge()
#        if args.b is not None:
#            cmd='samtools depth %s -a -b %s -o %s' % (filenames.limited_b, insbed.fn, filenames.depth_ins)
#        else:
#            cmd='samtools depth %s --reference %s -a -b %s -o %s' % (filenames.limited_c, args.fa, insbed.fn, filenames.depth_ins)
#        log.logger.debug('samtools depth command = `'+ cmd +'`')
#        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
#        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
#        if not out.returncode == 0:
#            log.logger.error('Error occurred during samtools depth running.')
#            exit(1)
        # load depth
        dep={}
        with open(filenames.depth_ins) as infile:
            for line in infile:
                ls=line.split()
                if not ls[0] in dep:
                    dep[ls[0]]={}
                dep[ls[0]][int(ls[1])]=int(ls[2])
        
        # calc depth at breakpoints
        def remove_1nt(lis):
            ave= sum(lis) / len(lis)
            diff_to_ave=[]
            for dep in lis:
                diff_to_ave.append(abs(ave - dep))
            max_diff=max(diff_to_ave)
            one_nt_removed=[]
            removed=False
            for dep,diff in zip(lis, diff_to_ave):
                if removed is False and diff == max_diff:
                    removed=True
                else:
                    one_nt_removed.append(dep)
            corrected_dep= sum(one_nt_removed) / len(one_nt_removed)
            return corrected_dep
            
        back_to_tsd_ratios={}
        only_tsd,only_del=[],[]
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.split()
                if ls[0] in dep:
                    left_pos=int(ls[4].split(',')[0].split('=')[1])
                    right_pos=int(ls[5].split(',')[0].split('=')[1])
                    if right_pos < left_pos:
                        # tsd
                        tsd_depth=[]
                        for pos in range(right_pos, left_pos):
                            tsd_depth.append(dep[ls[0]][pos])
                        if len(tsd_depth) >= params.min_tsd_len_to_remove_1nt:
                            tsd_depth=remove_1nt(tsd_depth)
                        else:
                            tsd_depth= sum(tsd_depth) / len(tsd_depth)
                        # left flank
                        left_depth=[]
                        for pos in range(right_pos - params.tsd_flank_len, right_pos):
                            left_depth.append(dep[ls[0]][pos])
                        left_depth=remove_1nt(left_depth)
                        # right flank
                        right_depth=[]
                        for pos in range(left_pos, left_pos + params.tsd_flank_len):
                            right_depth.append(dep[ls[0]][pos])
                        right_depth=remove_1nt(right_depth)
                        structure='TSD'
                    else:
                        if left_pos < right_pos:
                            # deleted region
                            tsd_depth=[]
                            for pos in range(left_pos, right_pos):
                                tsd_depth.append(dep[ls[0]][pos])
                            if len(tsd_depth) >= params.min_tsd_len_to_remove_1nt:
                                tsd_depth=remove_1nt(tsd_depth)
                            else:
                                tsd_depth= sum(tsd_depth) / len(tsd_depth)
                            structure='Del'
                        else:
                            # chimeric read num
                            left_num=int(ls[4].split(',')[1].split('=')[1])
                            right_num=int(ls[5].split(',')[1].split('=')[1])
                            tsd_depth= left_num + right_num
                            structure='no_TSD_no_Del'
                        # left flank
                        left_depth=[]
                        for pos in range(left_pos - params.tsd_flank_len, left_pos):
                            left_depth.append(dep[ls[0]][pos])
                        left_depth=remove_1nt(left_depth)
                        # right flank
                        right_depth=[]
                        for pos in range(right_pos, right_pos + params.tsd_flank_len):
                            right_depth.append(dep[ls[0]][pos])
                        right_depth=remove_1nt(right_depth)
                else:
                    left_depth,right_depth=0,0
                    left_num=int(ls[4].split(',')[1].split('=')[1])
                    right_num=int(ls[5].split(',')[1].split('=')[1])
                    tsd_depth= left_num + right_num
                    structure='no_background_read'
                if left_depth + right_depth > 0:
                    back_to_tsd_ratio= tsd_depth / ((left_depth + right_depth) / 2)
                else:
                    back_to_tsd_ratio=-1
                if back_to_tsd_ratio >= 0:
                    if structure == 'TSD' and 'confidence:high' in line:
                        only_tsd.append(back_to_tsd_ratio)
                    elif structure == 'Del' and 'confidence:high' in line:
                        only_del.append(back_to_tsd_ratio)
                back_to_tsd_ratios[ls[10]]=[back_to_tsd_ratio, structure]
        tsd_kernel=stats.gaussian_kde(only_tsd)
        tsd_x=np.linspace(0, 3, 600)
        tsd_y=tsd_kernel(tsd_x)
        del_kernel=stats.gaussian_kde(only_del)
        del_x=np.linspace(0, 3, 600)
        del_y=del_kernel(del_x)
        
        def find_threshold(xs, ys):
            highest=[-1, -1]
            peaks,bottoms=[],[]
            up=True
            prev_y=-1
            for x,y in zip(xs, ys):
                if y > highest[1]:
                    highest=[x, y]
                if y >= prev_y:
                    if up is False:
                        bottoms.append(x)
                    up=True
                else:
                    if up is True:
                        peaks.append(x)
                    up=False
                prev_y=y
            return peaks, bottoms, highest
        
        # find threshold
        xs=np.linspace(1, 2, 200)  # TSD
        ys=tsd_kernel(xs)
        peaks,bottoms,highest=find_threshold(xs, ys)
        if len(bottoms) >= 1:
            tsd_threshold=bottoms[-1]
            for peak in peaks:
                if peak < tsd_threshold:
                    mono_peak=peak
                else:
                    bi_peak=peak
                    break
            tsd_mono_high_conf_threshold= (mono_peak + tsd_threshold + tsd_threshold) / 3
            tsd_bi_high_conf_threshold=   (bi_peak   + tsd_threshold + tsd_threshold) / 3
        else:
            tsd_threshold= ((highest[0] - 1) * 1.5) + 1  # would be corner cases
            tsd_mono_high_conf_threshold= ((highest[0] - 1) * 1.333) + 1
            tsd_bi_high_conf_threshold=   ((highest[0] - 1) * 1.666) + 1
        log.logger.debug('tsd_kernel,peaks=%s,bottoms=%s,tsd_threshold=%f' % (peaks, bottoms, tsd_threshold))
        # find threshold, DEL
        xs=np.linspace(0, 1, 200)
        ys=del_kernel(xs)
        peaks,bottoms,highest=find_threshold(xs, ys)
        if len(bottoms) >= 1:
            del_threshold=bottoms[-1]
        else:
            del_threshold= 0.25  # would be corner cases
        del_mono_high_conf_threshold= del_threshold * params.del_mono_high_threshold_coeff  # 0.66
        del_bi_high_conf_threshold=   del_threshold * params.del_bi_high_threshold_coeff  # 1.33
        log.logger.debug('del_kernel,peaks=%s,bottoms=%s,del_threshold=%f' % (peaks, bottoms, del_threshold))
        # plot
        plt.figure(figsize=(3, 3))  # (x, y)
        gs=gridspec.GridSpec(2, 1)  # (y, x)
        gs.update(wspace=0.1)
        ax=plt.subplot(gs[0])  # TSD
        ax.fill([0]+ tsd_x.tolist() +[3], [0]+ tsd_y.tolist() +[0], c='dodgerblue', alpha=0.25, label='TSD, n=%d' % len(only_tsd))
        ax.axvline(x=tsd_threshold, linewidth=0.5, alpha=0.5, color='steelblue', linestyle='dashed', label='Threshold=%.2f' % tsd_threshold)
        ax.set_xlim(0, 3)
        ax.set_ylim(ymin=0)
        ax.set_ylabel('Density')
        ax.legend()
        ax=plt.subplot(gs[1])  # DEL
        ax.fill([0]+ del_x.tolist() +[3], [0]+ del_y.tolist() +[0], c='coral', alpha=0.25, label='Short_del, n=%d' % len(only_del))
        ax.axvline(x=del_threshold, linewidth=0.5, alpha=0.5, color='orangered', linestyle='dashed', label='Threshold=%.2f' % del_threshold)
        ax.set_xlim(0, 3)
        ax.set_ylim(ymin=0)
        ax.set_xlabel('Background depth to TSD/short-del depth ratio')
        ax.set_ylabel('Density')
        ax.legend()
        plt.suptitle('Gaussian kernel-density estimation')
        plt.savefig('./genotype_out/kde.pdf')
        
        # count allele count
        global cn_est_tsd_depth
        cn_est_tsd_depth={}
        for id in back_to_tsd_ratios:
            ratio,struc=back_to_tsd_ratios[id]
            allele=1
            if struc == 'TSD':
                if ratio < tsd_mono_high_conf_threshold:
                    allele='mono_high'
                elif tsd_mono_high_conf_threshold <= ratio < tsd_threshold:
                    allele='mono_low'
                elif tsd_threshold <= ratio < tsd_bi_high_conf_threshold:
                    allele='bi_low'
                elif tsd_bi_high_conf_threshold <= ratio < params.tsd_outlier:
                    allele='bi_high'
                else:
                    allele='outlier'
            elif struc == 'Del':
                if ratio < del_mono_high_conf_threshold:
                    allele='mono_high'
                elif del_mono_high_conf_threshold <= ratio < del_threshold:
                    allele='mono_low'
                elif del_threshold <= ratio < del_bi_high_conf_threshold:
                    allele='bi_low'
                elif del_bi_high_conf_threshold <= ratio < params.del_outlier:
                    allele='bi_high'
                else:
                    allele='outlier'
            else:
                allele='NA'
            cn_est_tsd_depth[id]=[allele, ratio, struc]
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def evaluate_spanning_read(args, params, filenames):
    log.logger.debug('started')
    try:
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
        
        def count_mismatch(cigar, mdz, ref_rel_start, ref_rel_end):
            # read dict, val = cigar
            read={}
            n=0
            tmp=''
            for c in cigar:
                if not c in cigar_op:
                    tmp += c
                elif c in cigar_read_retain:
                    for _ in range(int(tmp)):
                        read[n]=c
                        n += 1
                    tmp=''
                else:
                    tmp=''
            # add MD:Z info
            tmp=''
            mdz_l=[]
            for c in mdz:
                if c in nums:
                    tmp += c
                    mdz_del=False
                else:
                    if len(tmp) >= 1:
                        for _ in range(int(tmp)):
                            mdz_l.append('MI')
                    if c == '^':
                        mdz_del=True
                    if mdz_del is False:
                        mdz_l.append('X')
                    tmp=''
            for pos,c in zip(list(read.keys()), mdz_l):
                if c == 'X':
                    read[pos]='X'
            # ref_pos to read_pos dict
            ref_to_read={}
            ref_count,read_count=-1,-1
            tmp=''
            for c in cigar:
                if not c in cigar_op:
                    tmp += c
                else:
                    if c == 'M':
                        for _ in range(int(tmp)):
                            ref_count += 1
                            read_count += 1
                            ref_to_read[ref_count]=read_count
                    elif c == 'D':
                        for _ in range(int(tmp)):
                            ref_count += 1
                            ref_to_read[ref_count]=read_count
                    elif c == 'I':
                        for _ in range(int(tmp)):
                            read_count += 1
                    tmp=''
            ref_count += 1
            read_count += 1
            ref_to_read[ref_count]=read_count  # to use 0-based start, 1-based end data
            match=0
            for pos in range(ref_to_read[ref_rel_start], ref_to_read[ref_rel_end]):
                if read[pos] == 'M':
                    match += 1
            return match
        
        # load
        insbed=BedTool(args.ins_bed).slop(b=params.ins_slop_len, g=args.fai).sort().merge()
        if os.path.exists(filenames.limited_b +'.fai') is True:
            os.remove(filenames.limited_b +'.fai')
        if args.b is not None:
            pysam.index(filenames.limited_b)
            cmd='samtools view -@ %d %s -bh -M -L %s -o %s' % (args.p, filenames.limited_b, insbed.fn, filenames.tmp_bam)
        else:
            pysam.index(filenames.limited_c)
            cmd='samtools view -@ %d %s -T %s -bh -M -L %s -o %s' % (args.p, filenames.limited_c, args.fa, insbed.fn, filenames.tmp_bam)
        log.logger.debug('samtools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during samtools running.')
            exit(1)
        # load breakpoints
        bps={}
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.split()
                if not ls[0] in bps:
                    bps[ls[0]]=[]
                bps[ls[0]].append([int(ls[1]), int(ls[2]), ls[10]])
        # load reads
        infile=pysam.AlignmentFile(filenames.tmp_bam, 'rb')
        span_judge_true={}
        out_spanning=[]
        for line in infile:
            line=line.tostring()
            ls=line.strip().split('\t')
            if not 'S' in ls[5] and not 'H' in ls[5]:
                ref_len=calc_ref_len(ls[5])
                start= int(ls[3]) - 1  # 0-based
                end= start + ref_len
                for s,e,id in bps[ls[2]]:
                    if start < s < end:
                        if (start + params.min_len_overhang_for_spanning) < s < (end - params.min_len_overhang_for_spanning) and (start + params.min_len_overhang_for_spanning) < e < (end - params.min_len_overhang_for_spanning):
                            for info in ls[::-1]:
                                if 'MD:Z:' in info:
                                    mdz=info.replace('MD:Z:', '')
                                    break
                            left_match=count_mismatch(ls[5], mdz, 0, s - start)
                            right_match=count_mismatch(ls[5], mdz, e - start, end - start)
                            left_mismatch_ratio= (s - start - left_match) / (s - start)
                            right_mismatch_ratio= (end - e - right_match) / (end - e)
                            if left_match >= params.min_overhang_match_for_spanning and right_match >= params.min_overhang_match_for_spanning and left_mismatch_ratio <= params.max_overhang_mismatch_ratio_for_spanning and right_mismatch_ratio <= params.max_overhang_mismatch_ratio_for_spanning:
                                judge_span='True'
                                if not id in span_judge_true:
                                    span_judge_true[id]=0
                                span_judge_true[id] += 1
                            else:
                                judge_span='False'
                            b=bin(int(ls[1]))
                            if b[-1] == '1':  # paired-end
                                strand='/1' if b[-7] == '1' else '/2'
                            else:
                                strand=''
                            out_spanning.append('%s\t%s%s\t%s\t%d\t%d\t%d\t%d\n' % (id, ls[0], strand, judge_span, s - start, left_match, end - e, right_match))
        with gzip.open(filenames.out_spanning, 'wt') as outfile:
            outfile.write(''.join(out_spanning))
            outfile.flush()
            os.fdatasync(outfile.fileno())
        spanning_counts=[]
        for id in span_judge_true:
            spanning_counts.append(span_judge_true[id])
        spanning_count_median= statistics.median(spanning_counts)
        spanning_high_threshold= spanning_count_median * params.spanning_high_threshold_coeff
        spanning_zero_threshold= spanning_count_median * params.spanning_zero_threshold_coeff
        global cn_est_spanning
        cn_est_spanning={}
        spanning_outlier= args.cov * params.spanning_outlier_coeff
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.split()
                if ls[10] in span_judge_true:
                    if spanning_outlier <= span_judge_true[ls[10]]:
                        cn_est_spanning[ls[10]]=[span_judge_true[ls[10]], 'outlier']
                    elif spanning_high_threshold <= span_judge_true[ls[10]] < spanning_outlier:
                        cn_est_spanning[ls[10]]=[span_judge_true[ls[10]], 'bi_high']
                    elif spanning_zero_threshold <= span_judge_true[ls[10]] < spanning_high_threshold:
                        cn_est_spanning[ls[10]]=[span_judge_true[ls[10]], 'bi_low']
                    else:
                        cn_est_spanning[ls[10]]=[span_judge_true[ls[10]], 'mono_low']
                else:
                    cn_est_spanning[ls[10]]=[0, 'mono_high']
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
