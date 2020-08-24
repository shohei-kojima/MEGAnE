#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,gzip,subprocess,math,datetime
import pysam
import pybedtools
from pybedtools import BedTool
from scipy import stats
from scipy.optimize import curve_fit
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
        pybedtools.set_tempdir(args.pybedtools_tmp)
        
        def generate_slopbed(args, params, filenames, ins_slop_len, abs_slop_len):
            do_ins=False if args.only_abs is True else True
            do_abs=False if args.only_ins is True else True
            slopbed=[]
            if do_ins is True:
                ins_bed_found=False
                tmp=[]
                for bp_final in [filenames.bp_final_g, filenames.bp_final_p, filenames.bp_final_f, filenames.bp_final_u]:
                    if os.path.exists(bp_final) is True:
                        ins_bed_found=True
                        with open(bp_final) as infile:
                            log.logger.debug('%s loading.' % bp_final)
                            for line in infile:
                                tmp.append(line)
                if ins_bed_found is False:
                    log.logger.error('Available ins_bed not found.')
                    exit(1)
                insbed=BedTool(''.join(tmp), from_string=True).slop(b=ins_slop_len, g=args.fai)
                for line in insbed:
                    ls=str(line).strip().split('\t')
                    slopbed.append('\t'.join(ls[:3]))
            if do_abs is True:
                abs_bed_found=False
                tmp=[]
                for abs_bed_f in [args.abs_bed, args.abs_3t_bed]:
                    if os.path.exists(abs_bed_f) is True:
                        abs_bed_found=True
                        with open(abs_bed_f) as infile:
                            log.logger.debug('%s loading.' % abs_bed_f)
                            for line in infile:
                                tmp.append(line)
                absbed=BedTool(''.join(tmp), from_string=True).slop(b=abs_slop_len, g=args.fai)
                if abs_bed_found is False:
                    log.logger.error('Available abs_bed not found.')
                    exit(1)
                for line in absbed:
                    ls=str(line).strip().split('\t')
                    slopbed.append('\t'.join(ls[:3]))
            slopbed=BedTool('\n'.join(slopbed), from_string=True)
            slopbed=slopbed.sort().merge()
            # slopbed; sort at the order of cram file
            tmp={}
            for chr in args.main_chrs:
                tmp[chr]=[]
            for line in slopbed:
                line=str(line)
                ls=line.split()
                tmp[ls[0]].append(line)
            slopbed=[]
            for chr in tmp:
                slopbed.extend(tmp[chr])
            slopbed=BedTool(''.join(slopbed), from_string=True)
            return slopbed
        
        slopbed=generate_slopbed(args, params, filenames, params.ins_slop_len, params.abs_slop_len)
        if args.unsorted is False:
            if args.b is not None:
                if os.path.exists(args.b + '.bai') is False:
                    pysam.index(args.b)
                cmd='samtools view -@ %d %s -bh -M -L %s -o %s' % (args.p, args.b, slopbed.fn, filenames.limited_b)
            else:
                if os.path.exists(args.c + '.crai') is False:
                    pysam.index(args.c)
                cmd='samtools view -@ %d %s -T %s -Ch -M -L %s -o %s' % (args.p, args.c, args.fa, slopbed.fn, filenames.limited_c)
        else:
            if args.b is not None:
                if os.path.exists(args.b + '.bai') is False:
                    pysam.index(args.b)
                cmd='samtools view -@ %d %s -bh -L %s -o %s' % (args.p, args.b, slopbed.fn, filenames.limited_tb)
            else:
                if os.path.exists(args.c + '.crai') is False:
                    pysam.index(args.c)
                cmd='samtools view -@ %d %s -T %s -Ch -L %s -o %s' % (args.p, args.c, args.fa, slopbed.fn, filenames.limited_tc)
        log.logger.debug('samtools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during samtools running.')
            exit(1)
        if args.unsorted is True:
            if args.b is not None:
                cmd='samtools sort -@ %d %s -m %s -T %s -O BAM -o %s' % (args.p, filenames.limited_tb, params.bam_sort_maxmem, args.outdir, filenames.limited_b)
            else:
                cmd='samtools sort -@ %d %s -m %s -T %s -O CRAM -o %s' % (args.p, filenames.limited_tc, params.bam_sort_maxmem, args.outdir, filenames.limited_c)
            log.logger.debug('samtools command = `'+ cmd +'`')
            out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
            log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
            if not out.returncode == 0:
                log.logger.error('Error occurred during samtools running.')
                exit(1)
            else:
                if args.b is not None:
                    os.remove(filenames.limited_tb)
                else:
                    os.remove(filenames.limited_tc)            
        
        # convert to depth
        slopbed=generate_slopbed(args, params, filenames, params.tsd_flank_len + 1, params.abs_flank_len + 1)
        if args.b is not None:
            cmd='samtools depth %s -a -b %s -o %s' % (filenames.limited_b, slopbed.fn, filenames.depth)
        else:
            cmd='samtools depth %s --reference %s -a -b %s -o %s' % (filenames.limited_c, args.fa, slopbed.fn, filenames.depth)
        log.logger.debug('samtools depth command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during samtools depth running.')
            exit(1)
        pybedtools.cleanup()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


# main
def evaluate_tsd_depth(args, params, filenames):
    log.logger.debug('started')
    try:
        # load depth
        def generate_d_for_depth(filenames, slop_len):
            d={}
            for bp_final in [filenames.bp_final_g, filenames.bp_final_p, filenames.bp_final_f, filenames.bp_final_u]:
                if os.path.exists(bp_final) is True:
                    with open(bp_final) as infile:
                        for line in infile:
                            ls=line.split()
                            if not ls[0] in d:
                                d[ls[0]]={}
                            start= int(ls[1]) - slop_len
                            end  = int(ls[2]) + slop_len
                            if start < 0:
                                start=0
                            for pos in range(start, end):
                                d[ls[0]][pos]=0
            return d
        
        dep=generate_d_for_depth(filenames, params.tsd_flank_len + 1)
        with open(filenames.depth) as infile:
            for line in infile:
                ls=line.split()
                if ls[0] in dep:
                    pos= int(ls[1]) - 1
                    if pos in dep[ls[0]]:
                        dep[ls[0]][pos]=int(ls[2])
        
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
            
        ############### from here; for plot_tsd_dep in misc.py; not related to the main analysis
#        left_tsds,middles,right_tsds=[],[],[]
#        with open(args.ins_bed) as infile:
#            for line in infile:
#                ls=line.strip().split('\t')
#                if ls[0] in dep:
#                    left_pos=int(ls[4].split(',')[0].split('=')[1])
#                    right_pos=int(ls[5].split(',')[0].split('=')[1])
#                    if right_pos < left_pos:
#                        # tsd
#                        tsd_depth=[]
#                        for pos in range(right_pos, left_pos):
#                            tsd_depth.append(dep[ls[0]][pos] / args.cov)
#                        if np.mean(tsd_depth) < 3:
#                            middles.append(np.mean(tsd_depth))
#                            # left flank
#                            left_depth=[]
#                            for pos in range(right_pos - 10, right_pos):
#                                left_depth.append(dep[ls[0]][pos] / args.cov)
#                            left_tsds.append(left_depth)
#                            # right flank
#                            right_depth=[]
#                            for pos in range(left_pos, left_pos + 10):
#                                right_depth.append(dep[ls[0]][pos] / args.cov)
#                            right_tsds.append(right_depth)
#        import misc
#        misc.plot_tsd_dep(left_tsds, middles, right_tsds)
#        exit()
        ############### until here
        # load fai
        fai={}
        with open(args.fai) as infile:
            for line in infile:
                ls=line.split()
                fai[ls[0]]=int(ls[1])
        
        # calc. ratio
        back_to_tsd_ratios={}
        only_tsd,only_del=[],[]
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.strip().split('\t')
                left_pos=int(ls[4].split(',')[0].split('=')[1])
                right_pos=int(ls[5].split(',')[0].split('=')[1])
                min_pos=min(left_pos, right_pos)
                max_pos=max(left_pos, right_pos)
                if ls[0] in dep and min_pos >= params.tsd_flank_len and (max_pos + params.tsd_flank_len) <= fai[ls[0]]:
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
        if len(only_del) >= 2:
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
            tsd_outlier_low_threshold= 1 + ((mono_peak - 1) * params.tsd_outlier_low_coeff)
        else:
            tsd_threshold= ((highest[0] - 1) * 1.5) + 1  # would be corner cases
            tsd_mono_high_conf_threshold= ((highest[0] - 1) * 1.333) + 1
            tsd_bi_high_conf_threshold=   ((highest[0] - 1) * 1.666) + 1
            tsd_outlier_low_threshold= 1 + (highest[0] * params.tsd_outlier_low_coeff)
        log.logger.debug('tsd_kernel,peaks=%s,bottoms=%s,tsd_threshold=%f,tsd_mono_high_conf_threshold=%f,tsd_bi_high_conf_threshold=%f' % (peaks, bottoms, tsd_threshold, tsd_mono_high_conf_threshold, tsd_bi_high_conf_threshold))
        # find threshold, DEL
        if len(only_del) >= 2:
            xs=np.linspace(0, 1, 200)
            ys=del_kernel(xs)
            peaks,bottoms,highest=find_threshold(xs, ys)
            if len(bottoms) >= 1:
                del_threshold=bottoms[-1]
            else:
                del_threshold= 0.25  # would be corner cases
            del_mono_high_conf_threshold= del_threshold * params.del_mono_high_threshold_coeff  # 0.66
            del_bi_high_conf_threshold=   del_threshold * params.del_bi_high_threshold_coeff  # 1.33
            log.logger.debug('del_kernel,peaks=%s,bottoms=%s,del_threshold=%f,del_mono_high_conf_threshold=%f,del_bi_high_conf_threshold=%f' % (peaks, bottoms, del_threshold, del_mono_high_conf_threshold, del_bi_high_conf_threshold))
        else:
            log.logger.debug('Not enough MEI with small del at breakpoint found (n=%d).' % len(only_del))
        global tsd_thresholds, del_thresholds
        tsd_thresholds=[tsd_mono_high_conf_threshold, tsd_threshold, tsd_bi_high_conf_threshold, params.tsd_outlier, tsd_outlier_low_threshold]
        if len(only_del) >= 2:
            del_thresholds=[del_mono_high_conf_threshold, del_threshold, del_bi_high_conf_threshold, params.del_outlier]
        else:
            del_thresholds=[None, None, None, None]
        # plot; this is for debug
#        plt.figure(figsize=(3, 3))  # (x, y)
#        gs=gridspec.GridSpec(2, 1)  # (y, x)
#        gs.update(wspace=0.1)
#        ax=plt.subplot(gs[0])  # TSD
#        ax.fill([0]+ tsd_x.tolist() +[3], [0]+ tsd_y.tolist() +[0], c='dodgerblue', alpha=0.25, label='TSD, n=%d' % len(only_tsd))
#        ax.axvline(x=tsd_threshold, linewidth=0.5, alpha=0.5, color='steelblue', linestyle='dashed', label='Threshold=%.2f' % tsd_threshold)
#        ax.set_xlim(0, 3)
#        ax.set_ylim(ymin=0)
#        ax.set_ylabel('Density')
#        ax.legend()
#        ax=plt.subplot(gs[1])  # DEL
#        ax.fill([0]+ del_x.tolist() +[3], [0]+ del_y.tolist() +[0], c='coral', alpha=0.25, label='Short_del, n=%d' % len(only_del))
#        ax.axvline(x=del_threshold, linewidth=0.5, alpha=0.5, color='orangered', linestyle='dashed', label='Threshold=%.2f' % del_threshold)
#        ax.set_xlim(0, 3)
#        ax.set_ylim(ymin=0)
#        ax.set_xlabel('Background depth to TSD/short-del depth ratio')
#        ax.set_ylabel('Density')
#        ax.legend()
#        plt.suptitle('Gaussian kernel-density estimation')
#        plt.savefig('./genotype_out/kde.pdf')
#        plt.close()
        
        # count allele count
        global cn_est_tsd_depth
        cn_est_tsd_depth={}
        for id in back_to_tsd_ratios:
            ratio,struc=back_to_tsd_ratios[id]
            if struc == 'TSD':
                if ratio < tsd_outlier_low_threshold:
                    allele='outlier'
                elif tsd_outlier_low_threshold <= ratio < tsd_mono_high_conf_threshold:
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
                if len(only_del) >= 2:
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
            else:
                allele='NA'
            cn_est_tsd_depth[id]=[allele, ratio, struc]
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def evaluate_spanning_read(args, params, filenames):
    log.logger.debug('started')
    try:
        pybedtools.set_tempdir(args.pybedtools_tmp)
        
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
        if os.path.exists(filenames.limited_b +'.bai') is True:
            os.remove(filenames.limited_b +'.bai')
        elif os.path.exists(filenames.limited_c +'.crai') is True:
            os.remove(filenames.limited_c +'.crai')
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
                ls=line.strip().split('\t')
                if not ls[0] in bps:
                    bps[ls[0]]=[]
                bps[ls[0]].append([int(ls[1]), int(ls[2]), ls[10]])
        # load reads
        infile=pysam.AlignmentFile(filenames.tmp_bam, 'rb')
        span_judge_true={}
        outfile=gzip.open(filenames.out_spanning, 'wt')
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
                            outfile.write('%s\t%s%s\t%s\t%d\t%d\t%d\t%d\n' % (id, ls[0], strand, judge_span, s - start, left_match, end - e, right_match))
        outfile.flush()
        os.fdatasync(outfile.fileno())
        outfile.close()
        spanning_counts=[]
        for id in span_judge_true:
            spanning_counts.append(span_judge_true[id])
        spanning_count_median= np.median(spanning_counts)
        spanning_high_threshold= spanning_count_median * params.spanning_high_threshold_coeff
        spanning_zero_threshold= spanning_count_median * params.spanning_zero_threshold_coeff
        global cn_est_spanning
        cn_est_spanning={}
        spanning_outlier= args.cov * params.spanning_outlier_coeff
        log.logger.debug('spanning_zero_threshold=%f,spanning_high_threshold=%f,spanning_outlier=%f' % (spanning_zero_threshold, spanning_high_threshold, spanning_outlier))
        global spanning_thresholds
        spanning_thresholds=[spanning_zero_threshold, spanning_high_threshold, spanning_outlier]
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.strip().split('\t')
                if ls[10] in span_judge_true:
                    if spanning_outlier <= span_judge_true[ls[10]]:
                        cn_est_spanning[ls[10]]=['outlier', span_judge_true[ls[10]]]
                    elif spanning_high_threshold <= span_judge_true[ls[10]] < spanning_outlier:
                        cn_est_spanning[ls[10]]=['mono_high', span_judge_true[ls[10]]]
                    elif spanning_zero_threshold <= span_judge_true[ls[10]] < spanning_high_threshold:
                        cn_est_spanning[ls[10]]=['mono_low', span_judge_true[ls[10]]]
                    else:
                        cn_est_spanning[ls[10]]=['bi_low', span_judge_true[ls[10]]]
                else:
                    cn_est_spanning[ls[10]]=['bi_high', 0]
        pybedtools.cleanup()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def evaluate_discordant(args, params, filenames):
    log.logger.debug('started')
    try:
        def gaussian_func_biallelics(coeff):
            def gaussian_func_biallelic(x, a, mu, sigma):
                return ((1-coeff)*a*np.exp(-(x-mu)**2/(2*sigma**2))) + (coeff*a*np.exp(-(x-(2*mu))**2/(4*sigma**2)))
            return gaussian_func_biallelic
        
        def gaussian_func(x, a, mu, sigma):
            return a*np.exp(-(x-mu)**2/(2*sigma**2))

        def fit_gaussian(list_support_read_count):
            x,y=[],[]
            support_read_bin= int(np.ceil(args.cov / 25))
            for i in range(0, max(list_support_read_count), support_read_bin):
                x.append(i + ((support_read_bin - 1) / 2))
                y.append(sum([ list_support_read_count.count(i + j) for j in range(support_read_bin) ]))
            init_param_a=max(y)
            x_few,y_few=[],[]
            for i,j in zip(x,y):
                if j < (init_param_a / 2):
                    x_few.append(i)
                    y_few.append(j)
                if j == init_param_a:
                    init_param_mu=i
                    break
            init_param=[init_param_a, init_param_mu, args.cov * params.fit_gaussian_init_sigma_coeff]
            log.logger.debug('init_param_a=%d,init_param_mu=%f,init_param_sigma=%f' %(init_param_a, init_param_mu, float(args.cov * params.fit_gaussian_init_sigma_coeff)))
            fits=[]
            for coeff in range(100):
                coeff= coeff / 100
                func=gaussian_func_biallelics(coeff)
                try:
                    popt,pcov=curve_fit(func, x, y, init_param)
                    residuals= y - func(x, *popt)  # all x
                    rss=np.sum(residuals**2)
                    tss=np.sum((y - np.mean(y))**2)
                    r_squared= 1 - (rss / tss)
                    residuals= y_few - func(x_few, *popt)  # few x
                    rss=np.sum(residuals**2)
                    tss=np.sum((y_few - np.mean(y_few))**2)
                    r_squared_few= 1 - (rss / tss)
                    fits.append([r_squared, r_squared_few, coeff, popt, pcov])
                except RuntimeError:
                    log.logger.debug('curve_fit,RuntimeError,coeff=%f' % coeff)
            fits=sorted(fits)
            if len(fits) >= 1:
                r_squared=fits[-1][0]
                r_squared_few=fits[-1][1]
                coeff=fits[-1][2]
                popt=fits[-1][3]
                pcov=fits[-1][4]
                log.logger.debug('r_squared=%f,biallelic_coeff=%f' %(r_squared, coeff))
                return x, y, popt, pcov, r_squared, coeff
            else:
                return False, False, False, False, False, False

        # main
        all_mei_count_range=[]
        for_gaussian_fitting=[]
        disc_read_d={}
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.strip().split('\t')
                discordant_l=int(ls[4].split(',')[1].replace('chimeric=', '')) + int(ls[4].split(',')[2].replace('hybrid=', ''))
                discordant_r=int(ls[5].split(',')[1].replace('chimeric=', '')) + int(ls[5].split(',')[2].replace('hybrid=', ''))
                count= discordant_l + discordant_r
                if ls[6] == 'confidence:high' and 'unique:yes' in ls[7]:
                    for_gaussian_fitting.append(count)
                all_mei_count_range.append(count)
                disc_read_d[ls[10]]=count
        global disc_thresholds, cn_est_disc
        if len(for_gaussian_fitting) < 10:
            log.logger.warning('Not enough data found. Will skip allele count estimation from discordant read number. Will use other evidences.')
            disc_thresholds=False
            cn_est_disc=False
        else:
            if args.b is not None:
                input_sample=os.path.basename(args.b)
            else:
                input_sample=os.path.basename(args.c)
            input_bed=os.path.basename(args.ins_bed)
            x,y,popt,pcov,r_squared,coeff=fit_gaussian(for_gaussian_fitting)  # gaussian fitting
            
            if x is False:
                log.logger.warning('Gaussian curve fitting failed. Will use other evidences.')
                disc_thresholds=False
                cn_est_disc=False
            else:
                log.logger.debug('popt=%s,pcov=%s' %(str(popt), str(pcov)))
                xd=np.arange(0, math.ceil(max(all_mei_count_range)) + 1)
                estimated_curve=gaussian_func_biallelics(coeff)(xd, popt[0], popt[1], popt[2])
                estimated_curve_single_allele=gaussian_func(xd, (1-coeff) * popt[0], popt[1], popt[2])
                estimated_curve_bi_allele=gaussian_func(xd, coeff * popt[0], 2 * popt[1], 1.414 * popt[2])
                # determine threshold
                # if threshold = zscore
#                mono_sd= popt[2] ** 2  # standard deviation; sigma ** 2
#                mono_zscore= (xd - popt[1]) / mono_sd  # zscore = (value - mu) / stdev
#                di_zscore= (xd - (popt[1] * 2)) / (mono_sd * 2)
#                for pos,mono,di in zip(xd,mono_zscore,di_zscore):
#                    if abs(mono) < abs(di):
#                        dosage[pos]=1
#                    else:
#                        dosage[pos]=2
                # if threshold = curve bottom
                large_mono=False
                for tmp_x,mono_y,bi_y in zip(xd, estimated_curve_single_allele, estimated_curve_bi_allele):
                    if mono_y >= bi_y:
                        large_mono=True
                    if mono_y < bi_y and large_mono is True:
                        disc_threshold=tmp_x
                        break
                if popt[1] < disc_threshold < (popt[1] * 2):
                   disc_mono_high_conf_threshold= (popt[1] + disc_threshold + disc_threshold) / 3
                   disc_di_high_conf_threshold= ((popt[1] * 2) + disc_threshold + disc_threshold) / 3
                else:  # would be corner cases
                    disc_threshold= popt[1] * 1.5
                    disc_mono_high_conf_threshold= popt[1] * 1.33
                    disc_di_high_conf_threshold= popt[1] * 1.66
                disc_outlier_threshold= popt[1] * params.discordant_outlier_coeff  # theoretically 3x depth
                log.logger.debug('disc_threshold=%f,disc_mono_high_conf_threshold=%f,disc_di_high_conf_threshold=%f,disc_outlier_threshold=%f' % (float(disc_threshold), disc_mono_high_conf_threshold, disc_di_high_conf_threshold, disc_outlier_threshold))
                disc_thresholds=[disc_mono_high_conf_threshold, disc_threshold, disc_di_high_conf_threshold, disc_outlier_threshold]
                # prep for plot
                mono_x,mono_y, di_x,di_y, outlier_x,outlier_y=[],[], [],[], [],[]
                for xval,yval in zip(x,y):
                    if xval < disc_threshold:
                        mono_x.append(xval)
                        mono_y.append(yval)
                    elif disc_threshold <= xval < disc_outlier_threshold:
                        di_x.append(xval)
                        di_y.append(yval)
                    else:
                        outlier_x.append(xval)
                        outlier_y.append(yval)
                # plot
                fig=plt.figure(figsize=(3,3))
                ax=fig.add_subplot(111)
                ax.scatter(mono_x, mono_y, s=5, c='dodgerblue', linewidths=0.5, alpha=0.5, label='CN=1')
                ax.scatter(di_x, di_y, s=5, c='coral', linewidths=0.5, alpha=0.5, label='CN=2')
                ax.scatter(outlier_x, outlier_y, s=5, c='grey', linewidths=0.5, alpha=0.5, label='CN=outlier')
                ax.plot(xd, estimated_curve_single_allele, color='grey', alpha=0.5)
                ax.plot(xd, estimated_curve_bi_allele, color='grey', alpha=0.5)
                ax.plot(xd, estimated_curve, label='Gaussian curve fitting', color='red', alpha=0.5)
                ax.set_xlim(0, popt[1] * 4)
                ax.set_xlabel('Number of discordant reads per breakpoint')
                ax.set_ylabel('Number of MEI')
                ax.legend()
                plt.suptitle('sample=%s;%s,\nn=%d, r_squared=%f,' % (input_sample, input_bed, len(for_gaussian_fitting), r_squared))  # popt[1] = mean, popt[2] = sigma
                plt.savefig(filenames.disc_read_pdf)
                plt.close()
                log.logger.debug('gaussian_fitting_n=%d,r_squared=%f' %(len(for_gaussian_fitting), r_squared))
                # count allele count
                cn_est_disc={}
                for id in disc_read_d:
                    if disc_read_d[id] < disc_mono_high_conf_threshold:
                        cn_est_disc[id]=['mono_high', disc_read_d[id]]
                    elif disc_mono_high_conf_threshold <= disc_read_d[id] < disc_threshold:
                        cn_est_disc[id]=['mono_low', disc_read_d[id]]
                    elif disc_threshold <= disc_read_d[id] < disc_di_high_conf_threshold:
                        cn_est_disc[id]=['bi_low', disc_read_d[id]]
                    elif disc_di_high_conf_threshold <= disc_read_d[id] < disc_outlier_threshold:
                        cn_est_disc[id]=['bi_high', disc_read_d[id]]
                    else:
                        cn_est_disc[id]=['outlier', disc_read_d[id]]
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


