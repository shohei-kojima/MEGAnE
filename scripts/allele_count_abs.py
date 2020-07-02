#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,gzip,subprocess
from pybedtools import BedTool
import pysam
import numpy as np
from scipy import stats
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
        absbed=BedTool(args.abs_bed).slop(b=params.abs_slop_len, g=args.fai).sort().merge()
        if os.path.exists(filenames.limited_b +'.bai') is True:
            os.remove(filenames.limited_b +'.bai')
        elif os.path.exists(filenames.limited_c +'.crai') is True:
            os.remove(filenames.limited_c +'.crai')
        if args.b is not None:
            pysam.index(filenames.limited_b)
            cmd='samtools view -@ %d %s -bh -M -L %s -o %s' % (args.p, filenames.limited_b, absbed.fn, filenames.tmp_bam)
        else:
            pysam.index(filenames.limited_c)
            cmd='samtools view -@ %d %s -T %s -bh -M -L %s -o %s' % (args.p, filenames.limited_c, args.fa, absbed.fn, filenames.tmp_bam)
        log.logger.debug('samtools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during samtools running.')
            exit(1)
        # load breakpoints
        id=0
        bps={}
        with open(args.abs_bed) as infile:
            for line in infile:
                ls=line.strip().split('\t')
                if not ls[0] in bps:
                    bps[ls[0]]=[]
                bps[ls[0]].append([int(ls[1]), int(ls[2]), 'ID=%d' % id])
                id += 1
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
                        if (start + params.min_len_overhang_for_spanning_abs) < s < (end - params.min_len_overhang_for_spanning_abs):
                            for info in ls[::-1]:
                                if 'MD:Z:' in info:
                                    mdz=info.replace('MD:Z:', '')
                                    break
                            left_match=count_mismatch(ls[5], mdz, 0, s - start)
                            right_match=count_mismatch(ls[5], mdz, s - start, end - start)
                            left_mismatch_ratio= (s - start - left_match) / (s - start)
                            right_mismatch_ratio= (end - s - right_match) / (end - s)
                            if left_match >= params.min_overhang_match_for_spanning_abs and right_match >= params.min_overhang_match_for_spanning_abs and left_mismatch_ratio <= params.max_overhang_mismatch_ratio_for_spanning_abs and right_mismatch_ratio <= params.max_overhang_mismatch_ratio_for_spanning_abs:
                                judge_span='True'
                                if not id in span_judge_true:
                                    span_judge_true[id]=[0, 0]
                                span_judge_true[id][0] += 1
                            else:
                                judge_span='False'
                            b=bin(int(ls[1]))
                            if b[-1] == '1':  # paired-end
                                strand='/1' if b[-7] == '1' else '/2'
                            else:
                                strand=''
                            out_spanning.append('%s\t%s%s\t%s\t%d\t%d\t%d\t%d\n' % (id, ls[0], strand, judge_span, s - start, left_match, end - s, right_match))
                    elif start < e < end:
                        if (start + params.min_len_overhang_for_spanning_abs) < e < (end - params.min_len_overhang_for_spanning_abs):
                            for info in ls[::-1]:
                                if 'MD:Z:' in info:
                                    mdz=info.replace('MD:Z:', '')
                                    break
                            left_match=count_mismatch(ls[5], mdz, 0, e - start)
                            right_match=count_mismatch(ls[5], mdz, e - start, end - start)
                            left_mismatch_ratio= (e - start - left_match) / (e - start)
                            right_mismatch_ratio= (end - e - right_match) / (end - e)
                            if left_match >= params.min_overhang_match_for_spanning_abs and right_match >= params.min_overhang_match_for_spanning_abs and left_mismatch_ratio <= params.max_overhang_mismatch_ratio_for_spanning_abs and right_mismatch_ratio <= params.max_overhang_mismatch_ratio_for_spanning_abs:
                                judge_span='True'
                                if not id in span_judge_true:
                                    span_judge_true[id]=[0, 0]
                                span_judge_true[id][1] += 1
                            else:
                                judge_span='False'
                            b=bin(int(ls[1]))
                            if b[-1] == '1':  # paired-end
                                strand='/1' if b[-7] == '1' else '/2'
                            else:
                                strand=''
                            out_spanning.append('%s\t%s%s\t%s\t%d\t%d\t%d\t%d\n' % (id, ls[0], strand, judge_span, e - start, left_match, end - e, right_match))
        with gzip.open(filenames.out_spanning, 'wt') as outfile:
            outfile.write(''.join(out_spanning))
            outfile.flush()
            os.fdatasync(outfile.fileno())
        spanning_counts=[]
        for id in span_judge_true:
            spanning_counts.append(sum(span_judge_true[id]))
        spanning_count_median= np.median(spanning_counts)
        spanning_high_threshold= spanning_count_median * params.spanning_high_threshold_coeff_abs
        spanning_zero_threshold= spanning_count_median * params.spanning_zero_threshold_coeff_abs
        global cn_est_spanning
        cn_est_spanning={}
        spanning_outlier= args.cov * params.spanning_outlier_coeff_abs
        log.logger.debug('spanning_zero_threshold=%f,spanning_high_threshold=%f,spanning_outlier=%f' % (spanning_zero_threshold, spanning_high_threshold, spanning_outlier))
        global spanning_thresholds
        spanning_thresholds=[spanning_zero_threshold, spanning_high_threshold, spanning_outlier]
        n=0
        with open(args.abs_bed) as infile:
            for line in infile:
                id='ID=%d' % n
                ls=line.strip().split('\t')
                if id in span_judge_true:
                    if spanning_outlier <= sum(span_judge_true[id]):
                        cn_est_spanning[id]=['outlier', span_judge_true[id][0], span_judge_true[id][1]]
                    elif spanning_high_threshold <= span_judge_true[id] < spanning_outlier:
                        cn_est_spanning[id]=['bi_high', span_judge_true[id][0], span_judge_true[id][1]]
                    elif spanning_zero_threshold <= span_judge_true[id] < spanning_high_threshold:
                        cn_est_spanning[id]=['bi_low', span_judge_true[id][0], span_judge_true[id][1]]
                    else:
                        cn_est_spanning[id]=['mono_low', span_judge_true[id][0], span_judge_true[id][1]]
                else:
                    cn_est_spanning[id]=['mono_high', 0, 0]
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def evaluate_bp_depth(args, params, filenames):
    log.logger.debug('started')
    try:
        # convert to depth
#        absbed=BedTool(args.abs_bed).slop(b=params.abs_flank_len + 1, g=args.fai).sort().merge()
#        if args.b is not None:
#            cmd='samtools depth %s -a -b %s -o %s' % (filenames.limited_b, absbed.fn, filenames.depth_abs)
#        else:
#            cmd='samtools depth %s --reference %s -a -b %s -o %s' % (filenames.limited_c, args.fa, absbed.fn, filenames.depth_abs)
#        log.logger.debug('samtools depth command = `'+ cmd +'`')
#        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
#        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
#        if not out.returncode == 0:
#            log.logger.error('Error occurred during samtools depth running.')
#            exit(1)
        
        # load depth
        dep={}
        with open(filenames.depth_abs) as infile:
            for line in infile:
                ls=line.split()
                if not ls[0] in dep:
                    dep[ls[0]]={}
                dep[ls[0]][int(ls[1]) - 1]=int(ls[2])
        
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
        
        global back_to_bp_ratios
        back_to_bp_ratios={}
        min_ratios=[]
        n=0
        with open(args.abs_bed) as infile:
            for line in infile:
                id='ID=%d' % n
                ls=line.strip().split('\t')
                if ls[0] in dep:
                    left_pos=int(ls[1])
                    right_pos=int(ls[2])
                    # left
                    left_flank=[]
                    for pos in range(left_pos - params.abs_flank_len, left_pos):
                        left_flank.append(dep[ls[0]][pos])
                    left_flank=remove_1nt(left_flank)
                    left_del=[]
                    for pos in range(left_pos, left_pos + params.abs_flank_len):
                        left_del.append(dep[ls[0]][pos])
                    left_del=remove_1nt(left_del)
                    # right
                    right_flank=[]
                    for pos in range(right_pos, right_pos + params.abs_flank_len):
                        right_flank.append(dep[ls[0]][pos])
                    right_flank=remove_1nt(right_flank)
                    right_del=[]
                    for pos in range(right_pos - params.abs_flank_len, right_pos):
                        right_del.append(dep[ls[0]][pos])
                    right_del=remove_1nt(right_del)
                    # judge
                    if left_flank > 0 and right_flank > 0:
                        left_ratio= left_del / left_flank
                        right_ratio= right_del / right_flank
                        min_ratios.append(min(left_ratio, right_ratio))
                    elif left_flank > 0:
                        left_ratio= left_del / left_flank
                        min_ratios.append(left_ratio)
                    elif right_flank > 0:
                        right_ratio= right_del / right_flank
                        min_ratios.append(right_ratio)
                    else:
                        left_ratio=0
                        right_ratio=0
                else:
                    left_ratio=0
                    right_ratio=0
                back_to_bp_ratios[id]=[left_ratio, right_ratio]
        abs_kernel=stats.gaussian_kde(min_ratios)
        abs_x=np.linspace(0, 2, 400)
        abs_y=abs_kernel(abs_x)
        
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
        xs=np.linspace(0, 1, 200)
        ys=abs_kernel(xs)
        peaks,bottoms,highest=find_threshold(xs, ys)
        plausibles=[]
        for peak in peaks:
            if 0.5 <= peak < 1:
                plausibles.append(peak)
        if len(plausibles) == 1:
            mono_peak=plausibles[0]
        else:
            mono_peak=params.mono_peak_notfound
        abs_mono_high_conf_threshold= mono_peak / 3
        log.logger.debug('abs_kernel,peaks=%s,bottoms=%s,highest=%s,abs_mono_high_conf_threshold=%f' % (peaks, bottoms, highest, abs_mono_high_conf_threshold))
        global abs_thresholds
        abs_thresholds=[abs_mono_high_conf_threshold, params.abs_depth_outlier]
        # plot
        plt.figure(figsize=(3, 2))  # (x, y)
        gs=gridspec.GridSpec(1, 1)  # (y, x)
        ax=plt.subplot(gs[0])  # DEL
        ax.fill([0]+ abs_x.tolist() +[2], [0]+ abs_y.tolist() +[0], c='coral', alpha=0.25, label='Absent, n=%d' % len(min_ratios))
        ax.axvline(x=abs_mono_high_conf_threshold, linewidth=0.5, alpha=0.5, color='orangered', linestyle='dashed', label='Threshold=%.2f' % abs_mono_high_conf_threshold)
        ax.set_xlim(0, 2)
        ax.set_ylim(ymin=0)
        ax.set_xlabel('Background depth to abs depth ratio')
        ax.set_ylabel('Density')
        ax.legend()
        plt.suptitle('Gaussian kernel-density estimation')
        plt.savefig('./genotype_out/kde_abs.pdf')
        plt.close()
        exit()
        # count allele count
        global cn_est_tsd_depth
        cn_est_tsd_depth={}
        for id in back_to_bp_ratios:
            ratio,struc=back_to_bp_ratios[id]
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








def find_abs(args, params, filenames):  # deprecated
    log.logger.debug('started')
    try:
        min_chimeric_num=params.min_chimeric_num
        mes,_=load_me_classification(filenames.reshaped_rep)
        
        # load chimeric reads
        d={}
        with gzip.open(filenames.abs_gz) as infile:
            for line in infile:
                line=line.decode()
                ls=line.split()
                if ls[1] in args.main_chrs_set:
                    R= int(ls[6].split('-')[1])
                    L= int(ls[7].split('-')[0])
                    id='\t'.join(ls[1:4])
                    if not id in d:
                        d[id]=[[],[]]
                    d[id][0].append(R - L)
                    d[id][1].append(ls[0])
        
        # identify hotspots
        bed_high_cov=[]
        for id in d:
            if len(d[id][0]) >= min_chimeric_num:
                m=Counter(d[id][0]).most_common(1)[0][0]
                bed_high_cov.append('%s\t%s\t%d\n' % (id, ';'.join(d[id][1]), m))
        del(d)
        bed_high_cov=BedTool(''.join(bed_high_cov), from_string=True)

        # identify unfixed TEs
        bed_te=BedTool(filenames.repout_bed)
        bed_te_intersect=bed_te.intersect(bed_high_cov, wa=True, wb=True)
        d={}
        for line in bed_te_intersect:
            ls=str(line).split()
            id='\t'.join(ls[5:8])
            if not id in d:
                d[id]=[]
            d[id].append([int(ls[1]), int(ls[2]), ls[3], ls[8], ls[4], ls[9]])

        outfile_abs=open(filenames.abs_res, 'w')
        outfile_transd=open(filenames.transd_res, 'w')
        abs_id=0
        for id in d:
            chr,start,end=id.split('\t')
            start,end=int(start),int(end)
            ss,es=[],[]
            for s,e,_,_,_,_ in d[id]:
                ss.append(abs(start - s))
                es.append(abs(end   - e))
            te_start=min(ss)
            te_end  =min(es)
            if (te_start <= params.breakpoint_annotation_gap) and (te_end <= params.breakpoint_annotation_gap):
                te=[]
                for s,e,n,_,_,_ in d[id]:
                    te.append('%s\t%d\t%d\n' % (chr, s, e))
                te=sorted(te, key=lambda x:(x[1], x[2]))
                te=BedTool(''.join(te), from_string=True).merge()
                abs_bed=BedTool('%s\t%d\t%d\n' % (chr, start, end), from_string=True)
                intersect=te.intersect(abs_bed)
                count=0
                for line in intersect:
                    ls=str(line).split()
                    count += int(ls[2]) - int(ls[1])
                te_ratio= count / (end - start)
                if te_ratio >= params.abs_len_to_te_ratio:
                    te_names=[]
                    non_ME_len=0
                    for s,e,n,r,_,t in d[id]:
                        if start <= s:
                            c=(e - s) if e <= end else (end - s)
                        else:
                            c=(e - start) if e <= end else (end - start)
                        if (c / (e-s)) >= params.len_te_for_abs_ratio:
                            te_names.append(n)
                        te_name=n.split(':')[0]
                        if not te_name in mes:
                            non_ME_len += c
                    if (len(te_names) >= 1) and ((non_ME_len / (end - start)) <= params.non_ME_len_ratio):
                        te_names=sorted(list(set(te_names)))
                        outfile_abs.write('%s\t%d\t%d\t%s\t%s\tTSD_len=%s\tid=%d\n' % (chr, start, end, ';'.join(te_names), r, t, abs_id))
                        abs_id += 1
            elif (te_start <= params.breakpoint_annotation_gap) or (te_end <= params.breakpoint_annotation_gap):
                transd=False
                if te_start <= params.breakpoint_annotation_gap:
                    for s,e,n,r,strand,t in d[id]:
                        if (abs(start - s) == te_start) and (strand == '+') and (e < end):
                            te_name=n.split(':')[0]
                            if te_name in mes:
                                te_clas=mes[te_name]
                                if te_clas in args.rep_with_pA:
                                    tail=BedTool('%s\t%d\t%d\n' % (chr, end - params.transduction_pA_len, end), from_string=True)
                                    tail=tail.sequence(fi=args.fa)
                                    fa=parse_fasta(tail.seqfn)
                                    seq=list(fa.values())[0]
                                    if seq.count('A') >= (params.transduction_pA_len * params.transduction_pA_ratio):
                                        transd=True
                elif te_end <= params.breakpoint_annotation_gap:
                    for s,e,n,r,strand,t in d[id]:
                        if (abs(end - e) == te_end) and (strand == '-') and (start < s):
                            te_name=n.split(':')[0]
                            if te_name in mes:
                                te_clas=mes[te_name]
                                if te_clas in args.rep_with_pA:
                                    tail=BedTool('%s\t%d\t%d\n' % (chr, start, start + params.transduction_pA_len), from_string=True)
                                    tail=tail.sequence(fi=args.fa)
                                    fa=parse_fasta(tail.seqfn)
                                    seq=list(fa.values())[0]
                                    if seq.count('T') >= (params.transduction_pA_len * params.transduction_pA_ratio):
                                        transd=True
                if transd is True:
                    outfile_transd.write('%s\t%d\t%d\t%s\t%s\tTSD_len=%s\tid=%d\n' % (chr, start, end, n, r, t, abs_id))
                    abs_id += 1
        outfile_abs.flush()
        outfile_transd.flush()
        os.fdatasync(outfile_abs.fileno())
        os.fdatasync(outfile_transd.fileno())
        outfile_abs.close()
        outfile_transd.close()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
