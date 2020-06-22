#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,gzip,subprocess
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
            cmd='samtools view -@ %d %s -bh -M -L %s -o %s' % (args.p, args.b, slopbed.fn, filenames.limited_c)
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
            peaks,bottoms=[],[]
            up=True
            prev_y=-1
            for x,y in zip(xs, ys):
                if y >= prev_y:
                    if up is False:
                        bottoms.append(x)
                    up=True
                else:
                    if up is True:
                        peaks.append(x)
                    up=False
                prev_y=y
            return peaks, bottoms
        
        # find threshold
        xs=np.linspace(1, 2, 200)  # TSD
        ys=tsd_kernel(xs)
        peaks,bottoms=find_threshold(xs, ys)
        tsd_threshold=bottoms[-1]
        log.logger.debug('tsd_kernel,peaks=%s,bottoms=%s,tsd_threshold=%f' % (peaks, bottoms, tsd_threshold))
        # find threshold, DEL
        xs=np.linspace(0, 1, 200)
        ys=del_kernel(xs)
        peaks,bottoms=find_threshold(xs, ys)
        del_threshold=bottoms[-1]
        log.logger.debug('del_kernel,peaks=%s,bottoms=%s,del_threshold=%f' % (peaks, bottoms, del_threshold))
        # plot
        plt.figure(figsize=(3, 3))  # (x, y)
        gs=gridspec.GridSpec(2, 1)   # y, x
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
                if ratio >= tsd_threshold:
                    allele=2
            elif struc == 'Del':
                if ratio >= del_threshold:
                    allele=2
            else:
                allele='NA'
            cn_est_tsd_depth[id]=[allele, ratio, struc]
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def evaluate_spanning_read():
    log.logger.debug('started')
    try:
        pass
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
