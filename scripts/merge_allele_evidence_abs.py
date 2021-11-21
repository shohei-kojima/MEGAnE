#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,math
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('pdf')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import log,traceback

matplotlib.rcParams['lines.linewidth']=0.5
matplotlib.rcParams['axes.linewidth']=0.5
matplotlib.rcParams['xtick.major.width']=0.5
matplotlib.rcParams['ytick.major.width']=0.5
matplotlib.rcParams['font.size']=5


def plot_merged(args, params, filenames, data):
    log.logger.debug('started')
    try:
        # load genotyping results
        spanning_threshold_for_merge= math.ceil(args.cov * params.spanning_threshold_coeff_for_merge) * 2
        plt.figure(figsize=(4, 4))  # (x, y)
        gs=gridspec.GridSpec(5, 5, height_ratios=[0.05, 0.3, 0.1, 0.05, 0.3], width_ratios=[0.3, 0.05, 0.1, 0.3, 0.05])  # (y, x)
        gs.update(hspace=0.05, wspace=0.05)
        
        # skipping
        x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
        for id in data.cn_est_spanning:
            x.append(data.cn_est_spanning[id][2])
            y.append(data.cn_est_spanning[id][1])
            if data.merged_res[id][1] == 'PASS':
                if data.merged_res[id][0] == '1':
                    x_mono.append(data.cn_est_spanning[id][2])
                    y_mono.append(data.cn_est_spanning[id][1])
                else:
                    x_bi.append(data.cn_est_spanning[id][2])
                    y_bi.append(data.cn_est_spanning[id][1])
            else:
                x_failed.append(data.cn_est_spanning[id][2])
                y_failed.append(data.cn_est_spanning[id][1])
        
        ax=plt.subplot(gs[0,0])
        sns_x=[ i for i in x if i < data.spanning_thresholds[2] ]
        sns.violinplot(sns_x, orient='h', color='steelblue')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_xlim(-5, data.spanning_thresholds[2])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[1,1])
        sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
        sns.violinplot(sns_x, orient='v', color='steelblue')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_ylim(-5, data.spanning_thresholds[2])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[1,0])
        ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
        ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.1)
        ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.1)
        ax.axvline(x=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
        ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
        ax.set_xlabel('# spanning read, right breakpoint')
        ax.set_ylabel('# spanning read, left breakpoint')
        ax.set_xlim(-5, args.cov * params.spanning_outlier_coeff_abs)
        ax.set_ylim(-5, args.cov * params.spanning_outlier_coeff_abs)
        
        # relative depth to flanking
        x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
        for id in data.cn_est_depth:
            if data.cn_est_depth[id][1] >= 0 and data.cn_est_depth[id][2] >= 0:
                x.append(data.cn_est_depth[id][2])
                y.append(data.cn_est_depth[id][1])
                if data.merged_res[id][1] == 'PASS':
                    if data.merged_res[id][0] == '1':
                        x_mono.append(data.cn_est_depth[id][2])
                        y_mono.append(data.cn_est_depth[id][1])
                    else:
                        x_bi.append(data.cn_est_depth[id][2])
                        y_bi.append(data.cn_est_depth[id][1])
                else:
                    x_failed.append(data.cn_est_depth[id][2])
                    y_failed.append(data.cn_est_depth[id][1])
        
        ax=plt.subplot(gs[3,0])
        sns_x=[ i for i in x if i < 1.5 ]
        sns.violinplot(sns_x, orient='h', color='steelblue')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_xlim(-0.2, 1.5)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[4,1])
        sns_x=[ i for i in y if i < 1.5 ]
        sns.violinplot(sns_x, orient='v', color='steelblue')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_ylim(-0.2, 1.5)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[4,0])
        ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
        ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.1)
        ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.1)
        ax.axvline(x=data.abs_thresholds[0], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
        ax.axvline(x=data.abs_thresholds[1], linewidth=0.5, alpha=0.50, color='silver', linestyle='dashed')
        ax.axhline(y=data.abs_thresholds[0], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
        ax.axhline(y=data.abs_thresholds[1], linewidth=0.5, alpha=0.50, color='silver', linestyle='dashed')
        ax.set_xlabel('Relative depth to flanking, right breakpoint')
        ax.set_ylabel('Relative depth to flanking, left breakpoint')
        ax.set_xlim(-0.2, 1.5)
        ax.set_ylim(-0.2, 1.5)
        
        # spanning vs depth
        x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
        for id in data.cn_est_depth:
            if data.cn_est_depth[id][3] >= 0:
                x.append(data.cn_est_depth[id][3])
                y.append(data.cn_est_spanning[id][1] + data.cn_est_spanning[id][2])
                if data.merged_res[id][1] == 'PASS':
                    if data.merged_res[id][0] == '1':
                        x_mono.append(data.cn_est_depth[id][3])
                        y_mono.append(data.cn_est_spanning[id][1] + data.cn_est_spanning[id][2])
                    else:
                        x_bi.append(data.cn_est_depth[id][3])
                        y_bi.append(data.cn_est_spanning[id][1] + data.cn_est_spanning[id][2])
                else:
                    x_failed.append(data.cn_est_depth[id][3])
                    y_failed.append(data.cn_est_spanning[id][1] + data.cn_est_spanning[id][2])

        ax=plt.subplot(gs[0,3])
        sns_x=[ i for i in x if i < 1.5 ]
        sns.violinplot(sns_x, orient='h', color='steelblue')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_xlim(-0.2, 1.5)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[1,4])  # tsd, y=spanning
        sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
        sns.violinplot(sns_x, orient='v', color='steelblue')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_ylim(-5, data.spanning_thresholds[2])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[1,3])  # del, x=tsd, y=spanning
        ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
        ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.2)
        ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.2)
        ax.axvline(x=data.abs_thresholds[0], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
        ax.axvline(x=data.abs_thresholds[1], linewidth=0.5, alpha=0.50, color='silver', linestyle='dashed')
        ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
        ax.set_xlabel('Relative depth to flanking, min(left, right)')
        ax.set_ylabel('# spanning read, left + right')
        ax.set_xlim(-0.2, 1.5)
        ax.set_ylim(-5, args.cov * params.spanning_outlier_coeff_abs)
        
        plt.suptitle('Genotyping result for absent MEs')
        if args.no_pdf is False:
            plt.savefig(filenames.merged_pdf_abs)
        plt.close()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def merge(args, params, filenames, data):
    log.logger.debug('started')
    try:
        global merged_res
        merged_res={}
            
        def judge_outlier(id, data):
            tmp=[]
            if 'ID=3T' in id:
                tmp.append('3')
            if data.cn_est_depth[id][0] == 'outlier':
                tmp.append('D')
            if data.cn_est_spanning[id][0] == 'outlier':
                tmp.append('S')
            if len(tmp) == 0:
                tmp='PASS'
            else:
                tmp=';'.join(tmp)
            return tmp
            
        spanning_threshold_for_merge= math.ceil(args.cov * params.spanning_threshold_coeff_for_merge) * 2
        for id in data.cn_est_spanning:
            outlier_judge=judge_outlier(id, data)
            spanning_num= data.cn_est_spanning[id][1] + data.cn_est_spanning[id][2]
            if spanning_num >= spanning_threshold_for_merge:
                allele_count='1'
            elif data.cn_est_depth[id][3] < data.abs_thresholds[0]:
                allele_count='2'
            else:
                allele_count='1'
            merged_res[id]=[allele_count, outlier_judge]
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def merge_wo_discordant(args, params, filenames, data):
    log.logger.debug('started')
    try:
        global merged_res
        merged_res={}
            
        def judge_outlier(id, data):
            tmp=[]
            if 'ID=3T' in id:
                tmp.append('3')
            if data.cn_est_depth[id][0] == 'outlier':
                tmp.append('D')
            if data.cn_est_spanning[id][0] == 'outlier':
                tmp.append('S')
            if 'a%s' % id in data.disc_ids:
                if data.disc_ids['a%s' % id] == 'half':
                    tmp.append('H')
            else:
                tmp.append('R')
            if len(tmp) == 0:
                tmp='PASS'
            else:
                tmp=';'.join(tmp)
            return tmp
            
        outlier_threshold= args.cov * params.spanning_outlier_coeff_for_precall_abs
        spanning_threshold_for_merge= math.ceil(args.cov * params.spanning_threshold_coeff_for_merge) * 2
        for id in data.cn_est_spanning:
            outlier_judge=judge_outlier(id, data)
            if data.cn_est_spanning[id][1] >= outlier_threshold:
                allele_count='0'
            elif data.cn_est_spanning[id][2] >= outlier_threshold:
                allele_count='0'
            elif data.cn_est_depth[id][3] < data.abs_thresholds[0]:
                allele_count='2'
            else:
                allele_count='1'
            merged_res[id]=[allele_count, outlier_judge]
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
