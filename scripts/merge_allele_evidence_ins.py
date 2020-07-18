#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,math
from scipy import stats
import seaborn as sns
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import log,traceback

matplotlib.rcParams['lines.linewidth']=0.5
matplotlib.rcParams['axes.linewidth']=0.5
matplotlib.rcParams['xtick.major.width']=0.5
matplotlib.rcParams['ytick.major.width']=0.5
matplotlib.rcParams['font.size']=5


def plot_orig(args, params, filenames, data):
    log.logger.debug('started')
    try:
        if not data.disc_thresholds is False:
            plt.figure(figsize=(6.2, 4))  # (x, y)
            gs=gridspec.GridSpec(5, 8, height_ratios=[0.05, 0.3, 0.1, 0.05, 0.3], width_ratios=[0.3, 0.05, 0.1, 0.3, 0.05, 0.1, 0.3, 0.05])  # (y, x)
            gs.update(hspace=0.05, wspace=0.05)
            
            # 1
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
            
            ax=plt.subplot(gs[0,0])  # tsd, x=tsd
            sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,1])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,0])  # tsd, x=tsd, y=spanning
            ax.scatter(x, y, s=5, c='dodgerblue', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.tsd_thresholds[0], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[2], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[0], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[1], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, TSD')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])
            
            # 2
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_disc[id][1])
                    
            ax=plt.subplot(gs[3,0])  # tsd, x=tsd
            sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[4,1])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.disc_thresholds[3] ]
            sns.violinplot(sns_x, orient='v', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(0, data.disc_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[4,0])  # tsd, x=tsd, y=discordant reads
            ax.scatter(x, y, s=5, c='dodgerblue', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.tsd_thresholds[0], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[2], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, TSD')
            ax.set_ylabel('# discordant read')
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.set_ylim(0, data.disc_thresholds[3])
            
            # 3
            if not data.del_thresholds[0] is None:
                x,y=[],[]
                for id in data.cn_est_tsd_depth:
                    if data.cn_est_tsd_depth[id][2] == 'Del':
                        x.append(data.cn_est_tsd_depth[id][1])
                        y.append(data.cn_est_spanning[id][1])
                
                ax=plt.subplot(gs[0,3])  # del, x=tsd
                sns_x=[ i for i in x if i < data.del_thresholds[3] ]
                sns.violinplot(sns_x, orient='h', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                ax=plt.subplot(gs[1,4])  # tsd, y=spanning
                sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
                sns.violinplot(sns_x, orient='v', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_ylim(-5, data.spanning_thresholds[2])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])

                ax=plt.subplot(gs[1,3])  # del, x=tsd, y=spanning
                ax.scatter(x, y, s=5, c='coral', linewidths=0.5, alpha=0.1)
                ax.axvline(x=data.del_thresholds[0], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[2], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.axhline(y=data.spanning_thresholds[0], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axhline(y=data.spanning_thresholds[1], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.set_xlabel('Relative depth, Del')
                ax.set_ylabel('# spanning read')
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.set_ylim(-5, data.spanning_thresholds[2])
                
                # 4
                x,y=[],[]
                for id in data.cn_est_tsd_depth:
                    if data.cn_est_tsd_depth[id][2] == 'Del':
                        x.append(data.cn_est_tsd_depth[id][1])
                        y.append(data.cn_est_disc[id][1])
                
                ax=plt.subplot(gs[3,3])  # del, x=tsd
                sns_x=[ i for i in x if i < data.del_thresholds[3] ]
                sns.violinplot(sns_x, orient='h', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                ax=plt.subplot(gs[4,4])  # tsd, y=spanning
                sns_x=[ i for i in y if i < data.disc_thresholds[3] ]
                sns.violinplot(sns_x, orient='v', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_ylim(0, data.disc_thresholds[3])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])

                ax=plt.subplot(gs[4,3])  # del, x=tsd, y=discordant reads
                ax.scatter(x, y, s=5, c='coral', linewidths=0.5, alpha=0.1)
                ax.axvline(x=data.del_thresholds[0], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[2], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.axhline(y=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axhline(y=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axhline(y=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axhline(y=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.set_xlabel('Relative depth, Del')
                ax.set_ylabel('# discordant read')
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.set_ylim(0, data.disc_thresholds[3])
            
            # 5
            x,y=[],[]
            for id in data.cn_est_disc:
                x.append(data.cn_est_disc[id][1])
                y.append(data.cn_est_spanning[id][1])
            
            ax=plt.subplot(gs[0,6])  # del, x=tsd
            sns_x=[ i for i in x if i < data.disc_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='darkgreen')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.disc_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,7])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='darkgreen')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])

            ax=plt.subplot(gs[1,6])  # x=discordant, y=spanning reads
            ax.scatter(x, y, s=5, c='seagreen', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='forestgreen', linestyle='dashed')
            ax.axvline(x=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='forestgreen', linestyle='dashed')
            ax.axvline(x=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='forestgreen', linestyle='dashed')
            ax.axvline(x=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[0], linewidth=0.5, alpha=0.50, color='forestgreen', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[1], linewidth=0.5, alpha=0.25, color='forestgreen', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('# discordant read')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(0, data.disc_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])
        else:
            log.logger.warning('No discordant read stat available. Will use other evidences.')
            plt.figure(figsize=(4.1, 2))  # (x, y)
            gs=gridspec.GridSpec(2, 5, height_ratios=[0.05, 0.3], width_ratios=[0.3, 0.05, 0.1, 0.3, 0.05])  # (y, x)
            gs.update(hspace=0.05, wspace=0.05)
            
            # 1
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
            
            ax=plt.subplot(gs[0,0])  # tsd, x=tsd
            sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,1])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,0])  # tsd, x=tsd, y=spanning
            ax.scatter(x, y, s=5, c='dodgerblue', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.tsd_thresholds[0], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[2], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[0], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[1], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, TSD')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])

            # 2
            if not data.del_thresholds[0] is None:
                x,y=[],[]
                for id in data.cn_est_tsd_depth:
                    if data.cn_est_tsd_depth[id][2] == 'Del':
                        x.append(data.cn_est_tsd_depth[id][1])
                        y.append(data.cn_est_spanning[id][1])
                
                ax=plt.subplot(gs[0,3])  # del, x=tsd
                sns_x=[ i for i in x if i < data.del_thresholds[3] ]
                sns.violinplot(sns_x, orient='h', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                ax=plt.subplot(gs[1,4])  # tsd, y=spanning
                sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
                sns.violinplot(sns_x, orient='v', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_ylim(-5, data.spanning_thresholds[2])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])

                ax=plt.subplot(gs[1,3])  # del, x=tsd, y=spanning
                ax.scatter(x, y, s=5, c='coral', linewidths=0.5, alpha=0.1)
                ax.axvline(x=data.del_thresholds[0], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[2], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.axhline(y=data.spanning_thresholds[0], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axhline(y=data.spanning_thresholds[1], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.set_xlabel('Relative depth, Del')
                ax.set_ylabel('# spanning read')
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.set_ylim(-5, data.spanning_thresholds[2])
        plt.suptitle('This is a figure for debugging')
        plt.savefig(filenames.debug_pdf1)
        plt.close()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def plot_merged(args, params, filenames, data):
    log.logger.debug('started')
    try:
        # load mei search filter result
        global mei_filter
        mei_filter={}
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.strip().split()
                if ls[6] == 'confidence:high':
                    mei_filter[ls[10]]=True
                else:
                    mei_filter[ls[10]]=False
        # load genotyping results
        spanning_threshold_for_merge= 1 + math.ceil(args.cov * params.spanning_threshold_coeff_for_merge)
        if not data.disc_thresholds is False:
            plt.figure(figsize=(6.2, 4))  # (x, y)
            gs=gridspec.GridSpec(5, 8, height_ratios=[0.05, 0.3, 0.1, 0.05, 0.3], width_ratios=[0.3, 0.05, 0.1, 0.3, 0.05, 0.1, 0.3, 0.05])  # (y, x)
            gs.update(hspace=0.05, wspace=0.05)
            
            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
                    if mei_filter[id] is True:
                        if data.merged_res[id][0] == 1:
                            x_mono.append(data.cn_est_tsd_depth[id][1])
                            y_mono.append(data.cn_est_spanning[id][1])
                        else:
                            x_bi.append(data.cn_est_tsd_depth[id][1])
                            y_bi.append(data.cn_est_spanning[id][1])
                    else:
                        x_failed.append(data.cn_est_tsd_depth[id][1])
                        y_failed.append(data.cn_est_spanning[id][1])
            
            # 1
            ax=plt.subplot(gs[0,0])  # tsd, x=tsd
            sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,1])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,0])  # tsd, x=tsd, y=spanning
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.1)
            ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.tsd_thresholds[4], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, TSD')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])
            
            # 2
            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_disc[id][1])
                    if mei_filter[id] is True:
                        if data.merged_res[id][0] == 1:
                            x_mono.append(data.cn_est_tsd_depth[id][1])
                            y_mono.append(data.cn_est_disc[id][1])
                        else:
                            x_bi.append(data.cn_est_tsd_depth[id][1])
                            y_bi.append(data.cn_est_disc[id][1])
                    else:
                        x_failed.append(data.cn_est_tsd_depth[id][1])
                        y_failed.append(data.cn_est_disc[id][1])
                    
            ax=plt.subplot(gs[3,0])  # tsd, x=tsd
            sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[4,1])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.disc_thresholds[3] ]
            sns.violinplot(sns_x, orient='v', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(0, data.disc_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[4,0])  # tsd, x=tsd, y=discordant reads
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.1)
            ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.tsd_thresholds[4], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, TSD')
            ax.set_ylabel('# discordant read')
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.set_ylim(0, data.disc_thresholds[3])
            
            # 3
            if not data.del_thresholds[0] is None:
                x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
                for id in data.cn_est_tsd_depth:
                    if data.cn_est_tsd_depth[id][2] == 'Del':
                        x.append(data.cn_est_tsd_depth[id][1])
                        y.append(data.cn_est_spanning[id][1])
                        if mei_filter[id] is True:
                            if data.merged_res[id][0] == 1:
                                x_mono.append(data.cn_est_tsd_depth[id][1])
                                y_mono.append(data.cn_est_spanning[id][1])
                            else:
                                x_bi.append(data.cn_est_tsd_depth[id][1])
                                y_bi.append(data.cn_est_spanning[id][1])
                        else:
                            x_failed.append(data.cn_est_tsd_depth[id][1])
                            y_failed.append(data.cn_est_spanning[id][1])
                
                ax=plt.subplot(gs[0,3])  # del, x=tsd
                sns_x=[ i for i in x if i < data.del_thresholds[3] ]
                sns.violinplot(sns_x, orient='h', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                ax=plt.subplot(gs[1,4])  # tsd, y=spanning
                sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
                sns.violinplot(sns_x, orient='v', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_ylim(-5, data.spanning_thresholds[2])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])

                ax=plt.subplot(gs[1,3])  # del, x=tsd, y=spanning
                ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
                ax.scatter(x_mono, y_mono, s=5, c='gold', linewidths=0.5, alpha=0.2)
                ax.scatter(x_bi, y_bi, s=5, c='darkred', linewidths=0.5, alpha=0.2)
                ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.set_xlabel('Relative depth, Del')
                ax.set_ylabel('# spanning read')
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.set_ylim(-5, data.spanning_thresholds[2])
            
                # 4
                x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
                for id in data.cn_est_tsd_depth:
                    if data.cn_est_tsd_depth[id][2] == 'Del':
                        x.append(data.cn_est_tsd_depth[id][1])
                        y.append(data.cn_est_disc[id][1])
                        if mei_filter[id] is True:
                            if data.merged_res[id][0] == 1:
                                x_mono.append(data.cn_est_tsd_depth[id][1])
                                y_mono.append(data.cn_est_disc[id][1])
                            else:
                                x_bi.append(data.cn_est_tsd_depth[id][1])
                                y_bi.append(data.cn_est_disc[id][1])
                        else:
                            x_failed.append(data.cn_est_tsd_depth[id][1])
                            y_failed.append(data.cn_est_disc[id][1])

                
                ax=plt.subplot(gs[3,3])  # del, x=tsd
                sns_x=[ i for i in x if i < data.del_thresholds[3] ]
                sns.violinplot(sns_x, orient='h', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                ax=plt.subplot(gs[4,4])  # tsd, y=spanning
                sns_x=[ i for i in y if i < data.disc_thresholds[3] ]
                sns.violinplot(sns_x, orient='v', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_ylim(0, data.disc_thresholds[3])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])

                ax=plt.subplot(gs[4,3])  # del, x=tsd, y=discordant reads
                ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
                ax.scatter(x_mono, y_mono, s=5, c='gold', linewidths=0.5, alpha=0.2)
                ax.scatter(x_bi, y_bi, s=5, c='darkred', linewidths=0.5, alpha=0.2)
                ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.axhline(y=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axhline(y=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axhline(y=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
                ax.axhline(y=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.set_xlabel('Relative depth, Del')
                ax.set_ylabel('# discordant read')
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.set_ylim(0, data.disc_thresholds[3])
            
            # 5
            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_disc:
                x.append(data.cn_est_disc[id][1])
                y.append(data.cn_est_spanning[id][1])
                if mei_filter[id] is True:
                    if data.merged_res[id][0] == 1:
                        x_mono.append(data.cn_est_disc[id][1])
                        y_mono.append(data.cn_est_spanning[id][1])
                    else:
                        x_bi.append(data.cn_est_disc[id][1])
                        y_bi.append(data.cn_est_spanning[id][1])
                else:
                    x_failed.append(data.cn_est_disc[id][1])
                    y_failed.append(data.cn_est_spanning[id][1])
            
            ax=plt.subplot(gs[0,6])  # del, x=tsd
            sns_x=[ i for i in x if i < data.disc_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='darkgreen')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.disc_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,7])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='darkgreen')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])

            ax=plt.subplot(gs[1,6])  # x=discordant, y=spanning reads
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='lightgreen', linewidths=0.5, alpha=0.1)
            ax.scatter(x_bi, y_bi, s=5, c='darkolivegreen', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='forestgreen', linestyle='dashed')
            ax.axvline(x=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='forestgreen', linestyle='dashed')
            ax.axvline(x=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='forestgreen', linestyle='dashed')
            ax.axvline(x=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='forestgreen', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('# discordant read')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(0, data.disc_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])
        else:
            log.logger.warning('No discordant read stat available. Will use other evidences.')
            plt.figure(figsize=(4.1, 2))  # (x, y)
            gs=gridspec.GridSpec(2, 5, height_ratios=[0.05, 0.3], width_ratios=[0.3, 0.05, 0.1, 0.3, 0.05])  # (y, x)
            gs.update(hspace=0.05, wspace=0.05)
            
            # 1
            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
                    if mei_filter[id] is True:
                        if data.merged_res[id][0] == 1:
                            x_mono.append(data.cn_est_tsd_depth[id][1])
                            y_mono.append(data.cn_est_spanning[id][1])
                        else:
                            x_bi.append(data.cn_est_tsd_depth[id][1])
                            y_bi.append(data.cn_est_spanning[id][1])
                    else:
                        x_failed.append(data.cn_est_tsd_depth[id][1])
                        y_failed.append(data.cn_est_spanning[id][1])
            
            ax=plt.subplot(gs[0,0])  # tsd, x=tsd
            sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,1])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,0])  # tsd, x=tsd, y=spanning
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.1)
            ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, TSD')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])
            
            # 2
            if not data.del_thresholds[0] is None:
                x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
                for id in data.cn_est_tsd_depth:
                    if data.cn_est_tsd_depth[id][2] == 'Del':
                        x.append(data.cn_est_tsd_depth[id][1])
                        y.append(data.cn_est_spanning[id][1])
                        if mei_filter[id] is True:
                            if data.merged_res[id][0] == 1:
                                x_mono.append(data.cn_est_tsd_depth[id][1])
                                y_mono.append(data.cn_est_spanning[id][1])
                            else:
                                x_bi.append(data.cn_est_tsd_depth[id][1])
                                y_bi.append(data.cn_est_spanning[id][1])
                        else:
                            x_failed.append(data.cn_est_tsd_depth[id][1])
                            y_failed.append(data.cn_est_spanning[id][1])
                
                ax=plt.subplot(gs[0,3])  # del, x=tsd
                sns_x=[ i for i in x if i < data.del_thresholds[3] ]
                sns.violinplot(sns_x, orient='h', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                ax=plt.subplot(gs[1,4])  # tsd, y=spanning
                sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
                sns.violinplot(sns_x, orient='v', color='darkred')
                plt.setp(ax.collections, alpha=0.25)
                ax.set_ylim(-5, data.spanning_thresholds[2])
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])

                ax=plt.subplot(gs[1,3])  # del, x=tsd, y=spanning
                ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
                ax.scatter(x_mono, y_mono, s=5, c='gold', linewidths=0.5, alpha=0.2)
                ax.scatter(x_bi, y_bi, s=5, c='darkred', linewidths=0.5, alpha=0.2)
                ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
                ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
                ax.set_xlabel('Relative depth, Del')
                ax.set_ylabel('# spanning read')
                ax.set_xlim(-0.25, data.del_thresholds[3])
                ax.set_ylim(-5, data.spanning_thresholds[2])
        plt.suptitle('Genotyping result for insertions')
        plt.savefig(filenames.merged_pdf)
        plt.close()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def merge(args, params, filenames, data):
    log.logger.debug('started')
    try:
        global merged_res
        merged_res={}
        if not data.disc_thresholds is False:
            
            def judge_outlier(id, data):
                tmp=[]
                if data.cn_est_spanning[id][0] == 'outlier':
                    tmp.append('S')
                if data.cn_est_tsd_depth[id][0] == 'outlier':
                    tmp.append('D')
#                if data.cn_est_disc[id][0] == 'outlier':
#                    tmp.append('R')
                if len(tmp) == 0:
                    tmp='PASS'
                else:
                    tmp=''.join(tmp)
                return tmp
                
            spanning_threshold_for_merge= 1 + math.ceil(args.cov * params.spanning_threshold_coeff_for_merge)
            for id in data.cn_est_spanning:
                outlier_judge=judge_outlier(id, data)
                if data.cn_est_spanning[id][1] >= spanning_threshold_for_merge:
                    allele_count=1
                    if data.cn_est_tsd_depth[id][2] == 'TSD':
                        if data.cn_est_tsd_depth[id][1] < data.tsd_thresholds[1] and data.cn_est_disc[id][1] < data.disc_thresholds[2]:
                            status='H'
                        else:
                            status='L'
                    elif data.cn_est_tsd_depth[id][2] == 'Del':
                        if not data.del_thresholds[0] is None:
                            if data.cn_est_tsd_depth[id][1] >= data.del_thresholds[1] and data.cn_est_disc[id][1] >= data.disc_thresholds[0]:
                                status='H'
                            else:
                                status='L'
                        else:
                            status='L'
                    else:  # 'no_TSD_no_Del'
                        status='L'
                else:
                    if data.cn_est_tsd_depth[id][2] == 'TSD':
                        if data.cn_est_tsd_depth[id][1] >= data.tsd_thresholds[1] and data.cn_est_disc[id][1] >= data.disc_thresholds[0]:
                            allele_count=2
                            status='H'
                        else:
                            allele_count=1
                            status='L'
                    elif data.cn_est_tsd_depth[id][2] == 'Del':
                        if not data.del_thresholds[0] is None:
                            if data.cn_est_tsd_depth[id][1] < data.del_thresholds[1] and data.cn_est_disc[id][1] >= data.disc_thresholds[0]:
                                allele_count=2
                                status='H'
                            else:
                                allele_count=1
                                status='L'
                        else:
                            status='L'
                            if data.cn_est_disc[id][1] >= data.disc_thresholds[0]:
                                allele_count=2
                            else:
                                allele_count=1
                    else:  # 'no_TSD_no_Del' and 'no_background_read'
                        if data.cn_est_disc[id][1] >= data.disc_thresholds[2]:
                            allele_count=2
                            status='L'
                        else:
                            allele_count=1
                            status='L'
                if not outlier_judge == 'PASS':
                    status='L'
                merged_res[id]=[allele_count, outlier_judge, status]
        else:
            log.logger.warning('No discordant read stat available. Will use other evidences.')
            def judge_outlier(id, data):
                tmp=[]
                if data.cn_est_spanning[id][0] == 'outlier':
                    tmp.append('S')
                if data.cn_est_tsd_depth[id][0] == 'outlier':
                    tmp.append('D')
                tmp.append('R')
                tmp=''.join(tmp)
                return tmp
            
            status='L'
            spanning_threshold_for_merge= 1 + math.ceil(args.cov * params.spanning_threshold_coeff_for_merge)
            for id in data.cn_est_spanning:
                outlier_judge=judge_outlier(id, data)
                if data.cn_est_spanning[id][1] >= spanning_threshold_for_merge:
                    allele_count=1
                else:
                    if data.cn_est_tsd_depth[id][2] == 'TSD':
                        if data.cn_est_tsd_depth[id][1] >= data.tsd_thresholds[1]:
                            allele_count=2
                        else:
                            allele_count=1
                    elif data.cn_est_tsd_depth[id][2] == 'Del':
                        if not data.del_thresholds[0] is None:
                            if data.cn_est_tsd_depth[id][1] < data.del_thresholds[1]:
                                allele_count=2
                            else:
                                allele_count=1
                        else:
                            allele_count=1
                    else:  # 'no_TSD_no_Del' and 'no_background_read'
                        allele_count=2
                merged_res[id]=[allele_count, outlier_judge, status]
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def plot_gt(args, params, filenames, data):
    log.logger.debug('started')
    try:
        # load mei search filter result
        global mei_filter
        mei_filter={}
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.strip().split()
                if ls[6] == 'confidence:high':
                    mei_filter[ls[10]]=True
                else:
                    mei_filter[ls[10]]=False
        # load GT
        gt=set()
        with open('./genotype_out/gt_intersect.bed') as infile:
            for line in infile:
                ls=line.split()
                gt.add(ls[10])
        print(len(gt))
        # load genotyping results
        spanning_threshold_for_merge= 1 + math.ceil(args.cov * params.spanning_threshold_coeff_for_merge)
        plt.figure(figsize=(6.2, 4))  # (x, y)
        gs=gridspec.GridSpec(5, 8, height_ratios=[0.05, 0.3, 0.1, 0.05, 0.3], width_ratios=[0.3, 0.05, 0.1, 0.3, 0.05, 0.1, 0.3, 0.05])  # (y, x)
        gs.update(hspace=0.05, wspace=0.05)
        
        x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
        for id in data.cn_est_tsd_depth:
            if data.cn_est_tsd_depth[id][2] == 'TSD':
                x.append(data.cn_est_tsd_depth[id][1])
                y.append(data.cn_est_spanning[id][1])
                if id in gt:
                    if data.merged_res[id][0] == 1:
                        x_mono.append(data.cn_est_tsd_depth[id][1])
                        y_mono.append(data.cn_est_spanning[id][1])
                    else:
                        x_bi.append(data.cn_est_tsd_depth[id][1])
                        y_bi.append(data.cn_est_spanning[id][1])
                else:
                    x_failed.append(data.cn_est_tsd_depth[id][1])
                    y_failed.append(data.cn_est_spanning[id][1])
        
        ax=plt.subplot(gs[0,0])  # tsd, x=tsd
        sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
        sns.violinplot(sns_x, orient='h', color='steelblue')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_xlim(0, data.tsd_thresholds[3])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[1,1])  # tsd, y=spanning
        sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
        sns.violinplot(sns_x, orient='v', color='steelblue')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_ylim(-5, data.spanning_thresholds[2])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[1,0])  # tsd, x=tsd, y=spanning
        ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
        ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.1)
        ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.1)
        ax.axvline(x=data.tsd_thresholds[4], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
        ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
        ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
        ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
        ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
        ax.set_xlabel('Relative depth, TSD')
        ax.set_ylabel('# spanning read')
        ax.set_xlim(0, data.tsd_thresholds[3])
        ax.set_ylim(-5, data.spanning_thresholds[2])
        
        x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
        for id in data.cn_est_tsd_depth:
            if data.cn_est_tsd_depth[id][2] == 'TSD':
                x.append(data.cn_est_tsd_depth[id][1])
                y.append(data.cn_est_disc[id][1])
                if id in gt:
                    if data.merged_res[id][0] == 1:
                        x_mono.append(data.cn_est_tsd_depth[id][1])
                        y_mono.append(data.cn_est_disc[id][1])
                    else:
                        x_bi.append(data.cn_est_tsd_depth[id][1])
                        y_bi.append(data.cn_est_disc[id][1])
                else:
                    x_failed.append(data.cn_est_tsd_depth[id][1])
                    y_failed.append(data.cn_est_disc[id][1])
                
        ax=plt.subplot(gs[3,0])  # tsd, x=tsd
        sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
        sns.violinplot(sns_x, orient='h', color='steelblue')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_xlim(0, data.tsd_thresholds[3])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[4,1])  # tsd, y=spanning
        sns_x=[ i for i in y if i < data.disc_thresholds[3] ]
        sns.violinplot(sns_x, orient='v', color='steelblue')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_ylim(0, data.disc_thresholds[3])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[4,0])  # tsd, x=tsd, y=discordant reads
        ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
        ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.1)
        ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.1)
        ax.axvline(x=data.tsd_thresholds[4], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
        ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
        ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
        ax.axhline(y=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
        ax.axhline(y=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
        ax.axhline(y=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
        ax.axhline(y=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
        ax.set_xlabel('Relative depth, TSD')
        ax.set_ylabel('# discordant read')
        ax.set_xlim(0, data.tsd_thresholds[3])
        ax.set_ylim(0, data.disc_thresholds[3])

        x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
        for id in data.cn_est_tsd_depth:
            if data.cn_est_tsd_depth[id][2] == 'Del':
                x.append(data.cn_est_tsd_depth[id][1])
                y.append(data.cn_est_spanning[id][1])
                if id in gt:
                    if data.merged_res[id][0] == 1:
                        x_mono.append(data.cn_est_tsd_depth[id][1])
                        y_mono.append(data.cn_est_spanning[id][1])
                    else:
                        x_bi.append(data.cn_est_tsd_depth[id][1])
                        y_bi.append(data.cn_est_spanning[id][1])
                else:
                    x_failed.append(data.cn_est_tsd_depth[id][1])
                    y_failed.append(data.cn_est_spanning[id][1])
        
        ax=plt.subplot(gs[0,3])  # del, x=tsd
        sns_x=[ i for i in x if i < data.del_thresholds[3] ]
        sns.violinplot(sns_x, orient='h', color='darkred')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_xlim(-0.25, data.del_thresholds[3])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[1,4])  # tsd, y=spanning
        sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
        sns.violinplot(sns_x, orient='v', color='darkred')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_ylim(-5, data.spanning_thresholds[2])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

        ax=plt.subplot(gs[1,3])  # del, x=tsd, y=spanning
        ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
        ax.scatter(x_mono, y_mono, s=5, c='gold', linewidths=0.5, alpha=0.2)
        ax.scatter(x_bi, y_bi, s=5, c='darkred', linewidths=0.5, alpha=0.2)
        ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
        ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
        ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
        ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
        ax.set_xlabel('Relative depth, Del')
        ax.set_ylabel('# spanning read')
        ax.set_xlim(-0.25, data.del_thresholds[3])
        ax.set_ylim(-5, data.spanning_thresholds[2])
        
        x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
        for id in data.cn_est_tsd_depth:
            if data.cn_est_tsd_depth[id][2] == 'Del':
                x.append(data.cn_est_tsd_depth[id][1])
                y.append(data.cn_est_disc[id][1])
                if id in gt:
                    if data.merged_res[id][0] == 1:
                        x_mono.append(data.cn_est_tsd_depth[id][1])
                        y_mono.append(data.cn_est_disc[id][1])
                    else:
                        x_bi.append(data.cn_est_tsd_depth[id][1])
                        y_bi.append(data.cn_est_disc[id][1])
                else:
                    x_failed.append(data.cn_est_tsd_depth[id][1])
                    y_failed.append(data.cn_est_disc[id][1])

        
        ax=plt.subplot(gs[3,3])  # del, x=tsd
        sns_x=[ i for i in x if i < data.del_thresholds[3] ]
        sns.violinplot(sns_x, orient='h', color='darkred')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_xlim(-0.25, data.del_thresholds[3])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[4,4])  # tsd, y=spanning
        sns_x=[ i for i in y if i < data.disc_thresholds[3] ]
        sns.violinplot(sns_x, orient='v', color='darkred')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_ylim(0, data.disc_thresholds[3])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

        ax=plt.subplot(gs[4,3])  # del, x=tsd, y=discordant reads
        ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
        ax.scatter(x_mono, y_mono, s=5, c='gold', linewidths=0.5, alpha=0.2)
        ax.scatter(x_bi, y_bi, s=5, c='darkred', linewidths=0.5, alpha=0.2)
        ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
        ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
        ax.axhline(y=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
        ax.axhline(y=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
        ax.axhline(y=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
        ax.axhline(y=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
        ax.set_xlabel('Relative depth, Del')
        ax.set_ylabel('# discordant read')
        ax.set_xlim(-0.25, data.del_thresholds[3])
        ax.set_ylim(0, data.disc_thresholds[3])

        x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
        for id in data.cn_est_disc:
            x.append(data.cn_est_disc[id][1])
            y.append(data.cn_est_spanning[id][1])
            if id in gt:
                if data.merged_res[id][0] == 1:
                    x_mono.append(data.cn_est_disc[id][1])
                    y_mono.append(data.cn_est_spanning[id][1])
                else:
                    x_bi.append(data.cn_est_disc[id][1])
                    y_bi.append(data.cn_est_spanning[id][1])
            else:
                x_failed.append(data.cn_est_disc[id][1])
                y_failed.append(data.cn_est_spanning[id][1])
        
        ax=plt.subplot(gs[0,6])  # del, x=tsd
        sns_x=[ i for i in x if i < data.disc_thresholds[3] ]
        sns.violinplot(sns_x, orient='h', color='darkgreen')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_xlim(0, data.disc_thresholds[3])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        ax=plt.subplot(gs[1,7])  # tsd, y=spanning
        sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
        sns.violinplot(sns_x, orient='v', color='darkgreen')
        plt.setp(ax.collections, alpha=0.25)
        ax.set_ylim(-5, data.spanning_thresholds[2])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

        ax=plt.subplot(gs[1,6])  # x=discordant, y=spanning reads
        ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
        ax.scatter(x_mono, y_mono, s=5, c='lightgreen', linewidths=0.5, alpha=0.1)
        ax.scatter(x_bi, y_bi, s=5, c='darkolivegreen', linewidths=0.5, alpha=0.1)
        ax.axvline(x=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='forestgreen', linestyle='dashed')
        ax.axvline(x=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='forestgreen', linestyle='dashed')
        ax.axvline(x=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='forestgreen', linestyle='dashed')
        ax.axvline(x=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
        ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='forestgreen', linestyle='dashed')
        ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
        ax.set_xlabel('# discordant read')
        ax.set_ylabel('# spanning read')
        ax.set_xlim(0, data.disc_thresholds[3])
        ax.set_ylim(-5, data.spanning_thresholds[2])
        plt.suptitle('Genotyping result for insertions')
        plt.savefig(filenames.merged_pdf.replace('.pdf', '_gt.pdf'))
        plt.close()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def plot_custom(args, params, filenames, data):
    log.logger.debug('started')
    try:
        # load mei search filter result
        global mei_filter
        mei_filter={}
        with open(args.ins_bed) as infile:
            for line in infile:
                ls=line.strip().split()
                if ls[6] == 'confidence:high':
                    mei_filter[ls[10]]=True
                else:
                    mei_filter[ls[10]]=False
        # load genotyping results
        spanning_threshold_for_merge= 1 + math.ceil(args.cov * params.spanning_threshold_coeff_for_merge)
        if not data.disc_thresholds is False:
            plt.figure(figsize=(6.2, 4))  # (x, y)
            gs=gridspec.GridSpec(5, 8, height_ratios=[0.05, 0.3, 0.1, 0.05, 0.3], width_ratios=[0.3, 0.05, 0.1, 0.3, 0.05, 0.1, 0.3, 0.05])  # (y, x)
            gs.update(hspace=0.05, wspace=0.05)
            
            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
                    if mei_filter[id] is True:
                        if data.merged_res[id][0] == 1:
                            x_mono.append(data.cn_est_tsd_depth[id][1])
                            y_mono.append(data.cn_est_spanning[id][1])
                        else:
                            x_bi.append(data.cn_est_tsd_depth[id][1])
                            y_bi.append(data.cn_est_spanning[id][1])
                    else:
                        x_failed.append(data.cn_est_tsd_depth[id][1])
                        y_failed.append(data.cn_est_spanning[id][1])
            
            ax=plt.subplot(gs[0,0])  # tsd, x=tsd
            sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,1])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,0])  # tsd, x=tsd, y=spanning
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.1)
            ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.1)
#            ax.axvline(x=data.tsd_thresholds[4], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
#            ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
#            ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, TSD')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])
            
            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    if data.cn_est_spanning[id][1] < spanning_threshold_for_merge:
                        x.append(data.cn_est_tsd_depth[id][1])
                        y.append(data.cn_est_disc[id][1])
                        if mei_filter[id] is True:
                            if data.merged_res[id][0] == 1:
                                x_mono.append(data.cn_est_tsd_depth[id][1])
                                y_mono.append(data.cn_est_disc[id][1])
                            else:
                                x_bi.append(data.cn_est_tsd_depth[id][1])
                                y_bi.append(data.cn_est_disc[id][1])
                        else:
                            x_failed.append(data.cn_est_tsd_depth[id][1])
                            y_failed.append(data.cn_est_disc[id][1])
                    
            ax=plt.subplot(gs[3,0])  # tsd, x=tsd
            sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[4,1])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.disc_thresholds[3] ]
            sns.violinplot(sns_x, orient='v', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(0, data.disc_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[4,0])  # tsd, x=tsd, y=discordant reads
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.1)
            ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.1)
#            ax.axvline(x=data.tsd_thresholds[4], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
#            ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
#            ax.axhline(y=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
#            ax.axhline(y=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
#            ax.axhline(y=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, TSD')
            ax.set_ylabel('# discordant read')
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.set_ylim(0, data.disc_thresholds[3])

            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'Del':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
                    if mei_filter[id] is True:
                        if data.merged_res[id][0] == 1:
                            x_mono.append(data.cn_est_tsd_depth[id][1])
                            y_mono.append(data.cn_est_spanning[id][1])
                        else:
                            x_bi.append(data.cn_est_tsd_depth[id][1])
                            y_bi.append(data.cn_est_spanning[id][1])
                    else:
                        x_failed.append(data.cn_est_tsd_depth[id][1])
                        y_failed.append(data.cn_est_spanning[id][1])
            
            ax=plt.subplot(gs[0,3])  # del, x=tsd
            sns_x=[ i for i in x if i < data.del_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='darkred')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(-0.25, data.del_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,4])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='darkred')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])

            ax=plt.subplot(gs[1,3])  # del, x=tsd, y=spanning
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='gold', linewidths=0.5, alpha=0.2)
            ax.scatter(x_bi, y_bi, s=5, c='darkred', linewidths=0.5, alpha=0.2)
            ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
            ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, Del')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(-0.25, data.del_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])
            
            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'Del':
                    if data.cn_est_spanning[id][1] < spanning_threshold_for_merge:
                        x.append(data.cn_est_tsd_depth[id][1])
                        y.append(data.cn_est_disc[id][1])
                        if mei_filter[id] is True:
                            if data.merged_res[id][0] == 1:
                                x_mono.append(data.cn_est_tsd_depth[id][1])
                                y_mono.append(data.cn_est_disc[id][1])
                            else:
                                x_bi.append(data.cn_est_tsd_depth[id][1])
                                y_bi.append(data.cn_est_disc[id][1])
                        else:
                            x_failed.append(data.cn_est_tsd_depth[id][1])
                            y_failed.append(data.cn_est_disc[id][1])

            
            ax=plt.subplot(gs[3,3])  # del, x=tsd
            sns_x=[ i for i in x if i < data.del_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='darkred')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(-0.25, data.del_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[4,4])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.disc_thresholds[3] ]
            sns.violinplot(sns_x, orient='v', color='darkred')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(0, data.disc_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])

            ax=plt.subplot(gs[4,3])  # del, x=tsd, y=discordant reads
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='gold', linewidths=0.5, alpha=0.2)
            ax.scatter(x_bi, y_bi, s=5, c='darkred', linewidths=0.5, alpha=0.2)
            ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
            ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, Del')
            ax.set_ylabel('# discordant read')
            ax.set_xlim(-0.25, data.del_thresholds[3])
            ax.set_ylim(0, data.disc_thresholds[3])

            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_disc:
                x.append(data.cn_est_disc[id][1])
                y.append(data.cn_est_spanning[id][1])
                if mei_filter[id] is True:
                    if data.merged_res[id][0] == 1:
                        x_mono.append(data.cn_est_disc[id][1])
                        y_mono.append(data.cn_est_spanning[id][1])
                    else:
                        x_bi.append(data.cn_est_disc[id][1])
                        y_bi.append(data.cn_est_spanning[id][1])
                else:
                    x_failed.append(data.cn_est_disc[id][1])
                    y_failed.append(data.cn_est_spanning[id][1])
            
            ax=plt.subplot(gs[0,6])  # del, x=tsd
            sns_x=[ i for i in x if i < data.disc_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='darkgreen')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.disc_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,7])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='darkgreen')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])

            ax=plt.subplot(gs[1,6])  # x=discordant, y=spanning reads
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='lightgreen', linewidths=0.5, alpha=0.1)
            ax.scatter(x_bi, y_bi, s=5, c='darkolivegreen', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.disc_thresholds[0], linewidth=0.5, alpha=0.25, color='forestgreen', linestyle='dashed')
            ax.axvline(x=data.disc_thresholds[1], linewidth=0.5, alpha=0.50, color='forestgreen', linestyle='dashed')
            ax.axvline(x=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='forestgreen', linestyle='dashed')
            ax.axvline(x=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='forestgreen', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('# discordant read')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(0, data.disc_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])
        else:
            log.logger.warning('No discordant read stat available. Will use other evidences.')
            plt.figure(figsize=(4.1, 2))  # (x, y)
            gs=gridspec.GridSpec(2, 5, height_ratios=[0.05, 0.3], width_ratios=[0.3, 0.05, 0.1, 0.3, 0.05])  # (y, x)
            gs.update(hspace=0.05, wspace=0.05)
            
            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
                    if mei_filter[id] is True:
                        if data.merged_res[id][0] == 1:
                            x_mono.append(data.cn_est_tsd_depth[id][1])
                            y_mono.append(data.cn_est_spanning[id][1])
                        else:
                            x_bi.append(data.cn_est_tsd_depth[id][1])
                            y_bi.append(data.cn_est_spanning[id][1])
                    else:
                        x_failed.append(data.cn_est_tsd_depth[id][1])
                        y_failed.append(data.cn_est_spanning[id][1])
            
            ax=plt.subplot(gs[0,0])  # tsd, x=tsd
            sns_x=[ i for i in x if i < data.tsd_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,1])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='steelblue')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,0])  # tsd, x=tsd, y=spanning
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='lightskyblue', linewidths=0.5, alpha=0.1)
            ax.scatter(x_bi, y_bi, s=5, c='steelblue', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, TSD')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])
            
            x,y, x_mono,y_mono, x_bi,y_bi, x_failed,y_failed=[],[], [],[], [],[], [],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'Del':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
                    if mei_filter[id] is True:
                        if data.merged_res[id][0] == 1:
                            x_mono.append(data.cn_est_tsd_depth[id][1])
                            y_mono.append(data.cn_est_spanning[id][1])
                        else:
                            x_bi.append(data.cn_est_tsd_depth[id][1])
                            y_bi.append(data.cn_est_spanning[id][1])
                    else:
                        x_failed.append(data.cn_est_tsd_depth[id][1])
                        y_failed.append(data.cn_est_spanning[id][1])
            
            ax=plt.subplot(gs[0,3])  # del, x=tsd
            sns_x=[ i for i in x if i < data.del_thresholds[3] ]
            sns.violinplot(sns_x, orient='h', color='darkred')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_xlim(-0.25, data.del_thresholds[3])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            
            ax=plt.subplot(gs[1,4])  # tsd, y=spanning
            sns_x=[ i for i in y if i < data.spanning_thresholds[2] ]
            sns.violinplot(sns_x, orient='v', color='darkred')
            plt.setp(ax.collections, alpha=0.25)
            ax.set_ylim(-5, data.spanning_thresholds[2])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])

            ax=plt.subplot(gs[1,3])  # del, x=tsd, y=spanning
            ax.scatter(x_failed, y_failed, s=5, c='silver', linewidths=0.5, alpha=0.1)
            ax.scatter(x_mono, y_mono, s=5, c='gold', linewidths=0.5, alpha=0.2)
            ax.scatter(x_bi, y_bi, s=5, c='darkred', linewidths=0.5, alpha=0.2)
            ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
            ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=spanning_threshold_for_merge, linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
            ax.axhline(y=data.spanning_thresholds[2], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, Del')
            ax.set_ylabel('# spanning read')
            ax.set_xlim(-0.25, data.del_thresholds[3])
            ax.set_ylim(-5, data.spanning_thresholds[2])
        plt.suptitle('Genotyping result for insertions')
        plt.savefig(filenames.merged_pdf, transparent=True)
        plt.close()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
