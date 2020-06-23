#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os
from scipy import stats
import numpy as np
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


def plot_orig_simple(args, params, filenames, data):
    log.logger.debug('started')
    try:
        if not data.disc_thresholds is False:
            plt.figure(figsize=(6.2, 4))  # (x, y)
            gs=gridspec.GridSpec(2, 3)  # (y, x)
            gs.update(hspace=0.3, wspace=0.3)
            
            ax=plt.subplot(gs[0,0])  # tsd, x=tsd, y=spanning
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
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
            
            ax=plt.subplot(gs[1,0])  # tsd, x=tsd, y=discordant reads
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_disc[id][1])
            ax.scatter(x, y, s=5, c='dodgerblue', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.tsd_thresholds[0], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[1], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[2], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axvline(x=data.tsd_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[0], linewidth=0.5, alpha=0.50, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[1], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='steelblue', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, TSD')
            ax.set_ylabel('# discordant read')
            ax.set_xlim(0, data.tsd_thresholds[3])
            ax.set_ylim(0, data.disc_thresholds[3])

            ax=plt.subplot(gs[0,1])  # del, x=tsd, y=spanning
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'Del':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
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

            ax=plt.subplot(gs[1,1])  # del, x=tsd, y=discordant reads
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'Del':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_disc[id][1])
            ax.scatter(x, y, s=5, c='coral', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.del_thresholds[0], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
            ax.axvline(x=data.del_thresholds[1], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
            ax.axvline(x=data.del_thresholds[2], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
            ax.axvline(x=data.del_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[0], linewidth=0.5, alpha=0.50, color='orangered', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[1], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[2], linewidth=0.5, alpha=0.25, color='orangered', linestyle='dashed')
            ax.axhline(y=data.disc_thresholds[3], linewidth=0.5, alpha=0.25, color='grey', linestyle='dashed')
            ax.set_xlabel('Relative depth, Del')
            ax.set_ylabel('# discordant read')
            ax.set_xlim(-0.25, data.del_thresholds[3])
            ax.set_ylim(0, data.disc_thresholds[3])

            ax=plt.subplot(gs[0,2])  # x=discordant, y=spanning reads
            x,y=[],[]
            for id in data.cn_est_disc:
                x.append(data.cn_est_disc[id][1])
                y.append(data.cn_est_spanning[id][1])
            ax.scatter(x, y, s=5, c='seagreen', linewidths=0.5, alpha=0.1)
            ax.axvline(x=data.disc_thresholds[0], linewidth=0.5, alpha=0.50, color='forestgreen', linestyle='dashed')
            ax.axvline(x=data.disc_thresholds[1], linewidth=0.5, alpha=0.25, color='forestgreen', linestyle='dashed')
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
            plt.figure(figsize=(2, 4))  # (x, y)
            gs=gridspec.GridSpec(2, 1)  # (y, x)
            
            ax=plt.subplot(gs[0])  # tsd, x=tsd, y=spanning
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
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
            ax.set_ylim(0, data.spanning_thresholds[2])
            
            ax=plt.subplot(gs[1])  # del, x=tsd, y=spanning
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'Del':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
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
            ax.set_xlim(0, data.del_thresholds[3])
            ax.set_ylim(0, data.spanning_thresholds[2])
        plt.suptitle('This is a figure for debug')
        plt.savefig(filenames.debug_pdf1)
        plt.close()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def plot_orig(args, params, filenames, data):
    log.logger.debug('started')
    try:
        if not data.disc_thresholds is False:
            plt.figure(figsize=(6.2, 4))  # (x, y)
            gs=gridspec.GridSpec(5, 8, height_ratios=[0.05, 0.3, 0.1, 0.05, 0.3], width_ratios=[0.3, 0.05, 0.1, 0.3, 0.05, 0.1, 0.3, 0.05])  # (y, x)
            gs.update(hspace=0.1, wspace=0.1)

            
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
            plt.figure(figsize=(2, 4))  # (x, y)
            gs=gridspec.GridSpec(2, 1)  # (y, x)
            
            ax=plt.subplot(gs[0])  # tsd, x=tsd, y=spanning
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'TSD':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
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
            ax.set_ylim(0, data.spanning_thresholds[2])
            
            ax=plt.subplot(gs[1])  # del, x=tsd, y=spanning
            x,y=[],[]
            for id in data.cn_est_tsd_depth:
                if data.cn_est_tsd_depth[id][2] == 'Del':
                    x.append(data.cn_est_tsd_depth[id][1])
                    y.append(data.cn_est_spanning[id][1])
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
            ax.set_xlim(0, data.del_thresholds[3])
            ax.set_ylim(0, data.spanning_thresholds[2])
        plt.suptitle('This is a figure for debug')
        plt.savefig(filenames.debug_pdf1)
        plt.close()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
