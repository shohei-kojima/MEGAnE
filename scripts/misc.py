#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

# this is a miscellaneous scripts not related to the main workflow
# most are for debugging

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import log,traceback


def plot_tsd_dep(left, middle, right):  # 20-nt flank
    index=[]
    for n in range(10):
        index.append('left_%d' % (10 - n))
    df=pd.DataFrame(left, columns=index)
    df['tsd_mean']=middle
    index=[]
    for n in range(10):
        index.append('right_%d' % (n + 1))
    tmp=pd.DataFrame(right, columns=index)
    for k in tmp:
        df[k]=tmp[k]
#    print(df)
    # plot
    fig=plt.figure(figsize=(5, 3))  # (x, y)
    ax=fig.add_subplot(111)
#    sns.violinplot(data=df, width=1, fliersize=0, gridsize=1000, palette='cool')
#    plt.setp(ax.collections, alpha=0.25)
    sns.boxplot(data=df, width=0.75, showfliers=False, palette='cool', boxprops=dict(alpha=.3))
    ax.set_ylim(0, 2.5)
    ax.set_xlabel('Position relative to TSD')
    ax.set_ylabel('Depth relative to autosome')
    plt.xticks(rotation=90)
    plt.suptitle('misc/plot_tsd_dep')
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.9)
    plt.savefig('./genotype_out/plot_out_misc_tsd_dep.pdf')



