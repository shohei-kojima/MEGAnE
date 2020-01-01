#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os
from utils import parse_fasta
from pybedtools import BedTool
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['lines.linewidth']=0.5
matplotlib.rcParams['axes.linewidth']=0.5
matplotlib.rcParams['xtick.major.width']=0.5
matplotlib.rcParams['ytick.major.width']=0.5
matplotlib.rcParams['font.size']=5


def filter(args, params, filenames):
    if args.b is not None:
        input_sample=os.path.basename(args.b)
    elif args.c is not None:
        input_sample=os.path.basename(args.c)
    nts=['A', 'T']
    
    def gaussian_func_biallelic(x, a, mu, sigma):
        return (a*np.exp(-(x-mu)**2/(2*sigma**2))) + (0.333*a*np.exp(-(x-(2*mu))**2/(4*sigma**2)))

    def gaussian_func(x, a, mu, sigma):
        return a*np.exp(-(x-mu)**2/(2*sigma**2))

    def fit_gaussian(list_support_read_count):
        x,y=[],[]
        for i in range(max(list_support_read_count)):
            x.append(i)
            y.append(list_support_read_count.count(i))
        init_param=[args.cov * params.fit_gaussian_init_a_coeff, args.cov * params.fit_gaussian_init_a_coeff, args.cov * params.fit_gaussian_init_a_coeff]
        popt,pcov=curve_fit(gaussian_func_biallelic, x, y, init_param)
        CI99=norm.interval(alpha=params.fit_gaussian_CI_alpha, loc=popt[1], scale=popt[2])
        return x, y, popt, pcov, CI99

    def L1_filter(line, r_pos, l_pos):
        cand=True
        if pA_only is False:
            ls=line.split()
            bed=''
            if (r_pos - l_pos) >= params.L1_filter_min_TSD_len:
                bed += ls[0] +'\t'+ str(l_pos) +'\t'+ str(r_pos) +'\n'
                bed=BedTool(bed, from_string=True)
                fa=bed.sequence(fi=args.fa)
                fa=parse_fasta(fa.seqfn)
                for h in fa:
                    seq=fa[h].upper()
                    seqlen=len(seq)
                    total_AT=0
                    for nt in nts:
                        c=seq.count(nt)
                        if (100 * (c / seqlen)) >= params.L1_filter_A_or_T_perc:
                            cand=False
                        total_AT += c
                    if (100 * (total_AT / seqlen)) >= params.L1_filter_A_plus_T_perc:
                        cand=False
        if len(R_eval) >= 1:
            if not (min(R_eval) < params.L1_filter_eval_threshold):
                cand=False
        if len(L_eval) >= 1:
            if not (min(L_eval) < params.L1_filter_eval_threshold):
                cand=False
        return cand

    # determine threshold by fitting gaussian function
    for_gaussian_fitting=[]
    hybrid_num=[]
    with open(filenames.bp_merged_all) as infile:
        for line in infile:
            ls=line.split('\t')
            r_pos,r_num=ls[3].split(':')
            l_pos,l_num=ls[4].split(':')
            r_num=0 if r_num == 'NA' else int(r_num)
            l_num=0 if l_num == 'NA' else int(l_num)
            total_read_count= r_num + l_num
            R_eval,L_eval=[],[]
            vs=ls[8].split(';')
            vs=[ float(v) for v in vs if not (v == 'NA') and not (v == '') ]
            if len(vs) >= 1:
                R_eval.extend(vs)
            vs=ls[9].split(';')
            vs=[ float(v) for v in vs if not (v == 'NA') and not (v == '') ]
            if len(vs) >= 1:
                L_eval.extend(vs)
            if (int(ls[12]) >= (args.cov * params.hybrid_read_coeff_for_gaussian_fitting)) and (int(ls[13]) >= (args.cov * params.hybrid_read_coeff_for_gaussian_fitting)):
                if (len(R_eval) >= 1) and (len(L_eval) >= 1):
                    if (min(R_eval) < params.eval_threshold_for_gaussian_fitting) and (min(L_eval) < params.eval_threshold_for_gaussian_fitting):
                        for_gaussian_fitting.append(total_read_count)
                        hybrid_num.append(int(ls[12]) + int(ls[13]))
    if len(for_gaussian_fitting) >= 3:
        x,y,popt,pcov,CI99=fit_gaussian(for_gaussian_fitting)
        xd=np.arange(min(x), max(x), 0.5)
        estimated_curve=gaussian_func_biallelic(xd, popt[0], popt[1], popt[2])
        estimated_curve_single_allele=gaussian_func(xd, popt[0], popt[1], popt[2])
        estimated_curve_bi_allele=gaussian_func(xd, 0.333 * popt[0], 2 * popt[1], 1.414 * popt[2])
        # plot
        fig=plt.figure(figsize=(3,3))
        ax=fig.add_subplot(111)
        ax.scatter(x, y, s=5, c='dodgerblue', linewidths=0.5, alpha=0.5, label='Actual data')
        ax.plot(xd, estimated_curve_single_allele, color='grey', alpha=0.5)
        ax.plot(xd, estimated_curve_bi_allele, color='grey', alpha=0.5)
        ax.plot(xd, estimated_curve, label='Gaussian curve fitting', color='red', alpha=0.5)
        ax.axvline(x=round(CI99[0]), linewidth=1, alpha=0.5, color='green', linestyle='dashed', label='Auto-determined threshold')
        ax.set_xlim(0, popt[1] * 4)
        ax.set_xlabel('Number of support reads per MEI')
        ax.set_ylabel('Number of MEI')
        ax.legend()
        plt.suptitle('sample=%s,\nn=%d, mu=%f, sigma=%f,\nlowCI99=%f' % (input_sample, len(for_gaussian_fitting), popt[1], popt[2], CI99[0]))  # popt[1] = mean, popt[2] = sigma
        plt.savefig(filenames.gaussian_plot)
        plt.close()
        # parameter settings
        total_read_threshold= round(CI99[0])
        zero_hybrid_total_read_threshold= round(popt[1] * ((sum(for_gaussian_fitting) - sum(hybrid_num)) / sum(for_gaussian_fitting)))
    else:
        print('AIM-UP: Warning: Not enough data to automatically determine thresholds for MEI filtering. Please check if your data is enough big or contans discordant reads. Proceed anyway.')
        total_read_threshold=3
        zero_hybrid_total_read_threshold=5

    # main
    all=[]
    high=set()
    with open(filenames.bp_merged_all) as infile:
        for line in infile:
            ls=line.split('\t')
            r_pos,r_num=ls[3].split(':')
            l_pos,l_num=ls[4].split(':')
            r_num=0 if r_num == 'NA' else int(r_num)
            l_num=0 if l_num == 'NA' else int(l_num)
            r_pos=0 if r_pos == 'NA' else int(r_pos)
            l_pos=0 if l_pos == 'NA' else int(l_pos)
            total_read_count= r_num + l_num
            R_eval,L_eval=[],[]
            vs=ls[8].split(';')
            vs=[ float(v) for v in vs if not (v == 'NA') and not (v == '') ]
            if len(vs) >= 1:
                R_eval.extend(vs)
            vs=ls[9].split(';')
            vs=[ float(v) for v in vs if not (v == 'NA') and not (v == '') ]
            if len(vs) >= 1:
                L_eval.extend(vs)
            retain_count=False
            if total_read_count >= total_read_threshold:
                retain_count=True
            retain_eval=False
            for l in [R_eval, L_eval]:
                if len(l) >= 1:
                    vmin=min(l)
                    if vmin <= params.first_filter_eval_threshold:
                        retain_eval=True
            pA_only=True
            if ls[5] == 'NA':
                if ls[6] == 'NA':
                    pA_only=False
            elif not ls[6] == 'NA':
                pA_only=False

            # first filter
            if (retain_count is True) and (retain_eval is True):
                if (int(ls[12]) + int(ls[13])) >= params.first_filter_total_hybrid_read_num:
                    all.append(line)
                elif total_read_count >= zero_hybrid_total_read_threshold:
                    all.append(line)
            # second filter, 90% accuracy
            if (retain_count is True) and (retain_eval is True):
                if (int(ls[12]) >= params.second_filter_hybrid_read_num) and (int(ls[13]) >= params.second_filter_hybrid_read_num):
                    if 'L1' in line:
                        L1_judge=L1_filter(line, r_pos, l_pos)
                        if L1_judge is True:
                            high.add(line)
                    else:
                        high.add(line)
                elif (len(R_eval) >= 1) and (len(L_eval) >= 1):
                    if (min(R_eval) < params.second_filter_eval_threshold_for_few_hybrid) and (min(L_eval) < params.second_filter_eval_threshold_for_few_hybrid):
                        if 'L1' in line:
                            L1_judge=L1_filter(line, r_pos, l_pos)
                            if L1_judge is True:
                                high.add(line)
                        else:
                            high.add(line)
    
    with open(filenames.bp_merged_filt, 'w') as outfile:
        for line in all:
            if line in high:
                line=line.strip()
                outfile.write(line +'\thigh\n')
            else:
                line=line.strip()
                outfile.write(line +'\tlow\n')


def grouping(args, filenames):
    def merge(one_l, other_l):
        orig_len=len(one_l)
        tmp=[]
        tmp.extend(one_l)
        unmerged=[]
        for a in other_l:
            merged=False
            for line in one_l:
                if line in a:
                    tmp.extend(a)
                    merged=True
                    break
            if merged is False:
                unmerged.append(a)
        tmp=list(set(tmp))
        tmp=sorted(tmp)
        after_len=len(tmp)
        if after_len > orig_len:
            merged=True
        else:
            merged=False
        return merged, tmp, unmerged

    # group non-singletons
    L,R={},{}
    for chr in args.main_chrs_set:
        L[chr]={}
        R[chr]={}
    with open(filenames.bp_merged_filt) as infile:
        for line in infile:
            line=line.strip()
            ls=line.split()
            r=ls[3].split(':')[0]
            l=ls[4].split(':')[0]
            L[ls[0]][l]=line
            R[ls[0]][r]=line

    all_0=[]
    with open(filenames.overhang_MEI) as infile:
        for line in infile:
            ls=line.split()
            poss=ls[0].split(';')[:-1]
            tmp_s=[]
            for p in poss:
                chr,tmp=p.split(':', 1)
                start,tmp=tmp.split('-', 1)
                end,dir,_=tmp.split('/', 2)
                if dir == 'L':
                    if start in L[chr]:
                        tmp_s.append(L[chr][start])
                else:
                    if end in R[chr]:
                        tmp_s.append(R[chr][end])
            all_0.append(tmp_s)

    all_1=[]
    for a in all_0:
        if not a in all_1:
            all_1.append(a)
    del(all_0)

    unmerged=all_1
    unmerged_len=len(unmerged)
    final=[]
    while unmerged_len >= 1:
        b=True
        l=unmerged[0]
        unmerged=unmerged[1:]
        while b is True:
            b,l,unmerged=merge(l, unmerged)
        final.append(l)
        unmerged_len=len(unmerged)

    singletons,multis={},{}
    me_clas=set()
    global ins_n
    ins_n=0
    n=1
    for f in final:
        if len(f) == 1:
            ls=f[0].split()
            if ls[-1] == 'high':
                ins_n += 1
            ls.append('singleton')
            if not ls[7] in singletons:
                singletons[ls[7]]=[]
                me_clas.add(ls[7])
            singletons[ls[7]].append(ls)
        else:
            groupname='group%d' % n
            for line in f:
                ls=line.split()
                ls.append(groupname)
                if not ls[7] in multis:
                    multis[ls[7]]={}
                    me_clas.add(ls[7])
                if not groupname in multis[ls[7]]:
                    multis[ls[7]][groupname]=[]
                multis[ls[7]][groupname].append(ls)
            n += 1
    me_clas=sorted(list(me_clas))

    with open(filenames.bp_merged_group, 'w') as outfile:
        for m in me_clas:
            if m in singletons:
                singletons[m]=sorted(singletons[m], key=lambda x:(x[0], int(x[1]), int(x[2])))
                for ls in singletons[m]:
                    outfile.write('\t'.join(ls) +'\n')
        for m in me_clas:
            if m in multis:
                for g in multis[m]:
                    multis[m][g]=sorted(multis[m][g], key=lambda x:(x[0], int(x[1]), int(x[2])))
                    for ls in multis[m][g]:
                        outfile.write('\t'.join(ls) +'\n')

