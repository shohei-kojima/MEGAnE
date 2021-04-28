#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,math
from utils import parse_fasta
import pybedtools
from pybedtools import BedTool
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import log,traceback


matplotlib.rcParams['lines.linewidth']=0.5
matplotlib.rcParams['axes.linewidth']=0.5
matplotlib.rcParams['xtick.major.width']=0.5
matplotlib.rcParams['ytick.major.width']=0.5
matplotlib.rcParams['font.size']=5


def filter(args, params, filenames):
    log.logger.debug('started')
    try:
        pybedtools.set_tempdir(args.pybedtools_tmp)
        
        if args.b is not None:
            input_sample=os.path.basename(args.b)
        elif args.c is not None:
            input_sample=os.path.basename(args.c)
        nts=['A', 'T']
        
        # determine cutoff from actual percentile threshold
        def determine_cutoff(list_support_read_count):
            cutoff_pos= math.ceil(len(list_support_read_count) * params.actual_cutoff_rank)
            cutoff=sorted(list_support_read_count)[cutoff_pos]
            return cutoff

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
            if args.monoallelic is False:
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
                r_squared=fits[-1][0]
                r_squared_few=fits[-1][1]
                coeff=fits[-1][2]
                popt=fits[-1][3]
                pcov=fits[-1][4]
                log.logger.debug('r_squared=%f,biallelic_coeff=%f' %(r_squared, coeff))
            elif args.monoallelic is True:
                coeff=None
                popt,pcov=curve_fit(gaussian_func, x, y, init_param)
                residuals= y - gaussian_func(x, *popt)
                rss=np.sum(residuals**2)
                tss=np.sum((y - np.mean(y))**2)
                r_squared= 1 - (rss / tss)
            reject1perc=norm.interval(alpha=params.fit_gaussian_CI_alpha, loc=popt[1], scale=abs(popt[2]))
            return x, y, popt, pcov, reject1perc, r_squared, coeff

        def L1_filter(line, r_pos, l_pos, pA_only, R_eval, L_eval):
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
        hybrid_num_threshold= args.cov * params.hybrid_read_coeff_for_gaussian_fitting
        chimeric_num_threshold= args.cov * params.chimeric_read_coeff_for_gaussian_fitting
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
                if (int(ls[12]) >= hybrid_num_threshold) and (int(ls[13]) >= hybrid_num_threshold):
                    if (len(R_eval) >= chimeric_num_threshold) and (len(L_eval) >= chimeric_num_threshold):
                        if (min(R_eval) < params.eval_threshold_for_gaussian_fitting) and (min(L_eval) < params.eval_threshold_for_gaussian_fitting):
                            for_gaussian_fitting.append(total_read_count)
                            hybrid_num.append(int(ls[12]) + int(ls[13]))
                        elif not ls[5] == 'NA':
                            L_eval_min=min(L_eval)
                            if L_eval_min < params.eval_threshold_for_gaussian_fitting:
                                for_gaussian_fitting.append(total_read_count)
                                hybrid_num.append(int(ls[12]) + int(ls[13]))
                        elif not ls[6] == 'NA':
                            R_eval_min=min(R_eval)
                            if R_eval_min < params.eval_threshold_for_gaussian_fitting:
                                for_gaussian_fitting.append(total_read_count)
                                hybrid_num.append(int(ls[12]) + int(ls[13]))
                    elif not ls[5] == 'NA':
                        if len(L_eval) >= chimeric_num_threshold:
                            L_eval_min=min(L_eval)
                            if L_eval_min < params.eval_threshold_for_gaussian_fitting:
                                for_gaussian_fitting.append(total_read_count)
                                hybrid_num.append(int(ls[12]) + int(ls[13]))
                    elif not ls[6] == 'NA':
                        if len(R_eval) >= chimeric_num_threshold:
                            R_eval_min=min(R_eval)
                            if R_eval_min < params.eval_threshold_for_gaussian_fitting:
                                for_gaussian_fitting.append(total_read_count)
                                hybrid_num.append(int(ls[12]) + int(ls[13]))

        global gaussian_executed
        if len(for_gaussian_fitting) >= 3:
            gaussian_executed=True
            cutoff=determine_cutoff(for_gaussian_fitting)  # percentile cutoff
            log.logger.debug('total_support_read_num_cutoff=%d' % cutoff)
            
            x,y,popt,pcov,reject1perc,r_squared,coeff=fit_gaussian(for_gaussian_fitting)  # gaussian cutoff
            log.logger.debug('popt=%s,pcov=%s' %(str(popt), str(pcov)))
            xd=np.arange(min(x), max(x), 0.5)
            if args.monoallelic is False:
                estimated_curve=gaussian_func_biallelics(coeff)(xd, popt[0], popt[1], popt[2])
                estimated_curve_single_allele=gaussian_func(xd, (1-coeff) * popt[0], popt[1], popt[2])
                estimated_curve_bi_allele=gaussian_func(xd, coeff * popt[0], 2 * popt[1], 1.414 * popt[2])
            elif args.monoallelic is True:
                estimated_curve=gaussian_func(xd, popt[0], popt[1], popt[2])
            # plot
            total_read_thresholds=[reject1perc[0], cutoff]  # parameter setting
            if args.no_pdf is False:
                fig=plt.figure(figsize=(3,3))
                ax=fig.add_subplot(111)
                ax.scatter(x, y, s=5, c='dodgerblue', linewidths=0.5, alpha=0.5, label='Actual data')
                if args.monoallelic is False:
                    ax.plot(xd, estimated_curve_single_allele, color='grey', alpha=0.5)
                    ax.plot(xd, estimated_curve_bi_allele, color='grey', alpha=0.5)
                ax.plot(xd, estimated_curve, label='Gaussian curve fitting', color='red', alpha=0.5)
                ax.axvline(x=round(reject1perc[0]), linewidth=1, alpha=0.5, color='olive', linestyle='dashed', label='Gaussian_cutoff')
                ax.axvline(x=round(cutoff), linewidth=1, alpha=0.5, color='darkgreen', linestyle='dashed', label='Percentile_cutoff')
                ax.set_xlim(0, popt[1] * 4)
                ax.set_xlabel('Number of support reads per MEI')
                ax.set_ylabel('Number of MEI')
                ax.legend()
                plt.suptitle('sample=%s,\nn=%d, r_squared=%f,\ngaussian_cutoff=%f, percentile_cutoff=%d' % (input_sample, len(for_gaussian_fitting), r_squared, total_read_thresholds[0], int(total_read_thresholds[1])))  # popt[1] = mean, popt[2] = sigma
                plt.savefig(filenames.gaussian_plot)
                plt.close()
            zero_hybrid_total_read_threshold= round(popt[1] * ((sum(for_gaussian_fitting) - sum(hybrid_num)) / sum(for_gaussian_fitting)))  # parameter setting
            log.logger.debug('gaussian_fitting_n=%d,r_squared=%f,gaussian_cutoff=%f,percentile_cutoff=%f,zero_hybrid_total_read_threshold=%f' %(len(for_gaussian_fitting), r_squared, total_read_thresholds[0], total_read_thresholds[1], zero_hybrid_total_read_threshold))
        else:
            log.logger.warning('Not enough data found. Cannot automatically determine thresholds for MEI filtering. Please check if your data contains discordant reads enough for auto-filtering. Proceed anyway.')
            gaussian_executed=False
            total_read_thresholds=[3]
            zero_hybrid_total_read_threshold=5
            log.logger.debug('cutoff=%d,zero_hybrid_total_read_threshold=%d' %(total_read_thresholds[0], zero_hybrid_total_read_threshold))

        # main
        def main_filter(total_read_threshold, outfilename):
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
                            if args.L1_filt_off == False and ls[7] == 'L1':
                                L1_judge=L1_filter(line, r_pos, l_pos, pA_only, R_eval, L_eval)
                                pybedtools.cleanup()
                                if L1_judge is True:
                                    high.add(line)
                            else:
                                high.add(line)
                        elif (len(R_eval) >= 1) and (len(L_eval) >= 1):
                            proc_next=False
                            if args.lowdep is True or args.verylowdep is True:
                                if min(R_eval + L_eval) < params.second_filter_eval_threshold_for_few_hybrid and (int(ls[12]) + int(ls[13])) >= params.second_filter_hybrid_read_num:
                                    proc_next=True
                            else:
                                if (min(R_eval) < params.second_filter_eval_threshold_for_few_hybrid) and (min(L_eval) < params.second_filter_eval_threshold_for_few_hybrid):
                                    proc_next=True
                            if proc_next is True:
                                if args.L1_filt_off == False and ls[7] == 'L1':
                                    L1_judge=L1_filter(line, r_pos, l_pos, pA_only, R_eval, L_eval)
                                    pybedtools.cleanup()
                                    if L1_judge is True:
                                        high.add(line)
                                else:
                                    high.add(line)
                        elif args.lowdep is True or args.verylowdep is True:
                            RL_eval= R_eval + L_eval
                            proc_next=False
                            if args.lowdep is True:
                                if min(RL_eval) < params.second_filter_eval_threshold_for_few_hybrid and (int(ls[12]) + int(ls[13])) >= (2 * params.second_filter_hybrid_read_num):
                                    proc_next=True
                            elif args.verylowdep is True:
                                if min(RL_eval) < params.second_filter_eval_threshold_for_few_hybrid:
                                    proc_next=True
                            if proc_next is True:
                                if args.L1_filt_off == False and ls[7] == 'L1':
                                    L1_judge=L1_filter(line, r_pos, l_pos, pA_only, R_eval, L_eval)
                                    pybedtools.cleanup()
                                    if L1_judge is True:
                                        high.add(line)
                                else:
                                    high.add(line)
            
            with open(outfilename, 'w') as outfile:
                for line in all:
                    if line in high:
                        line=line.strip()
                        outfile.write(line +'\thigh\n')
                    else:
                        line=line.strip()
                        outfile.write(line +'\tlow\n')
                outfile.flush()
                os.fdatasync(outfile.fileno())

        # main filtering
        if gaussian_executed is True:
            main_filter(total_read_thresholds[0], filenames.bp_merged_filt_g)
            main_filter(total_read_thresholds[1], filenames.bp_merged_filt_p)
        else:
            main_filter(total_read_thresholds[0], filenames.bp_merged_filt_f)
        if args.threshold is not None:
            main_filter(args.threshold, filenames.bp_merged_filt_u)
        pybedtools.cleanup()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def grouping(args, filenames):
    log.logger.debug('started')
    try:
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

        # main grouping; group non-singletons
        def main_grouping(infilename, outfilename):
            L,R={},{}
            for chr in args.main_chrs_set:
                L[chr]={}
                R[chr]={}
            with open(infilename) as infile:
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

            with open(outfilename, 'w') as outfile:
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
                outfile.flush()
                os.fdatasync(outfile.fileno())
            return ins_n
                    
        global ins_ns
        ins_ns=[]
        if args.gaussian_executed is True:
            ins_n=main_grouping(filenames.bp_merged_filt_g, filenames.bp_merged_groupg)
            ins_ns.append(ins_n)
            ins_n=main_grouping(filenames.bp_merged_filt_p, filenames.bp_merged_groupp)
            ins_ns.append(ins_n)
        else:
            ins_n=main_grouping(filenames.bp_merged_filt_f, filenames.bp_merged_groupf)
            ins_ns.append(ins_n)
        if args.threshold is not None:
            ins_n=main_grouping(filenames.bp_merged_filt_u, filenames.bp_merged_groupu)
            ins_ns.append(ins_n)

    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

