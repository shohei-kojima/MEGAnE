#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,math,datetime
from utils import parse_fasta
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


def calc_dosage(args, params, filenames, infilename, plotfilename, outfilename):
    log.logger.debug('started')
    try:
        if args.monoallelic is True:
            log.logger.warning('Do not infer gene dosage, due to "-monoallelic" option. Exit.')
            exit(0)
        
        def gaussian_func_biallelics(coeff):
            def gaussian_func_biallelic(x, a, mu, sigma):
                return ((1-coeff)*a*np.exp(-(x-mu)**2/(2*sigma**2))) + (coeff*a*np.exp(-(x-(2*mu))**2/(4*sigma**2)))
            return gaussian_func_biallelic
        
        def gaussian_func(x, a, mu, sigma):
            return a*np.exp(-(x-mu)**2/(2*sigma**2))

        def fit_gaussian(list_support_read_count):
            x,y=[],[]
            support_read_bin= int(np.ceil(args.cov / 50))
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
        with open(infilename) as infile:
            for line in infile:
                ls=line.split('\t')
                chimeric_l=int(ls[4].split(',')[1].replace('chimeric=', '')) + int(ls[4].split(',')[2].replace('hybrid=', ''))
                chimeric_r=int(ls[5].split(',')[1].replace('chimeric=', '')) + int(ls[5].split(',')[2].replace('hybrid=', ''))
                if 'MEI_left_breakpoint=pT' in ls[8]:
                    count=chimeric_r
                elif 'MEI_right_breakpoint=pA' in ls[8]:
                    count=chimeric_l
                else:
                    count= math.ceil((chimeric_l + chimeric_r) / 2)
                if ls[6] == 'confidence:high' and 'unique:yes' in ls[7]:
                    for_gaussian_fitting.append(count)
                all_mei_count_range.append(count)
        if len(for_gaussian_fitting) < 10:
            log.logger.warning('Not enough data found. Cannot automatically determine thresholds for MEI filtering. Please check if your data contains MEIs enough for inferring gene dosage.')
        else:
            if args.b is not None:
                input_sample=os.path.basename(args.b)
            elif args.c is not None:
                input_sample=os.path.basename(args.c)
            input_bed=os.path.basename(infilename)
            x,y,popt,pcov,r_squared,coeff=fit_gaussian(for_gaussian_fitting)  # gaussian fitting
            if x is False:
                log.logger.debug('Gaussian curve fit failed. Exporting %s skipped.' % outfilename)
            else:
                log.logger.debug('popt=%s,pcov=%s' %(str(popt), str(pcov)))
                xd=np.arange(0, math.ceil(max(all_mei_count_range)) + 1)
                mono_sd= popt[2] ** 2  # standard deviation; sigma ** 2
                mono_zscore= (xd - popt[1]) / mono_sd  # zscore = (value - mu) / stdev
                di_zscore= (xd - (popt[1] * 2)) / (mono_sd * 2)
                dosage={}
                for pos,mono,di in zip(xd,mono_zscore,di_zscore):
                    if abs(mono) < abs(di):
                        dosage[pos]=1
                    else:
                        dosage[pos]=2
                mono_x,mono_y, di_x,di_y=[],[], [],[]
                for xval,yval in zip(x,y):
                    if dosage[math.ceil(xval)] == 1:
                        mono_x.append(xval)
                        mono_y.append(yval)
                    else:
                        di_x.append(xval)
                        di_y.append(yval)
                estimated_curve=gaussian_func_biallelics(coeff)(xd, popt[0], popt[1], popt[2])
                estimated_curve_single_allele=gaussian_func(xd, (1-coeff) * popt[0], popt[1], popt[2])
                estimated_curve_bi_allele=gaussian_func(xd, coeff * popt[0], 2 * popt[1], 1.414 * popt[2])
                
                # plot
                fig=plt.figure(figsize=(3,3))
                ax=fig.add_subplot(111)
                ax.scatter(mono_x, mono_y, s=5, c='dodgerblue', linewidths=0.5, alpha=0.5, label='Dosage=1')
                ax.scatter(di_x, di_y, s=5, c='coral', linewidths=0.5, alpha=0.5, label='Dosage=2')
                ax.plot(xd, estimated_curve_single_allele, color='grey', alpha=0.5)
                ax.plot(xd, estimated_curve_bi_allele, color='grey', alpha=0.5)
                ax.plot(xd, estimated_curve, label='Gaussian curve fitting', color='red', alpha=0.5)
                ax.set_xlim(0, popt[1] * 4)
                ax.set_xlabel('Number of chimeric + hybrid reads per breakpoint')
                ax.set_ylabel('Number of MEI')
                ax.legend()
                plt.suptitle('sample=%s;%s,\nn=%d, r_squared=%f,' % (input_sample, input_bed, len(for_gaussian_fitting), r_squared))  # popt[1] = mean, popt[2] = sigma
                plt.savefig(plotfilename)
                plt.close()
                log.logger.debug('gaussian_fitting_n=%d,r_squared=%f' %(len(for_gaussian_fitting), r_squared))
                
                # save results
                header=[]
                header.append('##fileformat=VCFv4.1\n')
                header.append('##fileDate=%s\n' % str(datetime.datetime.now()).split('.')[0])
                header.append('##source=MEI search version "%s"\n' % args.version)
                header.append('##reference=%s\n' % args.fa)
                with open(args.fai) as infile:
                    for line in infile:
                        ls=line.split()
                        if ls[0] in args.main_chrs_set:
                            header.append('##contig=<ID=%s,length=%s>\n' % (ls[0], ls[1]))
                header.append('##ALT=<ID=INS:ME,Description="Insertion of mobile element">\n')
                header.append('##FILTER=<ID=LC,Description="Low confidence">\n')
                header.append('##FILTER=<ID=NU,Description="Not unique">\n')
                header.append('##FILTER=<ID=S,Description="Shorter than 50-bp">\n')
                header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
                header.append('##INFO=<ID=MEPRED,Number=1,Type=String,Description="ME prediction status">\n')
                header.append('##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">\n')
                header.append('##INFO=<ID=MEI_rpos,Number=1,Type=Integer,Description="Reference poition of MEI right breakpoint">\n')
                header.append('##INFO=<ID=MEI,Number=.,Type=String,Description="Mobile element info">\n')
                header.append('##INFO=<ID=MEI_lbp,Number=.,Type=String,Description="Mobile element info of left breakpint">\n')
                header.append('##INFO=<ID=MEI_rbp,Number=.,Type=String,Description="Mobile element info of right breakpint">\n')
                header.append('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
                header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % input_sample)
                # retrieve ref seq
                bed=[]
                with open(infilename) as infile:
                    for line in infile:
                        ls=line.strip().split('\t')
                        bed.append('%s\t%s\t%d\t%s\n' % (ls[0], ls[1], int(ls[1]) + 1, ls[-1]))
                bed=BedTool(''.join(bed), from_string=True)
                fa=bed.sequence(fi=args.fa, name=True)
                tmp=parse_fasta(fa.seqfn)
                fa={}
                for h in tmp:
                    fa[h.split('::')[0]]=tmp[h].upper()
                # load ME length
                with open(outfilename, 'w') as outfile:  # vcf file
                    outfile.write(''.join(header))
                    with open(infilename) as infile:
                        for line in infile:
                            ls=line.strip().split('\t')
                            chimeric_l=int(ls[4].split(',')[1].replace('chimeric=', ''))
                            chimeric_r=int(ls[5].split(',')[1].replace('chimeric=', ''))
                            if 'MEI_left_breakpoint=pT' in ls[8]:
                                count=chimeric_r
                            elif 'MEI_right_breakpoint=pA' in ls[8]:
                                count=chimeric_l
                            else:
                                count= math.ceil((chimeric_l + chimeric_r) / 2)
                            start= int(ls[1]) + 1  # 0-based to 1-based
                            id=ls[-1].replace('ID=', '')
                            id=id[:-1] if id[-1] == ';' else id
                            seq=fa['>%s' % ls[-1]]
                            if ls[6] == 'confidence:high' and 'unique:yes' in ls[7]:
                                filt='PASS'
                            elif ls[6] == 'confidence:high' and 'unique:no' in ls[7]:
                                filt='NU'
                            elif ls[6] == 'confidence:low' and 'unique:yes' in ls[7]:
                                if '50bp_or_longer:yes' in ls[7]:
                                    filt='LC'
                                else:
                                    filt='LC;S'
                            else:
                                if '50bp_or_longer:yes' in ls[7]:
                                    filt='LC;NU'
                                else:
                                    filt='LC;NU;S'
                            left= int(ls[4].split(',')[0].replace('MEI_left:ref_pos=', ''))
                            right= int(ls[5].split(',')[0].replace('MEI_right:ref_pos=', ''))
                            homlen= (left - right) if left > right else 0
                            meinfo=ls[8].split(',', 1)[1].replace('MEI_left_breakpoint', 'MEI_lbp').replace('MEI_right_breakpoint', 'MEI_rbp')
                            mepred='PASS' if 'subfamily_pred:status=PASS' in ls[8] else 'FAILED'
                            info='SVTYPE=%s;MEPRED=%s;HOMLEN=%d;MEI_rpos=%s;%s' % (ls[3], mepred, homlen, ls[2], meinfo)
                            outline='%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n' % (ls[0], start, id, seq, '<INS:ME>', '.', filt, info, 'CN', dosage[count])
                            outfile.write(outline)
                    outfile.flush()
                    os.fdatasync(outfile.fileno())
                
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


