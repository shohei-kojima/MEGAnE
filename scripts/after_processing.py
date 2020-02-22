#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os
import math
from statistics import mean
import log,traceback


def grouped_mei_to_bed(filenames):
    log.logger.debug('started')
    
    try:
        def predict_shared_te(R_list, L_list, r_strand, l_strand):
            cands=[]
            R_tes,L_tes=set(),set()
            for l in R_list:
                R_tes.add(l[0])
            for l in L_list:
                L_tes.add(l[0])
            shared_tes= R_tes & L_tes
            if len(shared_tes) >= 1:
                tes=sorted(list(shared_tes))
                d={}
                for te in tes:
                    d[te]=[[], []]
                for te,s,_,_ in R_list:
                    d[te][0].append(int(s))
                for te,_,e,_ in L_list:
                    d[te][1].append(int(e))
                for te in d:
                    d[te][0]=round(mean(d[te][0]))
                    d[te][1]=round(mean(d[te][1]))
                    cands.append('%s,%d/%d,%s/%s' % (te, d[te][0], d[te][1], r_strand, l_strand))
            return cands

        def convert_to_bed(infilename, outfilename):
            with open(outfilename, 'w') as outfile:
                with open(infilename) as infile:
                    for line in infile:
                        ls=line.split()
                        # if breakpoints are not pA
                        if not (ls[8] == 'NA') and not (ls[9] == 'NA'):
                            # find shared TE between R and L
                            R_tes,L_tes=set(),set()
                            for tes in ls[10].split(';;'):
                                for te in tes.split(';'):
                                    if not te == '':
                                        ts=te.split(',')
                                        R_tes.add(ts[0])
                            for tes in ls[11].split(';;'):
                                for te in tes.split(';'):
                                    if not te == '':
                                        ts=te.split(',')
                                        L_tes.add(ts[0])
                            shared_tes= R_tes & L_tes
                            # if shared TEs exist
                            if len(shared_tes) >= 1:
                                evals_d={}
                                for te in shared_tes:
                                    evals_d[te]=[]
                                for eval,tes in zip(ls[8].split(';'), ls[10].split(';;')):
                                    if not eval == 'NA':
                                        for te in tes.split(';'):
                                            if not te == '':
                                                te_name=te.split(',')[0]
                                                if te_name in evals_d:
                                                    evals_d[te_name].append(float(eval))
                                for eval,tes in zip(ls[9].split(';'), ls[11].split(';;')):
                                    if not eval == 'NA':
                                        for te in tes.split(';'):
                                            if not te == '':
                                                te_name=te.split(',')[0]
                                                if te_name in evals_d:
                                                    evals_d[te_name].append(float(eval))
                                for te in evals_d:
                                    evals_d[te]=min(evals_d[te])
                                min_eval=min(list(evals_d.values()))
                                tes_min_eval=[]
                                for te in evals_d:
                                    if evals_d[te] == min_eval:
                                        tes_min_eval.append(te)
                                tes_min_eval=sorted(tes_min_eval)
                                R_plus,L_plus=[],[]
                                R_minus,L_minus=[],[]
                                # R breakpoint
                                for tes in ls[10].split(';;'):
                                    for te in tes.split(';'):
                                        if not te == '':
                                            ts=te.split(',')
                                            if ts[0] in tes_min_eval:
                                                if ts[3] == '+':
                                                    R_plus.append(ts)
                                                else:
                                                    R_minus.append(ts)
                                # L breakpoint
                                for tes in ls[11].split(';;'):
                                    for te in tes.split(';'):
                                        if not te == '':
                                            ts=te.split(',')
                                            if ts[0] in tes_min_eval:
                                                if ts[3] == '+':
                                                    L_plus.append(ts)
                                                else:
                                                    L_minus.append(ts)
                                cands=[]
                                if (len(R_plus) >= 1) and (len(L_plus) >= 1):
                                    tmp=predict_shared_te(R_plus, L_plus, '+', '+')
                                    cands.extend(tmp)
                                elif (len(R_minus) >= 1) and (len(L_minus) >= 1):
                                    tmp=predict_shared_te(R_minus, L_minus, '-', '-')
                                    cands.extend(tmp)
                                elif (len(R_plus) >= 1) and (len(L_minus) >= 1):
                                    tmp=predict_shared_te(R_plus, L_minus, '+', '-')
                                    cands.extend(tmp)
                                elif (len(R_minus) >= 1) and (len(L_plus) >= 1):
                                    tmp=predict_shared_te(R_minus, L_plus, '-', '+')
                                    cands.extend(tmp)
                                pred_status='PASS'
                                pred_res='MEI=' + '|'.join(cands)
                            # if not shared TE exist
                            else:
                                R_plus,L_plus=[],[]
                                R_minus,L_minus=[],[]
                                # R breakpoint
                                evals=[ float(i) for i in ls[8].split(';') ]
                                min_eval=min(evals)
                                for eval,tes in zip(ls[8].split(';'), ls[10].split(';;')):
                                    if not eval == 'NA':
                                        for te in tes.split(';'):
                                            if not te == '':
                                                ts=te.split(',')
                                                if float(eval) == min_eval:
                                                    if ts[3] == '+':
                                                        R_plus.append(ts)
                                                    else:
                                                        R_minus.append(ts)
                                # L breakpoint
                                evals=[ float(i) for i in ls[9].split(';') ]
                                min_eval=min(evals)
                                for eval,tes in zip(ls[9].split(';'), ls[11].split(';;')):
                                    if not eval == 'NA':
                                        for te in tes.split(';'):
                                            if not te == '':
                                                ts=te.split(',')
                                                if float(eval) == min_eval:
                                                    if ts[3] == '+':
                                                        L_plus.append(ts)
                                                    else:
                                                        L_minus.append(ts)
                                # reshape
                                R_str_l,L_str_l=[],[]
                                if len(R_plus) >= 1:
                                    d={}
                                    for te,s,_,_ in R_plus:
                                        if not te in d:
                                            d[te]=[]
                                        d[te].append(int(s))
                                    for te in d:
                                        breapoint=round(mean(d[te]))
                                        R_str_l.append('%s,%d,+' % (te, breapoint))
                                if len(R_minus) >= 1:
                                    d={}
                                    for te,s,_,_ in R_minus:
                                        if not te in d:
                                            d[te]=[]
                                        d[te].append(int(s))
                                    for te in d:
                                        breapoint=round(mean(d[te]))
                                        R_str_l.append('%s,%d,-' % (te, breapoint))
                                if len(L_plus) >= 1:
                                    d={}
                                    for te,_,e,_ in L_plus:
                                        if not te in d:
                                            d[te]=[]
                                        d[te].append(int(e))
                                    for te in d:
                                        breapoint=round(mean(d[te]))
                                        L_str_l.append('%s,%d,+' % (te, breapoint))
                                if len(L_minus) >= 1:
                                    d={}
                                    for te,_,e,_ in L_minus:
                                        if not te in d:
                                            d[te]=[]
                                        d[te].append(int(e))
                                    for te in d:
                                        breapoint=round(mean(d[te]))
                                        L_str_l.append('%s,%d,-' % (te, breapoint))
                                pred_status='Complex_structure'
                                pred_res='MEI_left_breakpoint=' + '|'.join(R_str_l) +';'+ 'MEI_right_breakpoint=' + '|'.join(L_str_l)
                        # if either end is pA
                        elif ls[9] == 'NA':
                            bp_plus,bp_minus=[],[]
                            # R breakpoint
                            evals=[ float(i) for i in ls[8].split(';') ]
                            min_eval=min(evals)
                            for eval,tes in zip(ls[8].split(';'), ls[10].split(';;')):
                                if not eval == 'NA':
                                    for te in tes.split(';'):
                                        if not te == '':
                                            ts=te.split(',')
                                            if float(eval) == min_eval:
                                                if ts[3] == '+':
                                                    bp_plus.append(ts)
                                                else:
                                                    bp_minus.append(ts)
                            # normal structure
                            if len(bp_plus) >= 1:
                                bp_str_l=[]
                                d={}
                                for te,s,_,_ in bp_plus:
                                    if not te in d:
                                        d[te]=[]
                                    d[te].append(int(s))
                                for te in d:
                                    breapoint=round(mean(d[te]))
                                    bp_str_l.append('%s,%d,+' % (te, breapoint))
                                pred_status='PASS'
                                pred_res='MEI_left_breakpoint=' + '|'.join(bp_str_l) +';'+ 'MEI_right_breakpoint=pA'
                            # complex structure
                            else:
                                bp_str_l=[]
                                d={}
                                for te,_,e,_ in bp_minus:
                                    if not te in d:
                                        d[te]=[]
                                    d[te].append(int(e))
                                for te in d:
                                    breapoint=round(mean(d[te]))
                                    bp_str_l.append('%s,%d,-' % (te, breapoint))
                                pred_status='complex_structure'
                                pred_res='MEI_left_breakpoint=' + '|'.join(bp_str_l) +';'+ 'MEI_right_breakpoint=pA'
                        elif ls[8] == 'NA':
                            bp_plus,bp_minus=[],[]
                            # L breakpoint
                            evals=[ float(i) for i in ls[9].split(';') ]
                            min_eval=min(evals)
                            for eval,tes in zip(ls[9].split(';'), ls[11].split(';;')):
                                if not eval == 'NA':
                                    for te in tes.split(';'):
                                        if not te == '':
                                            ts=te.split(',')
                                            if float(eval) == min_eval:
                                                if ts[3] == '+':
                                                    bp_plus.append(ts)
                                                else:
                                                    bp_minus.append(ts)
                            # normal structure
                            if len(bp_minus) >= 1:
                                bp_str_l=[]
                                d={}
                                for te,_,e,_ in bp_minus:
                                    if not te in d:
                                        d[te]=[]
                                    d[te].append(int(e))
                                for te in d:
                                    breapoint=round(mean(d[te]))
                                    bp_str_l.append('%s,%d,-' % (te, breapoint))
                                pred_status='PASS'
                                pred_res='MEI_left_breakpoint=pT' +';'+ 'MEI_right_breakpoint='  + '|'.join(bp_str_l)
                            # complex structure
                            else:
                                bp_str_l=[]
                                d={}
                                for te,s,_,_ in bp_plus:
                                    if not te in d:
                                        d[te]=[]
                                    d[te].append(int(s))
                                for te in d:
                                    breapoint=round(mean(d[te]))
                                    bp_str_l.append('%s,%d,+' % (te, breapoint))
                                pred_status='complex_structure'
                                pred_res='MEI_left_breakpoint=pT' +';'+ 'MEI_right_breakpoint='  + '|'.join(bp_str_l)

                        l_pos,l_num=ls[3].split(':')
                        r_pos,r_num=ls[4].split(':')
                        l_chim= int(l_num) - int(ls[12])
                        r_chim= int(r_num) - int(ls[13])
                        uniq='yes' if ls[15] == 'singleton' else 'no,%s' % ls[15]
                        tmp=[ls[0], ls[1], ls[2], ls[7], 'MEI_left:ref_pos=%s,chimeric=%d,hybrid=%s' % (l_pos, l_chim, ls[12]), 'MEI_right:ref_pos=%s,chimeric=%d,hybrid=%s' % (r_pos, r_chim, ls[13]), 'confidence:%s' % ls[14], 'unique:%s' % uniq, 'subfamily_pred:status=%s,%s' % (pred_status, pred_res)]
                        tmp= [ str(i) for i in tmp ]
                        outfile.write('\t'.join(tmp) +'\n')
                outfile.flush()
                os.fdatasync(outfile.fileno())

        # main
        if os.path.exists(filenames.bp_merged_groupg) is True:
            convert_to_bed(filenames.bp_merged_groupg, filenames.bp_final_g)
        if os.path.exists(filenames.bp_merged_groupp) is True:
            convert_to_bed(filenames.bp_merged_groupp, filenames.bp_final_p)
        if os.path.exists(filenames.bp_merged_groupf) is True:
            convert_to_bed(filenames.bp_merged_groupf, filenames.bp_final_f)

    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

