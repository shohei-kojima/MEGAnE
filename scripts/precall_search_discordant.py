#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os
import pysam
import pybedtools
from pybedtools import BedTool
import log,traceback


cigar_op={'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}
cigar_ref_retain={'M', 'D', 'N', '=', 'X'}
cigar_read_retain={'M', 'I', '=', 'X'}


def detect_discordant(args, params, filenames):
    log.logger.debug('started')
    try:
        pybedtools.set_tempdir(args.pybedtools_tmp)
        
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
        
        if not args.b is None:
            infile=pysam.AlignmentFile(filenames.limited_b, 'rb')
        else:
            infile=pysam.AlignmentFile(filenames.limited_c, 'rc', reference_filename=args.fa)
        disc_bed=[]
        for line in infile:
            line=line.tostring()
            ls=line.strip().split('\t')
            if int(ls[1]) < 2048:  # remove supplementary alignment
                b=bin(int(ls[1]))
                if b[-1] == '1':   # paired-end
                    if 'S' in ls[5] or 'H' in ls[5]:
                        ref_len=calc_ref_len(ls[5])
                        start= int(ls[3]) - 1   # 0-based
                        end= start + ref_len
                        disc_bed.append('%s\t%d\t%d\n' % (ls[2], start, end))
        disc_bed=BedTool(''.join(disc_bed), from_string=True)
        
        def generate_slopbed(args, params, filenames, ins_slop_len, ins_discard_flank_len, abs_slop_len):
            do_ins=False if args.only_abs is True else True
            do_abs=False if args.only_ins is True else True
            slopbed=[]
            count_ins,count_abs=0,0
            if do_ins is True:
                ins_bed_found=False
                iids=set()
                for bp_final in [filenames.bp_final_g, filenames.bp_final_p, filenames.bp_final_f, filenames.bp_final_u, filenames.bp_final_d]:
                    if os.path.exists(bp_final) is True:
                        ins_bed_found=True
                        with open(bp_final) as infile:
                            log.logger.debug('%s loading.' % bp_final)
                            for line in infile:
                                ls=line.split()
                                slopbed.append('%s\t%s\t%s\t%s\n' % (ls[0], int(ls[1]) - ins_slop_len, int(ls[1]) - ins_discard_flank_len, 'l%s' % ls[10]))
                                slopbed.append('%s\t%s\t%s\t%s\n' % (ls[0], int(ls[2]) + ins_discard_flank_len, int(ls[2]) + ins_slop_len, 'r%s' % ls[10]))
                                iids.add(ls[10])
                                count_ins += 1
                if ins_bed_found is False:
                    log.logger.error('Available ins_bed not found.')
                    exit(1)
            if do_abs is True:
                abs_bed_found=False
                tmpl,tmpr=[],[]
                id=0
                aids=set()
                for abs_bed_f in [args.abs_bed, args.abs_3t_bed]:
                    if not abs_bed_f is None:
                        if os.path.exists(abs_bed_f) is True:
                            abs_bed_found=True
                            with open(abs_bed_f) as infile:
                                log.logger.debug('%s loading.' % abs_bed_f)
                                for line in infile:
                                    ls=line.split()
                                    tmpl.append('%s\t%s\t%s\t%s\n' % (ls[0], ls[1], ls[1], 'laID=%d' % id))
                                    tmpr.append('%s\t%s\t%s\t%s\n' % (ls[0], ls[2], ls[2], 'raID=%d' % id))
                                    aids.add('aID=%d' % id)
                                    count_abs += 1
                                    id += 1
                absbedl=BedTool(''.join(tmpl), from_string=True).slop(l=abs_slop_len, r=0, g=args.fai)
                absbedr=BedTool(''.join(tmpr), from_string=True).slop(l=0, r=abs_slop_len, g=args.fai)
                if abs_bed_found is False:
                    log.logger.error('Available abs_bed not found.')
                    exit(1)
                for line in absbedl:
                    slopbed.append(str(line).strip())
                for line in absbedr:
                    slopbed.append(str(line).strip())
            slopbed=BedTool('\n'.join(slopbed), from_string=True)
            return slopbed,count_ins,count_abs,iids,aids
        
        slopbed,count_ins,count_abs,iids,aids=generate_slopbed(args, params, filenames, params.ins_slop_len_for_disc_detection, params.ins_discard_flank_len, params.abs_slop_len_for_disc_detection)
        log.logger.info('total ins=%d,abs=%d bed lines found.' % (count_ins, count_abs))
        
        intersect=slopbed.intersect(disc_bed, wa=True, nonamecheck=True)
        disc_ids_d={}
        for id in iids:
            disc_ids_d['l%s' % id]=0
            disc_ids_d['r%s' % id]=0
        for id in aids:
            disc_ids_d['l%s' % id]=0
            disc_ids_d['r%s' % id]=0
        for line in intersect:
            ls=str(line).split()
            disc_ids_d[ls[3]] += 1
        count_ins,count_abs=0,0
        global disc_ids
        disc_ids={}
        for id in aids:
            if disc_ids_d['l%s' % id] >= params.abs_disc_ids_threshold and disc_ids_d['r%s' % id] >= params.abs_disc_ids_threshold:
                disc_ids[id]='both'
                count_abs += 1
            elif disc_ids_d['l%s' % id] >= params.abs_disc_ids_threshold:
                disc_ids[id]='half'
                count_abs += 1
            elif disc_ids_d['r%s' % id] >= params.abs_disc_ids_threshold:
                disc_ids[id]='half'
                count_abs += 1
        ins_disc_ids_threshold= args.cov * params.ins_disc_ids_threshold_coeff
        for id in iids:
            if disc_ids_d['l%s' % id] >= ins_disc_ids_threshold and disc_ids_d['r%s' % id] >= ins_disc_ids_threshold:
                disc_ids[id]='both'
                count_ins += 1
            elif disc_ids_d['l%s' % id] >= ins_disc_ids_threshold:
                disc_ids[id]='half'
                count_ins += 1
            elif disc_ids_d['r%s' % id] >= ins_disc_ids_threshold:
                disc_ids[id]='half'
                count_ins += 1
        log.logger.info('total ins=%d,abs=%d bed lines with discordant read(s) found.' % (count_ins, count_abs))
        
        pybedtools.cleanup()
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

