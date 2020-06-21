#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,gzip,subprocess
import pysam
from pybedtools import BedTool
import log,traceback


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
def genotype_ins(args, params, filenames):
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
            
        chim_count={}
        back_to_tsd_ratios=[]
        back_to_chim_ratios=[]
        structures=[]
        chimerics=[]
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
                            structure='no_tsd'
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
                    structure='no_read'
                left_num=int(ls[4].split(',')[1].split('=')[1])
                right_num=int(ls[5].split(',')[1].split('=')[1])
                if left_depth + right_depth > 0:
                    back_to_tsd_ratio= tsd_depth / ((left_depth + right_depth) / 2)
                    back_to_chim_ratio= (left_num + right_num) / ((left_depth + right_depth) / 2)
                else:
                    back_to_tsd_ratio=-1
                    back_to_chim_ratio=-1
                back_to_tsd_ratios.append(back_to_tsd_ratio)
                back_to_chim_ratios.append(back_to_chim_ratio)
                chimerics.append(left_num + right_num)
                structures.append(structure)
        out=[]
        for r,c,ch,s in zip(back_to_tsd_ratios, back_to_chim_ratios, chimerics, structures):
            out.append('%f\t%f\t%d\t%s\n' % (r, c, ch, s))
        with open('./genotype_out/tmp_back_to_tsd_ratio.txt', 'w') as outfile:
            outfile.write(''.join(out))
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

