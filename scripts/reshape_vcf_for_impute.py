#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,gzip,argparse,collections
import log,traceback


def resolve_overlap_in_a_vcf(f, outfname_base):
    log.logger.debug('started')
    try:
        out_nonoverlap='%s_biallelic.vcf.gz' % outfname_base
        
        # resolve overlap
        if f[-7:] == '.vcf.gz' or f[-8:] == '.vcf.bgz':
            infile=gzip.open(f)
            gzip_judge=True
        elif f[-4:] == '.vcf':
            infile=open(f)
            gzip_judge=False
        else:
            log.logger.error('please specify .vcf file')
            exit(1)
        
        d={}
        filts={}
        for line in infile:
            if gzip_judge is True:
                line=line.decode()
            if not line[0] == '#':
                ls=line.split('\t', 10)
                if not ls[0] in d:
                    d[ls[0]]=collections.Counter()
                d[ls[0]][int(ls[1])] += 1
                filts[ls[2]]=ls[6]

        hits={}
        overlap_ins_ins=0
        for chr in d:
            for pos in d[chr]:
                if d[chr][pos] >= 2:
                    if not chr in hits:
                        hits[chr]=set()
                    hits[chr].add(str(pos))
                    overlap_ins_ins += d[chr][pos]

        if overlap_ins_ins >= 1:
            log.logger.info('Multi-allelic ME in %s detected, n = %d' % (f, overlap_ins_ins))
            
            infile.seek(0)
            hits_lines={}
            for chr in hits:
                hits_lines[chr]={}
                for pos in hits[chr]:
                    hits_lines[chr][pos]=[]
            for line in infile:
                if gzip_judge is True:
                    line=line.decode()
                if not line[0] == '#':
                    ls=line.split('\t', 10)
                    if ls[0] in hits:
                        if ls[1] in hits[ls[0]]:
                            hits_lines[ls[0]][ls[1]].append(line)
            
            keep={}
            for chr in hits:
                keep[chr]={}
                for pos in hits[chr]:
                    keep[chr][pos]=[]
            
            for chr in hits_lines:
                for pos in hits_lines[chr]:
                    acs=[]
                    for line in hits_lines[chr][pos]:
                        ls=line.split('\t', 10)
                        for info in ls[7].split(';')[::-1]:
                            if 'AC=' in info:
                                ac=int(info.replace('AC=', ''))
                                break
                        filt=filts[ls[2]]
                        filt= 1 if filt == 'PASS' else 0
                        acs.append([filt, ac, line, ls[2]])
                    acs=sorted(acs)
                    acs=[ l[1:] for l in acs ]
                    keep[chr][pos]=acs[-1][1]
                    outchar=[]
                    outchar.append('keep=%s,%d' % (acs[-1][2], acs[-1][0]))
                    for ac,_,id in acs[:-1]:
                        outchar.append('discard=%s,%d' % (id, ac))
                    outchar=';'.join(outchar)
                    log.logger.debug(outchar)
                    
            infile.seek(0)
            outfile=gzip.open(out_nonoverlap, 'wt')
            for line in infile:
                if gzip_judge is True:
                    line=line.decode()
                if not line[0] == '#':
                    ls=line.split('\t', 10)
                    if ls[0] in keep:
                        if ls[1] in keep[ls[0]]:
                            if len(keep[ls[0]][ls[1]]) > 0:
                                outfile.write(keep[ls[0]][ls[1]])
                                keep[ls[0]][ls[1]]=''
                        else:
                            outfile.write(line)
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)
            infile.close()
            outfile.close()
        else:
            log.logger.info('Multi-allelic ME in %s was not detected.' % f)
            
            infile.seek(0)
            outfile=gzip.open(out_nonoverlap, 'wt')
            for line in infile:
                if gzip_judge is True:
                    line=line.decode()
                if not line[0] == '#':
                    ls=line.split('\t', 10)
                    outfile.write(line)
                else:
                    outfile.write(line)
            infile.close()
            outfile.close()
        return out_nonoverlap
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def resolve_overlap_between_vcf_vcf(ins, abs, cohort):
    log.logger.debug('started')
    try:
        out_nonoverlap='%s_biallelic.vcf.gz' % cohort
        
        # resolve overlap
        f=ins
        if f[-7:] == '.vcf.gz' or f[-8:] == '.vcf.bgz':
            insfile=gzip.open(f)
            gzip_judge=True
        else:
            insfile=open(f)
            gzip_judge=False
        
        f=abs
        if f[-7:] == '.vcf.gz' or f[-8:] == '.vcf.bgz':
            absfile=gzip.open(f)
            gzip_judge=True
        else:
            absfile=open(f)
            gzip_judge=False
        
        d={}
        filts={}
        for line in insfile:
            if gzip_judge is True:
                line=line.decode()
            if not line[0] == '#':
                ls=line.split('\t', 10)
                if not ls[0] in d:
                    d[ls[0]]=collections.Counter()
                d[ls[0]][int(ls[1])] += 1
                filts[ls[2]]=ls[6]
        for line in absfile:
            if gzip_judge is True:
                line=line.decode()
            if not line[0] == '#':
                ls=line.split('\t', 10)
                if not ls[0] in d:
                    d[ls[0]]=collections.Counter()
                d[ls[0]][int(ls[1])] += 1
                filts[ls[2]]=ls[6]
        
        hits={}
        overlap_ins_abs=0
        for chr in d:
            for pos in d[chr]:
                if d[chr][pos] >= 2:
                    if not chr in hits:
                        hits[chr]=set()
                    hits[chr].add(str(pos))
                    overlap_ins_abs += d[chr][pos]
        
        if overlap_ins_abs >= 1:
            log.logger.info('Multi-allelic ME between %s and %s detected, n = %d' % (ins, abs, overlap_ins_abs))
            
            insfile.seek(0)
            absfile.seek(0)
            hits_lines={}
            for chr in hits:
                hits_lines[chr]={}
                for pos in hits[chr]:
                    hits_lines[chr][pos]=[]
            for line in insfile:
                if gzip_judge is True:
                    line=line.decode()
                if not line[0] == '#':
                    ls=line.split('\t', 10)
                    if ls[0] in hits:
                        if ls[1] in hits[ls[0]]:
                            hits_lines[ls[0]][ls[1]].append(line)
            for line in absfile:
                if gzip_judge is True:
                    line=line.decode()
                if not line[0] == '#':
                    ls=line.split('\t', 10)
                    if ls[0] in hits:
                        if ls[1] in hits[ls[0]]:
                            hits_lines[ls[0]][ls[1]].append(line)
            
            keep={}
            for chr in hits:
                keep[chr]={}
                for pos in hits[chr]:
                    keep[chr][pos]=[]
            
            for chr in hits_lines:
                for pos in hits_lines[chr]:
                    acs=[]
                    for line in hits_lines[chr][pos]:
                        ls=line.split('\t', 10)
                        for info in ls[7].split(';')[::-1]:
                            if 'AC=' in info:
                                ac=int(info.replace('AC=', ''))
                                break
                        filt=filts[ls[2]]
                        filt= 1 if filt == 'PASS' else 0
                        acs.append([filt, ac, line, ls[2]])
                    acs=sorted(acs)
                    acs=[ l[1:] for l in acs ]
                    keep[chr][pos]=acs[-1][1]
                    outchar=[]
                    outchar.append('keep=%s,%d' % (acs[-1][2], acs[-1][0]))
                    for ac,_,id in acs[:-1]:
                        outchar.append('discard=%s,%d' % (id, ac))
                    outchar=';'.join(outchar)
                    log.logger.debug(outchar)
                    
            insfile.seek(0)
            absfile.seek(0)
            outfile=gzip.open(out_nonoverlap, 'wt')
            for line in insfile:
                if gzip_judge is True:
                    line=line.decode()
                if not line[0] == '#':
                    ls=line.split('\t', 10)
                    if ls[0] in keep:
                        if ls[1] in keep[ls[0]]:
                            if ls[2] in keep[ls[0]][ls[1]]:
                                outfile.write(keep[ls[0]][ls[1]])
                        else:
                            outfile.write(line)
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)
            for line in absfile:
                if gzip_judge is True:
                    line=line.decode()
                if not line[0] == '#':
                    ls=line.split('\t', 10)
                    if ls[0] in keep:
                        if ls[1] in keep[ls[0]]:
                            if ls[2] in keep[ls[0]][ls[1]]:
                                outfile.write(keep[ls[0]][ls[1]])
                        else:
                            outfile.write(line)
                    else:
                        outfile.write(line)
            insfile.close()
            absfile.close()
            outfile.close()
        else:
            log.logger.info('Multi-allelic ME between %s and %s was not detected.' % (ins, abs))
            
            insfile.seek(0)
            absfile.seek(0)
            outfile=gzip.open(out_nonoverlap, 'wt')
            for line in insfile:
                if gzip_judge is True:
                    line=line.decode()
                if not line[0] == '#':
                    ls=line.split('\t', 10)
                    outfile.write(line)
                else:
                    outfile.write(line)
            for line in absfile:
                if gzip_judge is True:
                    line=line.decode()
                if not line[0] == '#':
                    ls=line.split('\t', 10)
                    outfile.write(line)
            insfile.close()
            absfile.close()
            outfile.close()
        return out_nonoverlap
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def convert_to_bed(final_outf):
    log.logger.debug('started')
    try:
        outfname= final_outf[:-7] + '.bed.gz'
        infile=gzip.open(final_outf)
        outfile=gzip.open(outfname, 'wt')
        
        for line in infile:
            line=line.decode()
            if not line[0] == '#':
                ls=line.split('\t', 10)
                for info in ls[7].split(';'):
                    if '0END=' in info:
                        end=info.replace('0END=', '')
                        break
                start= int(ls[1]) - 1
                outfile.write('%s\t%d\t%s\t%s\n' % (ls[0], start, end, ls[2]))
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


