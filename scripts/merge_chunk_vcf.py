#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020-2021 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,gzip
import log,traceback

def merge(args, params, filenames):
    log.logger.debug('started')
    try:
        # path check
        if os.path.exists(args.chunk_vcf_list) is False:
            log.logger.error('%s does not exist. Please check the path again.' % args.chunk_vcf_list)
            exit(1)
        paths=[]
        with open(args.chunk_vcf_list) as infile:
            for line in infile:
                path=line.strip()
                if os.path.exists(path) is False:
                    log.logger.error('%s does not exist. Please check the path again.' % path)
                    exit(1)
                if not '.vcf.gz' in path:
                    log.logger.error('%s may not be gziped VCF. Input files must be end with ".vcf.gz" extention.' % path)
                    exit(1)
                paths.append(path)
        log.logger.info('%d chunk vcf.gz files found.' % len(paths))
        
        # outfile obj
        if args.merge_mei is True:
            outfname=filenames.filled_vcf_ins
        else:
            outfname=filenames.filled_vcf_abs
        outfile=gzip.open(outfname, 'wt')
        
        # header
        line_n=0
        with gzip.open(paths[0], 'rt') as infile:
            for line in infile:
                line_n += 1
                if line[:6] == '#CHROM':
                    ls=line.split('\t', 10)
                    hs=ls[:9]
                    break
                outfile.write(line)
        
        # read all chunks
        infiles=[]
        for path in paths:
            infile=gzip.open(path, 'rt')
            for line in infile:
                if line[:6] == '#CHROM':
                    ls=line.strip().split()
                    hs.extend(ls[9:])
                    break
            infiles.append(infile)
        
        outfile.write('\t'.join(hs) +'\n')
        sample_ids_n= len(hs) - 9
        
        var_n=0
        for line in infiles[0]:
            line_n += 1
            tmp=[]
            tmp.append(line.strip().split('\t'))
            for infile in infiles[1:]:
                line=next(infile)
                tmp.append(line.strip().split('\t'))
            var=tmp[0][2]
            for ls in tmp[1:]:
                if not ls[2] == var:
                    log.logger.error('Input files have different variants at line %d. Please check the input files are generated with the same methods.' % line_n)
                    exit(1)
            new_ls=tmp[0][:9]
            for ls in tmp:
                new_ls.extend(ls[9:])
            
            ac=0
            missing=0
            zero=0
            for v in new_ls[9:]:
                if '.' in v:
                    missing += 1
                elif v == '1/1':
                    ac += 2
                elif v == '0/0':
                    zero += 1
                else:
                    ac += 1
            
            if args.merge_mei is True:
                change_to_nonpass=False
                if (ac / (sample_ids_n * 2)) < 0.05:
                    if missing >= (ac * 2):
                        change_to_nonpass=True
                elif missing >= ac or missing > (zero / 2):
                    change_to_nonpass=True
                if change_to_nonpass is True:
                    if new_ls[6] == 'PASS':
                        new_ls[6]='M'
                    else:
                        new_ls[6]='%s;M' % new_ls[6]
            if args.merge_absent_me is True:
                change_to_nonpass=False
                if (ac / (sample_ids_n * 2)) < 0.05:
                    if missing >= (ac * 2):
                        change_to_nonpass=True
                elif ((ac + missing) / (sample_ids_n * 2)) < 0.75:
                    if missing >= ac or missing > zero:
                        change_to_nonpass=True
                if change_to_nonpass is True:
                    if new_ls[6] == 'PASS':
                        new_ls[6]='M'
                    else:
                        new_ls[6]='%s;M' % new_ls[6]
            
            outfile.write('\t'.join(new_ls) +'\n')
            var_n += 1
            if (var_n % 200) == 0:
                log.logger.info('%d variants processed.' % var_n)
        outfile.close()
        log.logger.info('%d variants processed.' % var_n)
        for infile in infiles:
            infile.close()
    
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
