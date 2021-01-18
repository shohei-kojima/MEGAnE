#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,gzip,shutil
import log,traceback


class empclass:
    pass


def gzip_or_del(args, params, file):
    log.logger.debug('started,file=%s' % file)
    try:
        if args.keep is True:
            with open(file, 'rt') as f_in:
                with gzip.open(file +'.gz', 'wt', compresslevel=params.gzip_compresslevel) as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    f_out.flush()
                    os.fdatasync(f_out.fileno())
        os.remove(file)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def gzip_file(params, file):
    log.logger.debug('started,file=%s' % file)
    try:
        with open(file, 'rt') as f_in:
            with gzip.open(file +'.gz', 'wt', compresslevel=params.gzip_compresslevel) as f_out:
                shutil.copyfileobj(f_in, f_out)
                f_out.flush()
                os.fdatasync(f_out.fileno())
        os.remove(file)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def parse_fasta(path_to_file):
    tmp={}
    seq=''
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line and seq:
                tmp[header]=seq
                header=line.strip().replace(' ', '_')
                seq=''
            elif '>' in line and not seq:
                header=line.strip().replace(' ', '_')
            else:
                seq += line.strip()
        tmp[header]=seq
    return tmp


def parse_fasta_transd(path_to_file):
    tmp=[]
    seq=''
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line and seq:
                tmp.append(seq)
                seq=''
            elif '>' in line and not seq:
                pass
            else:
                seq += line.strip()
        tmp.append(seq)
    return tmp


def parse_fasta_for_merge_vcf(path_to_file):
    tmp={}
    seq=[]
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line and seq:
                tmp[header]=''.join(seq)
                header=line.strip()
                seq=[]
            elif '>' in line and not seq:
                header=line.strip()
            else:
                seq.append(line.strip())
        tmp[header]=''.join(seq)
    return tmp


def load_me_classification(path_to_file):
    mes,all_clas={},set()
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line:
                ls=line.strip().replace('>', '').split('\t')
                if len(ls) >= 2:
                    clas=ls[1]
                    if ls[0] == 'SVA2':  # exist in humrep.ref of RepBase24.01
                        clas='SINE'
                clas=clas.replace(' ', '_')
                mes[ls[0]]=clas
                all_clas.add(clas)
    return mes, all_clas


def gzip_d(gzfile):
    log.logger.debug('started,infile=%s' % gzfile)
    try:
        outfpath=gzfile[:-3]
        with open(outfpath, 'w') as outfile:
            with gzip.open(gzfile) as infile:
                for line in infile:
                    line=line.decode()
                    outfile.write(line)
            outfile.flush()
            os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def output_finish_comment(do_ins, do_abs, filenames):
    if do_ins is True and do_abs is True:
        if os.path.exists(filenames.bp_final_g) is True and os.path.exists(filenames.bp_final_p) is True:
            log.logger.info('\n\n\033[34mAll analysis finished! Thank you for using MEGAnE!\n\n  When you next generate joint VCF, please use:  MEI_final_gaussian_genotyped.vcf & absent_MEs_genotyped.vcf\n  When you analyze only this sample, please use: MEI_final_percentile_genotyped.vcf & absent_MEs_genotyped.vcf\n\033[0m\n')
        elif os.path.exists(filenames.bp_final_f) is True:
            log.logger.info('\n\n\033[34mAll analysis finished! Thank you for using MEGAnE!\n\n  We failed to call and genotype pMEIs. Please check again whether the input BAM/CRAM contains abundant pMEI.\n\033[0m\n')
    elif do_ins is True:
        if os.path.exists(filenames.bp_final_g) is True and os.path.exists(filenames.bp_final_p) is True:
            log.logger.info('\n\n\033[34mAll analysis finished! Thank you for using MEGAnE!\n\n  When you next generate joint VCF, please use:  MEI_final_gaussian_genotyped.vcf\n  When you analyze only this sample, please use: MEI_final_percentile_genotyped.vcf\n\033[0m\n')
        elif os.path.exists(filenames.bp_final_f) is True:
            log.logger.info('\n\n\033[34mAll analysis finished! Thank you for using MEGAnE!\n\n  We failed to call and genotype pMEIs. Please check again whether the input BAM/CRAM contains abundant pMEI.\n\033[0m\n')
    elif do_abs is True:
        log.logger.info('\n\n\033[34mAll analysis finished!\n\n  Please use following file for your analysis:  absent_MEs_genotyped.vcf\n\033[0m\n')


def output_finish_comment_merge_vcf(args, filenames):
    if args.merge_mei is True:
        log.logger.info('\n\n\033[34mAll analysis finished! Thank you for using MEGAnE!\n\n  Please use following file for your analysis:  %s\n\033[0m\n' % filenames.filled_vcf_ins)
    elif args.merge_absent_me is True:
        log.logger.info('\n\n\033[34mAll analysis finished! Thank you for using MEGAnE!\n\n  Please use following file for your analysis:  %s\n\033[0m\n' % filenames.filled_vcf_abs)

