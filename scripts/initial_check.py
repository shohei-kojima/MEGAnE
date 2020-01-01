#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,multiprocessing
from os.path import abspath,dirname,realpath,join

# http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    elif 'PATH' in os.environ:
        for path in os.environ['PATH'].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def check(args):
    # check python version
    version=sys.version_info
    if (version[0] >= 3) and (version[1] >= 7):
        pass
#        print('Python version check OK. You are using Python 3.7 or later.')
    else:
        print('Error: Please use Python 3.7 or later. Your Python is version %d.%d.' % (version[0], version[1]))
        exit()
    
    # check cpu num
    cpu_num=multiprocessing.cpu_count()
    if args.p > cpu_num:
        print('Error: Too many thread number. Please specify the number less than your cpu cores. You specified = %d, cpu cores = %d.' % (args.p, cpu_num))
        exit()
    
    # check PATH
    for i in ['blastn', 'bedtools']:
        if which(i) is None:
            print('Warning: %s not found in $PATH. Please check %s is installed and added to PATH. Please ignore when you are using job scheduler. Proceed anyway.' % (i, i))

    # check prerequisite modules
    from Bio.Seq import Seq
    import gzip
    from pybedtools import BedTool
    import matplotlib
    import pysam
    from Bio.Blast.Applications import NcbiblastnCommandline
