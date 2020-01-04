#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,logging

def start_log(args):
    global logger
    logger=logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    ch=logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter=logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    ch=logging.FileHandler(os.path.join(args.outdir, args.logfilename))
    formatter=logging.Formatter('%(asctime)s:%(levelname)s:file=%(filename)s:module=%(module)s:funcName=%(funcName)s:line=%(lineno)d:message=%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


