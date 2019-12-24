#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import sys
from os.path import abspath,dirname,realpath,join

def init():
    global base
    base=abspath(dirname(realpath(__file__)))
    sys.path.append(join(base, 'src'))

