#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension('extract_discordant_c', ['extract_discordant_c.pyx'])]   #assign new.pyx module in setup.py.
setup(
      name        = 'extract_discordant_c app',
      cmdclass    = {'build_ext':build_ext},
      ext_modules = ext_modules
      )

'''
cat extract_discordant.py > extract_discordant_c.pyx
python cython_setup.py build_ext --inplace
'''
