#!/usr/bin/env python

import os
from os.path import join
import sys
from glob import glob
import shutil

from scipy.distutils.misc_util import Configuration, get_path, dot_join


def configuration(parent_package='',parent_path=None):
    config = Configuration('special', parent_package, parent_path)
    
    define_macros = []
#    if sys.platform=='win32':
#        define_macros.append(('NOINFINITIES',None))
#        define_macros.append(('NONANS',None))


    # C libraries
    config.add_library('c_misc',sources=[join('c_misc','*.c')])
    config.add_library('cephes',sources=[join('cephes','*.c')],
                       macros=define_macros)

    # Fortran libraries
    config.add_library('mach',sources=[join('mach','*.f')])
    config.add_library('toms',sources=[join('amos','*.f')])
    config.add_library('amos',sources=[join('toms','*.f')])
    config.add_library('cdf',sources=[join('cdflib','*.f')])
    config.add_library('specfun',sources=[join('specfun','*.f')])

    # Extension _cephes
    sources = ['_cephesmodule.c', 'amos_wrappers.c', 'specfun_wrappers.c',
               'toms_wrappers.c','cdf_wrappers.c','ufunc_extras.c']
    config.add_extension('_cephes', sources=sources,
                         libraries=['amos','toms','c_misc','cephes','mach',
                                    'cdf', 'specfun'],
                         define_macros = define_macros
                         )
                         
    # Extension nc_specfun

    config.add_extension('specfun',
                         sources=['specfun.pyf'],
                         f2py_options=['--no-wrap-functions'],
                         define_macros=[],
                         libraries=['specfun'])
    return config

if __name__ == '__main__':
    from scipy.distutils.core import setup    
    setup(**configuration(parent_path=''))    
