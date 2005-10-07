#!/usr/bin/env python

from __future__ import nested_scopes
import os
from os.path import join
import glob

from scipy.distutils.misc_util import Configuration
from scipy.distutils.system_info import get_info

def configuration(parent_package='',parent_path=None):
    config = Configuration('integrate', parent_package, parent_path)

    blas_opt = get_info('blas_opt')
    if not blas_opt:
        raise NotFoundError,'no blas resources found'

    config.add_library('linpack_lite',
                       sources=[join('linpack_lite','*.f')])
    config.add_library('mach',
                       sources=[join('mach','*.f')])
    config.add_library('quadpack',
                       sources=[join('quadpack','*.f')])
    config.add_library('odepack',
                       sources=[join('odepack','*.f')])
    # should we try to weed through files and replace with calls to
    # LAPACK routines?
    # Yes, someday...

    

    # Extensions
    # quadpack:

    config.add_extension('_quadpack',
                         sources=['_quadpackmodule.c'],
                         libraries=['quadpack', 'linpack_lite', 'mach'])
    # odepack
    config.add_extension('_odepack',
                         sources=['_odepackmodule.c'],
                         libraries=['odepack','linpack_lite','mach'],
                         **blas_opt)
    
    # vode
    config.add_extension('vode',
                         sources=['vode.pyf'],
                         libraries=['odepack','linpack_lite','mach'])

    return config

if __name__ == '__main__':    
    from scipy.distutils.core import setup
    setup(**configuration(parent_path=''))
