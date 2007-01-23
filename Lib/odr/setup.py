#!/usr/bin/env python

import os,sys,re
from distutils import dep_util
from glob import glob
import warnings

from numpy.distutils.core import Extension
from numpy.distutils.misc_util import get_path, Configuration, dot_join

from numpy.distutils.system_info import get_info,dict_append,\
     AtlasNotFoundError,LapackNotFoundError,BlasNotFoundError,\
     LapackSrcNotFoundError,BlasSrcNotFoundError

def configuration(parent_package='', top_path=None):
    config = Configuration('odr', parent_package, top_path)

    libodr_files = ['d_odr.f',
                    'd_mprec.f',
                    'dlunoc.f']

    blas_info = get_info('blas_opt')
    if blas_info:
        libodr_files.append('d_lpk.f')
    else:
        warnings.warn(BlasNotFoundError.__doc__)
        libodr_files.append('d_lpkbls.f')

    libodr = [os.path.join('odrpack', x) for x in libodr_files]
    config.add_library('odrpack', sources=libodr)
    sources = ['__odrpack.c']
    libraries = ['odrpack'] + blas_info.pop('libraries', [])
    include_dirs = ['.'] + blas_info.pop('include_dirs', [])
    config.add_extension('__odrpack',
        sources=sources,
        libraries=libraries,
        include_dirs=include_dirs,
        **blas_info
    )

    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
