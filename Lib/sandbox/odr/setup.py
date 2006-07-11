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

    atlas_info = get_info('atlas')
    #atlas_info = {} # uncomment if ATLAS is available but want to use
                     # Fortran LAPACK/BLAS; useful for testing
    blas_libs = []
    if not atlas_info:
        warnings.warn(AtlasNotFoundError.__doc__)
        blas_info = get_info('blas')
        if blas_info:
            libodr_files.append('d_lpk.f')
            blas_libs.extend(blas_info['libraries'])
        else:
            warnings.warn(BlasNotFoundError.__doc__)
            libodr_files.append('d_lpkbls.f')
    else:
        libodr_files.append('d_lpk.f')
        blas_libs.extend(atlas_info['libraries'])

    libodr = [os.path.join('odrpack', x) for x in libodr_files]
    config.add_library('odrpack', sources=libodr)
    sources = ['__odrpack.c']
    config.add_extension('__odrpack',
                         sources=sources,
                         libraries=['odrpack']+blas_libs,
                         include_dirs=['.'],
                         library_dirs=atlas_info['library_dirs'],
                         )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
