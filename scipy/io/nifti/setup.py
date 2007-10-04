#!/usr/bin/env python
import os
from os.path import isfile, join, dirname
import sys
import numpy

nifti_wrapper_file = join('nifticlib.py')

# create an empty file to workaround crappy swig wrapper installation
if not isfile(nifti_wrapper_file):
    open(nifti_wrapper_file, 'w')

# find numpy headers
numpy_headers = join(dirname(numpy.__file__),'core','include')

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info

    config = Configuration('nifti',parent_package,top_path)
    config.add_data_dir('tests')

    include_dirs = [
      '.',
      './nifticlib/fsliolib',
      './nifticlib/niftilib',
      './nifticlib/znzlib']

    nifticlib_headers = ' -I'.join(include_dirs)
    swig_opts = ['-I'+nifticlib_headers, '-I'+numpy_headers]

    # Libraries
    config.add_library('fslio',
      sources=['./nifticlib/fsliolib/fslio.c'], include_dirs=include_dirs)
    config.add_library('niftiio',
      sources=['./nifticlib/niftilib/nifti1_io.c'], include_dirs=include_dirs)
    config.add_library('znz',
      sources=['./nifticlib/znzlib/znzlib.c'], include_dirs=include_dirs)

    # Extension
    config.add_extension('_nifticlib',
      sources = ['nifticlib.i'],
      include_dirs = include_dirs,
      libraries = ['niftiio', 'fslio', 'znz',],
      swig_opts = swig_opts)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
