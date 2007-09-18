#!/usr/bin/env python
from os.path import join
import sys

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info

    config = Configuration('nifti',parent_package,top_path)
    #config.add_data_dir('tests')

    include_dirs = ['.']

    # Libraries
    config.add_library('fslio', sources=['fslio.c'], include_dirs=include_dirs)
    config.add_library('niftiio', sources=['nifti1_io.c'], include_dirs=include_dirs)
    config.add_library('znz', sources=['znzlib.c'], include_dirs=include_dirs)

    # Extension
    config.add_extension('_nifticlib',
      sources = ['nifticlib_wrap.c'],
      include_dirs=include_dirs,
      libraries = ['niftiio', 'fslio', 'znz',])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
