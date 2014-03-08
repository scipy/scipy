#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import os
import sys
import glob
import subprocess

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('sparse',parent_package,top_path)

    config.add_data_dir('tests')
    config.add_data_dir('benchmarks')

    config.add_subpackage('linalg')
    config.add_subpackage('sparsetools')
    config.add_subpackage('csgraph')

    config.add_extension('_csparsetools',
                         sources=['_csparsetools.c'])

    def get_sparsetools_sources(ext, build_dir):
        # Defer generation of source files
        subprocess.call([sys.executable,
                         os.path.join(os.path.dirname(__file__),
                                      'generate_sparsetools.py'),
                         '--no-force'])
        return []

    depends = [os.path.join('sparsetools', 'sparsetools_gen.h')]
    depends += glob.glob('sparsetools/*.h')
    config.add_extension('sparsetools',
                         define_macros=[('__STDC_FORMAT_MACROS', 1)],
                         depends=depends,
                         include_dirs=['sparsetools'],
                         sources=[os.path.join('sparsetools', 'sparsetools.cxx'),
                                  get_sparsetools_sources]
                         )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
