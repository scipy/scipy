#!/usr/bin/env python
import numpy

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('convolve',parent_package,top_path,
                           package_path='lib')

    config.add_extension('_correlate',
                         sources=["src/_correlatemodule.c"],
                         define_macros = [('NUMPY', '1')],
                         include_dirs = [numpy.get_numarray_include()])
    config.add_extension('_lineshape',
                         sources=["src/_lineshapemodule.c"],
                         define_macros = [('NUMPY', '1')],
                         include_dirs = [numpy.get_numarray_include()])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    config = configuration(top_path='').todict()
    setup(author='Todd Miller',
          author_email = 'help@stsci.edu',
          description = 'image array convolution functions',
          version = '2.0',
          **config)
