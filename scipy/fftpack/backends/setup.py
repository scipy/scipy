#!/usr/bin/env python
# Created by Pearu Peterson, August 2002

from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('backends', parent_package, top_path)
        
    config.add_subpackage("fftw3")
    config.add_subpackage("fftw")
    config.add_subpackage("mkl")

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
