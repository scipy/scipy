#!/usr/bin/env python
# Created by Pearu Peterson, August 2002

from os.path import join

def build_backends(config):
    from numpy.distutils.system_info import get_info
    src = [join('src', i) for i in 
                             ['zfft.cxx', 'drfft.cxx', 'zfftnd.cxx']]
    info = get_info("fftw3")
    if info:
        config.add_library("fftw3_backend",
                           sources = src,
                           include_dirs = ["../common",
                                           info['include_dirs']])
        config.add_extension("_fftw3", sources = ["fftw3.pyf"], extra_info = info, libraries = ["fftw3_backend"])
    
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('fftw3', parent_package, top_path)
    build_backends(config)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
