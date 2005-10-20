#!/usr/bin/env python
# Created by Pearu Peterson, August 2002

import os
import sys
join = os.path.join

def configuration(parent_package='',top_path=None):
    from scipy.distutils.misc_util import Configuration
    from scipy.distutils.system_info import get_info,FFTWNotFoundError,\
         DJBFFTNotFoundError
    config = Configuration('fftpack',parent_package, top_path)
    fftw_info = get_info('fftw') or get_info('dfftw')
    if not fftw_info:
        print FFTWNotFoundError.__doc__
    djbfft_info = get_info('djbfft')
    if not djbfft_info:
        print DJBFFTNotFoundError.__doc__
    #djbfft_info = None
    #fftw_info = None
    
    config.add_data_dir('tests')

    config.add_library('dfftpack',
                       sources=[join('dfftpack','*.f')])

    sources = ['fftpack.pyf','src/zfft.c','src/drfft.c','src/zrfft.c',
               'src/zfftnd.c']
    extra_info = []
    if fftw_info:
        extra_info.append(fftw_info)
    if djbfft_info:
        extra_info.append(djbfft_info)
    config.add_extension('_fftpack',
                         sources=sources,
                         libraries=['dfftpack'],
                         extra_info = extra_info
                         )

    config.add_extension('convolve',
                         sources = ['convolve.pyf','src/convolve.c'],
                         libraries = ['dfftpack'],
                         extra_info = extra_info
                         )
    return config

if __name__ == '__main__':
    from scipy.distutils.core import setup
    from fftpack_version import fftpack_version
    setup(version=fftpack_version,
          description='fftpack - Discrete Fourier Transform package',
          author='Pearu Peterson',
          author_email = 'pearu@cens.ioc.ee',
          maintainer_email = 'scipy-dev@scipy.org',
          license = 'SciPy License (BSD Style)',
          **configuration(top_path='').todict())
