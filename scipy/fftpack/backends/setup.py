#!/usr/bin/env python
# Created by Pearu Peterson, August 2002

from os.path import join

def build_backends(config):
    from numpy.distutils.system_info import get_info
    # Build backends for fftpack and convolve
    backends_src = {}
    backends_src['djbfft'] = [join('djbfft/', i) for i in 
                              ['zfft.cxx', 'drfft.cxx', 'convolve.cxx']]
    backends_src['fftw3'] = [join('fftw3/', i) for i in 
                             ['zfft.cxx', 'drfft.cxx', 'zfftnd.cxx']]
    backends_src['fftw2'] = [join('fftw/', i) for i in 
                             ['zfft.cxx', 'drfft.cxx', 'zfftnd.cxx', 'convolve.cxx']]
    backends_src['fftpack'] = [join('fftpack/', i) for i in 
                             ['zfft.cxx', 'drfft.cxx', 'zfftnd.cxx', 'convolve.cxx']]
    backends_src['mkl'] = [join('mkl/', i) for i in 
                             ['zfft.cxx', 'zfftnd.cxx']]

    backends = ['mkl', 'djbfft', 'fftw3', 'fftw2']
    for backend in backends:
        info = get_info(backend)
        if info:
            config.add_library("%s_backend" % backend,
                               sources = backends_src[backend],
                               include_dirs = ["include",
                                               info['include_dirs']])

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('backends', parent_package, top_path)
    build_backends(config)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
