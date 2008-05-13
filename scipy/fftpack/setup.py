#!/usr/bin/env python
# Created by Pearu Peterson, August 2002

from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('fftpack',parent_package, top_path)

    backends = ['mkl', 'djbfft', 'fftw3', 'fftw', 'fftpack']
    info = dict([(k, False) for k in backends])
    info['fftpack'] = True

    djbfft_info = {}
    mkl_info = get_info('mkl')
    if mkl_info:
        mkl_info.setdefault('define_macros', []).append(('SCIPY_MKL_H', None))
        fft_opt_info = mkl_info
        info['mkl'] = True
    else:
        # Take the first in the list
        for b in ['fftw3', 'fftw']:
            tmp = get_info(b)
            if tmp:
                fft_opt_info = tmp
                info[b] = True
                break
        djbfft_info = get_info('djbfft')
        info['djbfft'] = True

    config.add_data_dir('tests')
    config.add_data_dir('benchmarks')

    config.add_library('dfftpack',
                       sources=[join('dfftpack','*.f')])

    backends_src = {}
    backends_src['djbfft'] = [join('src/djbfft/', i) for i in 
                              ['zfft.cxx', 'drfft.cxx']]
    backends_src['fftw3'] = [join('src/fftw3/', i) for i in 
                             ['zfft.cxx', 'drfft.cxx', 'zfftnd.cxx']]
    backends_src['fftw'] = [join('src/fftw/', i) for i in 
                             ['zfft.cxx', 'drfft.cxx', 'zfftnd.cxx']]
    backends_src['fftpack'] = [join('src/fftpack/', i) for i in 
                             ['zfft.cxx', 'drfft.cxx', 'zfftnd.cxx']]

    libs = ['dfftpack']
    for b in ['djbfft', 'fftw3', 'fftw', 'fftpack']:
        if info[b]:
            config.add_library('%s_backend' % b, 
                    sources = backends_src[b], 
                    include_dirs = ['src', djbfft_info['include_dirs'],
                                    fft_opt_info['include_dirs']])
            libs.append(b)

    sources = ['fftpack.pyf', 'src/fftpack.cxx', 'src/zrfft.c']

    config.add_extension('_fftpack',
        sources=sources,
        libraries = libs,
        extra_info=[fft_opt_info, djbfft_info],
        include_dirs = ['src'],
    )

    config.add_extension('convolve',
        sources = ['convolve.pyf','src/convolve.cxx'],
        libraries = libs,
        extra_info = [fft_opt_info, djbfft_info],
        include_dirs = ['src'],
    )
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    from fftpack_version import fftpack_version
    setup(version=fftpack_version,
          description='fftpack - Discrete Fourier Transform package',
          author='Pearu Peterson',
          author_email = 'pearu@cens.ioc.ee',
          maintainer_email = 'scipy-dev@scipy.org',
          license = 'SciPy License (BSD Style)',
          **configuration(top_path='').todict())
