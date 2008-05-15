#!/usr/bin/env python
# Created by Pearu Peterson, August 2002

from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('fftpack',parent_package, top_path)

    config.add_library('dfftpack',
                       sources=[join('DFFTPACK','*.f')])

    sources = ['fftpack.pyf', 'src/zrfft.c']
    for s in ["zfft.cxx", "zfftnd.cxx", "drfft.cxx"]:
        sources.append(join('src/fftpack', s))

    # Build the python extensions
    config.add_extension('_fftpack',
        sources=sources,
        libraries = ["dfftpack"],
        include_dirs = ['src'],
    )

    config.add_extension('convolve',
        sources = ['convolve.pyf', 'src/fftpack/convolve.cxx'],
        libraries = ["dfftpack"],
        include_dirs = ['src'],
    )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
