#!/usr/bin/env python
# Created by Pearu Peterson, August 2002

from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('fftpack',parent_package, top_path)

    config.add_sconscript('SConstruct')
    config.add_data_dir('tests')

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
