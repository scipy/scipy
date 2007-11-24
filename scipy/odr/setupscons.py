#!/usr/bin/env python

from os.path import join

def configuration(parent_package='', top_path=None):
    import warnings
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, BlasNotFoundError
    config = Configuration('odr', parent_package, top_path)

    config.add_sconscript('SConstruct')
    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
