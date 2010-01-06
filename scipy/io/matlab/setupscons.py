#!/usr/bin/env python

def configuration(parent_package='io',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('matlab', parent_package, top_path)
    config.add_sconscript('SConstruct')
    config.add_data_dir('tests')
    config.add_data_dir('benchmarks')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
