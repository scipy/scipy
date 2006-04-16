#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('fdfpack',parent_package,top_path)
    config.add_extension('fdf',
                         sources = ['src/pfdf.c','fdf.pyf'])
    config.add_data_dir('tests')
    config.add_data_dir('utils')
    config.add_data_files('src/32250.pdf')
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
