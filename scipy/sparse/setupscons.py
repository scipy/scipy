#!/usr/bin/env python

from os.path import join
import sys

def configuration(parent_package='',top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('sparse',parent_package,top_path,
                           setup_name = 'setupscons.py')

    config.add_data_dir('tests')
    config.add_subpackage('linalg')
    config.add_subpackage('sparsetools')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
