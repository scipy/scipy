#!/usr/bin/env python

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from os.path import join

def configuration(parent_package='', top_path=None):

    config = Configuration('maxentropy', parent_package, top_path)

    config.add_data_dir('tests')
    config.add_data_dir('examples')

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
