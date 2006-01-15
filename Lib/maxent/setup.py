#!/usr/bin/env python

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from os.path import join

def configuration(parent_package='', top_path=None):

    config = Configuration('maxent', parent_package, top_path)

    # calc_lwork:
    config.add_extension('logsumexp',
                         [join('src','logsumexp.f')]
                         )

    config.add_data_dir('tests')
    config.add_data_dir('examples')
    config.add_data_dir('doc')

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())

