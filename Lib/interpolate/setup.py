#!/usr/bin/env python

import os

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('interpolate', parent_package, top_path)

    config.add_library('fitpack',
                       sources=[os.path.join('fitpack', '*.f')],
                      )

    config.add_extension('_fitpack',
                         sources=['_fitpackmodule.c'],
                         libraries=['fitpack'],
                        )

    config.add_extension('dfitpack',
                         sources=['fitpack.pyf'],
                         libraries=['fitpack'],
                        )

    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
