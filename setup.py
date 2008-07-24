#!/usr/bin/env python

import os
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('', parent_package, top_path)


    config.add_extension('_interpolate',
                         ['_interpolate.cpp'],
                         include_dirs = ['.'],
                         depends = ['interpolate.h'])

    config.add_library('fitpack',
                       sources=[join('fitpack', '*.f')],
                      )

    config.add_extension('dfitpack',
                         sources=['fitpack.pyf'],
                         libraries=['fitpack'],
                        )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration().todict())