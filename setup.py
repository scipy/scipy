#!/usr/bin/env python

import os
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('', parent_package, top_path)


    # C++ extension for several basic interpolation types
    config.add_extension('_interpolate',
                         ['_interpolate.cpp'],
                         include_dirs = ['.'],
                         depends = ['interpolate.h'])

    # used by dfitpack extension
    config.add_library('_fitpack',
                       sources=[join('fitpack', '*.f')],
                      )

    # Fortran routines (collectively "FITPACK" for spline interpolation)
    config.add_extension('_dfitpack',
                         sources=['_fitpack.pyf'],
                         libraries=['_fitpack'],
                        )
                        
    # FIXME : add documentation files
    # config.add_data_dir(

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration().todict())