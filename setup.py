#!/usr/bin/env python

from numpy import get_include

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
                        
    # ND Image routines for ND interpolation
    config.add_extension('_nd_image',
                        sources=["ndimage/nd_image.c",
                                    "ndimage/ni_filters.c",
                                    "ndimage/ni_fourier.c",
                                    "ndimage/ni_interpolation.c",
                                    "ndimage/ni_measure.c",
                                    "ndimage/ni_morphology.c",
                                    "ndimage/ni_support.c"],
                        include_dirs=['ndimage']+[get_include()],
                        )
                        
    config.add_data_dir('docs')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration().todict())