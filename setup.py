#!/usr/bin/env python

from numpy import get_include

import os
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('', parent_package, top_path)


    # C++ extension for several basic interpolation types
    config.add_extension('_interpolate',
                        ['extensions/_interpolate.cpp'],
                        include_dirs = ['extensions'],
                        depends = ['extensions/interpolate.h'])

    # used by dfitpack extension
    config.add_library('_fitpack',
                        sources=[join('extensions/fitpack', '*.f')],
                      )

    # Fortran routines (collectively "FITPACK" for spline interpolation)
    config.add_extension('_dfitpack',
                        sources=['extensions/_fitpack.pyf'],
                        libraries=['_fitpack'],
                        )
                        
    # ND Image routines for ND interpolation
    config.add_extension('_nd_image',
                        sources = ["extensions/ndimage/nd_image.c",
                                        "extensions/ndimage/ni_filters.c",
                                        "extensions/ndimage/ni_fourier.c",
                                        "extensions/ndimage/ni_interpolation.c",
                                        "extensions/ndimage/ni_measure.c",
                                        "extensions/ndimage/ni_morphology.c",
                                        "extensions/ndimage/ni_support.c"],
                        include_dirs=['extensions/ndimage']+[get_include()],
                        )
    
    # implements algorithm 526 for 2D interpolation
    config.add_extension('_interp_526',
                        sources = ['extensions/interp_526a.pyf',
                                        'extensions/interp_526.f']
                        )
                        
    config.add_data_dir('docs')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration().todict())