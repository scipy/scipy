from __future__ import division, print_function, absolute_import

import os

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from numpy import get_include


def configuration(parent_package='', top_path=None):

    config = Configuration('ndimage', parent_package, top_path)

    include_dirs = ['src',
                    get_include(),
                    os.path.join(os.path.dirname(__file__), '..', '_lib', 'src')]

    config.add_extension("_nd_image",
        sources=["src/nd_image.c","src/ni_filters.c",
                 "src/ni_fourier.c","src/ni_interpolation.c",
                 "src/ni_measure.c",
                 "src/ni_morphology.c","src/ni_support.c"],
        include_dirs=include_dirs)

    # Cython wants the .c and .pyx to have the underscore.
    config.add_extension("_ni_label",
                         sources=["src/_ni_label.c",],
                         include_dirs=['src']+[get_include()])

    config.add_extension("_ctest",
                         sources=["src/_ctest.c"],
                         include_dirs=[get_include()])

    config.add_extension("_ctest_oldapi",
                         sources=["src/_ctest.c"],
                         include_dirs=[get_include()],
                         define_macros=[("OLDAPI", 1)])

    config.add_extension("_cytest",
                         sources=["src/_cytest.c"])
    
    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
