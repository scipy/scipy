from __future__ import division, print_function, absolute_import

from os.path import join
from numpy.distutils.system_info import get_info, get_standard_file, \
     BlasNotFoundError


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('csgraph', parent_package, top_path)

    config.add_data_dir('tests')

    config.add_extension('_shortest_path',
         sources=['_shortest_path.c'],
         include_dirs=[numpy.get_include()])

    config.add_extension('_traversal',
         sources=['_traversal.c'],
         include_dirs=[numpy.get_include()])

    config.add_extension('_min_spanning_tree',
         sources=['_min_spanning_tree.c'],
         include_dirs=[numpy.get_include()])

    config.add_extension('_tools',
         sources=['_tools.c'],
         include_dirs=[numpy.get_include()])

    return config
