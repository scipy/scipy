from os.path import join
from numpy.distutils.system_info import get_info, get_standard_file, \
     BlasNotFoundError


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('csgraph', parent_package, top_path)

    config.add_extension('graph_shortest_path',
         sources=['graph_shortest_path.c'],
         include_dirs=[numpy.get_include()])

    return config
