#!/usr/bin/env python

from os.path import join

def configuration(parent_package = '', top_path = None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    config = Configuration('spatial', parent_package, top_path)

    config.add_data_dir('tests')

    config.add_extension('ckdtree', sources=['ckdtree.c']) # FIXME: cython

    config.add_extension('_distance_wrap',
        sources=[join('src', 'distance_wrap.c'), join('src', 'distance.c')],
        include_dirs = [get_numpy_include_dirs()])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(maintainer = "SciPy Developers",
          author = "Anne Archibald",
          maintainer_email = "scipy-dev@scipy.org",
          description = "Spatial algorithms and data structures",
          url = "http://www.scipy.org",
          license = "SciPy License (BSD Style)",
          **configuration(top_path='').todict()
          )
