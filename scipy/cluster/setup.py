#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import sys

from os.path import join

if sys.version_info[0] >= 3:
    DEFINE_MACROS = [("SCIPY_PY3K", None)]
else:
    DEFINE_MACROS = []


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    config = Configuration('cluster', parent_package, top_path)

    config.add_data_dir('tests')

    config.add_extension('_vq',
        sources=[join('src', 'vq_module.c'), join('src', 'vq.c')],
        include_dirs=[get_numpy_include_dirs()],
        define_macros=DEFINE_MACROS)

    config.add_extension('_hierarchy_wrap',
        sources=[join('src', 'hierarchy_wrap.c'), join('src', 'hierarchy.c')],
        include_dirs=[get_numpy_include_dirs()],
        define_macros=DEFINE_MACROS)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(maintainer="SciPy Developers",
          author="Eric Jones",
          maintainer_email="scipy-dev@scipy.org",
          description="Clustering Algorithms (Information Theory)",
          url="http://www.scipy.org",
          license="SciPy License (BSD Style)",
          **configuration(top_path='').todict()
          )
