#!/usr/bin/env python

from os.path import join

def configuration(parent_package = '', top_path = None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    from numpy.distutils.system_info import get_info
    from distutils.sysconfig import get_python_inc

    config = Configuration('spatial', parent_package, top_path)

    config.add_data_dir('tests')

    qhull_src = ['geom2.c', 'geom.c', 'global.c', 'io.c', 'libqhull.c',
                 'mem.c', 'merge.c', 'poly2.c', 'poly.c', 'qset.c',
                 'random.c', 'rboxlib.c', 'stat.c', 'user.c', 'usermem.c',
                 'userprintf.c']

    config.add_library('qhull',
                       sources=[join('qhull', 'src', x) for x in qhull_src],
                       include_dirs=[get_python_inc(),
                                     get_numpy_include_dirs()],
                       # XXX: GCC dependency!
                       #extra_compiler_args=['-fno-strict-aliasing'],
                       )

    lapack = dict(get_info('lapack_opt'))
    try:
        libs = ['qhull'] + lapack.pop('libraries')
    except KeyError:
        libs = ['qhull']
    config.add_extension('qhull',
                         sources=['qhull.c'],
                         libraries=libs,
                         **lapack)

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
