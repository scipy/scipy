#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

from os.path import join


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('integrate', parent_package, top_path)

    blas_opt = get_info('blas_opt',notfound_action=2)

    linpack_lite_src = [join('linpack_lite','*.f')]
    mach_src = [join('mach','*.f')]
    quadpack_src = [join('quadpack','*.f')]
    odepack_src = [join('odepack','*.f')]
    dop_src = [join('dop','*.f')]

    config.add_library('linpack_lite', sources=linpack_lite_src)
    config.add_library('mach', sources=mach_src,
                       config_fc={'noopt':(__file__,1)})
    config.add_library('quadpack', sources=quadpack_src)
    config.add_library('odepack', sources=odepack_src)
    config.add_library('dop', sources=dop_src)
    # should we try to weed through files and replace with calls to
    # LAPACK routines?
    # Yes, someday...

    # Extensions
    # quadpack:

    config.add_extension('_quadpack',
                         sources=['_quadpackmodule.c'],
                         libraries=['quadpack', 'linpack_lite', 'mach'],
                         depends=(['quadpack.h','__quadpack.h']
                                  + quadpack_src + linpack_lite_src + mach_src))
    # odepack
    libs = ['odepack','linpack_lite','mach']

    # Remove libraries key from blas_opt
    if 'libraries' in blas_opt:    # key doesn't exist on OS X ...
        libs.extend(blas_opt['libraries'])
    newblas = {}
    for key in blas_opt:
        if key == 'libraries':
            continue
        newblas[key] = blas_opt[key]
    config.add_extension('_odepack',
                         sources=['_odepackmodule.c'],
                         libraries=libs,
                         depends=(['__odepack.h','multipack.h']
                                  + odepack_src + linpack_lite_src
                                  + mach_src),
                         **newblas)

    # vode
    config.add_extension('vode',
                         sources=['vode.pyf'],
                         libraries=libs,
                         depends=(odepack_src + linpack_lite_src
                                  + mach_src),
                         **newblas)

    # lsoda
    config.add_extension('lsoda',
                         sources=['lsoda.pyf'],
                         libraries=libs,
                         depends=(odepack_src + linpack_lite_src
                                  + mach_src),
                         **newblas)

    # dop
    config.add_extension('_dop',
                         sources=['dop.pyf'],
                         libraries=['dop'],
                         depends=dop_src)

    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
