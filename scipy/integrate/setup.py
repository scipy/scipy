#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

from os.path import join


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('integrate', parent_package, top_path)

    blas_opt = get_info('blas_opt',notfound_action=2)
    lapack_opt = get_info('lapack_opt',notfound_action=2)

    mach_src = [join('mach','*.f')]
    quadpack_src = [join('quadpack','*.f')]
    odepack_src = [join('odepack','*.f')]
    dop_src = [join('dop','*.f')]
    quadpack_test_src = [join('tests','_test_multivariate.c')]

    config.add_library('mach', sources=mach_src,
                       config_fc={'noopt':(__file__,1)})
    config.add_library('quadpack', sources=quadpack_src)
    config.add_library('odepack', sources=odepack_src)
    config.add_library('dop', sources=dop_src)

    # Extensions
    # quadpack:
    config.add_extension('_quadpack',
                         sources=['_quadpackmodule.c'],
                         libraries=(['quadpack', 'mach'] +
                                    lapack_opt['libraries']),
                         depends=(['quadpack.h','__quadpack.h']
                                  + quadpack_src + mach_src))

    # odepack
    libs = ['odepack','mach']

    # Remove libraries key from blas_opt
    if 'libraries' in blas_opt:    # key doesn't exist on OS X ...
        libs.extend(blas_opt['libraries'])
    libs.extend(lapack_opt['libraries'])
    newblas = {}
    for key in blas_opt:
        if key == 'libraries':
            continue
        newblas[key] = blas_opt[key]
    # TODO add LAPACK stuff to newblas
    config.add_extension('_odepack',
                         sources=['_odepackmodule.c'],
                         libraries=libs,
                         depends=(['__odepack.h','multipack.h']
                                  + odepack_src
                                  + mach_src),
                         **newblas)

    # vode
    config.add_extension('vode',
                         sources=['vode.pyf'],
                         libraries=libs,
                         depends=(odepack_src
                                  + mach_src),
                         **newblas)

    # lsoda
    config.add_extension('lsoda',
                         sources=['lsoda.pyf'],
                         libraries=libs,
                         depends=(odepack_src
                                  + mach_src),
                         **newblas)

    # dop
    config.add_extension('_dop',
                         sources=['dop.pyf'],
                         libraries=['dop'],
                         depends=dop_src)

    config.add_extension('_test_multivariate',
                         sources=quadpack_test_src)
    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
