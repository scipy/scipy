#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

from os.path import join


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('integrate', parent_package, top_path)

    # Get a local copy of lapack_opt_info
    lapack_opt = dict(get_info('lapack_opt',notfound_action=2))
    # Pop off the libraries list so it can be combined with
    # additional required libraries
    lapack_libs = lapack_opt.pop('libraries', [])

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
                         libraries=(['quadpack', 'mach'] + lapack_libs),
                         depends=(['quadpack.h','__quadpack.h']
                                  + quadpack_src + mach_src),
                         **lapack_opt)

    # odepack
    odepack_libs = ['odepack','mach'] + lapack_libs
    
    config.add_extension('_odepack',
                         sources=['_odepackmodule.c'],
                         libraries=odepack_libs,
                         depends=(['__odepack.h','multipack.h']
                                  + odepack_src
                                  + mach_src),
                         **lapack_opt)

    # vode
    config.add_extension('vode',
                         sources=['vode.pyf'],
                         libraries=odepack_libs,
                         depends=(odepack_src
                                  + mach_src),
                         **lapack_opt)

    # lsoda
    config.add_extension('lsoda',
                         sources=['lsoda.pyf'],
                         libraries=odepack_libs,
                         depends=(odepack_src
                                  + mach_src),
                         **lapack_opt)

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
