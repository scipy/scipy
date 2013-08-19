#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

from os.path import join

from scipy._build_utils import needs_g77_abi_wrapper


def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError
    from numpy.distutils.misc_util import Configuration

    config = Configuration('arpack',parent_package,top_path)

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    config = Configuration('arpack', parent_package, top_path)

    arpack_sources = [join('ARPACK','SRC', '*.f')]
    arpack_sources.extend([join('ARPACK','UTIL', '*.f')])
    arpack_sources.extend([join('ARPACK','LAPACK', '*.f')])

    if needs_g77_abi_wrapper(lapack_opt):
        arpack_sources += [join('ARPACK', 'FWRAPPERS', 'wrap_veclib_f.f'),
                           join('ARPACK', 'FWRAPPERS', 'wrap_veclib_c.c')]
    else:
        arpack_sources += [join('ARPACK', 'FWRAPPERS', 'wrap_dummy.f')]

    config.add_library('arpack_scipy', sources=arpack_sources,
                       include_dirs=[join('ARPACK', 'SRC')],
                       depends=[join('ARPACK', 'FWRAPPERS',
                                       'wrap_veclib_f.f'),
                                  join('ARPACK', 'FWRAPPERS',
                                       'wrap_veclib_c.c'),
                                  join('ARPACK', 'FWRAPPERS',
                                        'wrap_dummy.f')])

    config.add_extension('_arpack',
                         sources='arpack.pyf.src',
                         libraries=['arpack_scipy'],
                         extra_info=lapack_opt,
                         depends=arpack_sources,
                         )

    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
