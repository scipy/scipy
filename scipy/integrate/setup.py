#!/usr/bin/env python

from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('integrate', parent_package, top_path)

    blas_opt = get_info('blas_opt',notfound_action=2)

    config.add_library('linpack_lite',
                       sources=[join('linpack_lite','*.f')])
    config.add_library('mach',
                       sources=[join('mach','*.f')],
                       config_fc={'noopt':(__file__,1)})
    config.add_library('quadpack',
                       sources=[join('quadpack','*.f')])
    config.add_library('odepack',
                       sources=[join('odepack','*.f')])
    config.add_library('dop',
                       sources=[join('dop','*.f')])
    # should we try to weed through files and replace with calls to
    # LAPACK routines?
    # Yes, someday...



    # Extensions
    # quadpack:

    config.add_extension('_quadpack',
                         sources=['_quadpackmodule.c'],
                         libraries=['quadpack', 'linpack_lite', 'mach'],
                         depends=['quadpack.h','__quadpack.h'])
    # odepack
    libs = ['odepack','linpack_lite','mach']
    

    # Remove libraries key from blas_opt
    if 'libraries' in blas_opt:    # key doesn't exist on OS X ...
        libs.extend(blas_opt['libraries'])
    newblas = {}
    for key in blas_opt.keys():
        if key == 'libraries':
            continue
        newblas[key] = blas_opt[key]
    config.add_extension('_odepack',
                         sources=['_odepackmodule.c'],
                         libraries=libs,
                         depends=['__odepack.h','multipack.h'],
                         **newblas)

    # vode
    config.add_extension('vode',
                         sources=['vode.pyf'],
                         libraries=libs,
                         **newblas)

    # dop
    config.add_extension('_dop',
                         sources=['dop.pyf'],
                         libraries=['dop'])

    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
