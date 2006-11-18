#!/usr/bin/env python

import os

def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError
    from numpy.distutils.misc_util import Configuration
    
    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError,'no lapack/blas resources found'

    config = Configuration('arpack', parent_package, top_path)

    arpack_sources=[os.path.join('ARPACK','SRC', '*.f')]
    arpack_sources.extend([os.path.join('ARPACK','UTIL', '*.f')])
#    arpack_sources.extend([os.path.join('ARPACK','BLAS', '*.f')])
    arpack_sources.extend([os.path.join('ARPACK','LAPACK', '*.f')])

    config.add_library('arpack', sources=arpack_sources)


    config.add_extension('_arpack',
                         sources='arpack.pyf.src',
                         libraries=['arpack'],
                         extra_info = lapack_opt
                        )

    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
