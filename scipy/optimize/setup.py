#!/usr/bin/env python

from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('optimize',parent_package, top_path)

    config.add_library('minpack',sources=[join('minpack','*f')])
    config.add_extension('_minpack',
                         sources=['_minpackmodule.c'],
                         libraries=['minpack'],
                         depends=["minpack.h","__minpack.h"])

    config.add_library('rootfind',
                       sources=[join('Zeros','*.c')],
                       headers=[join('Zeros','zeros.h')])

    config.add_extension('_zeros',
                         sources=['zeros.c'],
                         libraries=['rootfind'])

    lapack = get_info('lapack_opt')
    sources=['lbfgsb.pyf','routines.f']
    config.add_extension('_lbfgsb',
                         sources=[join('lbfgsb',x) for x in sources],
                         **lapack)

    sources=['moduleTNC.c','tnc.c']
    config.add_extension('moduleTNC',
                         sources=[join('tnc',x) for x in sources],
                         depends=[join('tnc','tnc.h')])

    config.add_extension('_cobyla',
                         sources=[join('cobyla',x) for x in ['cobyla.pyf',
                                                             'cobyla2.f',
                                                             'trstlp.f']])
    sources = ['minpack2.pyf', 'dcsrch.f', 'dcstep.f']
    config.add_extension('minpack2',
                         sources=[join('minpack2',x) for x in sources])

    sources = ['slsqp.pyf', 'slsqp_optmz.f']
    config.add_extension('_slsqp', sources=[join('slsqp', x) for x in sources])

    config.add_extension('_nnls', sources=[join('nnls', x) \
                                          for x in ["nnls.f","nnls.pyf"]])

    config.add_data_dir('tests')
    config.add_data_dir('benchmarks')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
