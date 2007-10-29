#!/usr/bin/env python
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import os
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, dict_append
    confgr = Configuration('pyloess',parent_package,top_path)
    # Configuration of LOWESS
    confgr.add_extension('_lowess',
                         sources=[join('src', 'f_lowess.pyf'),
                                  join('src', 'lowess.f'),]
                         )
    # Configuration of STL
    confgr.add_extension('_stl',
                         sources=[join('src', 'f_stl.pyf'),
                                  join('src', 'stl.f')],
                         )
    # Configuration of LOESS
    f_sources = ('loessf.f', 'linpack_lite.f')
    confgr.add_library('floess',
                       sources = [join('src',x) for x in f_sources])
    blas_info = get_info('blas_opt')
    build_info = {}
    dict_append(build_info, **blas_info)
    dict_append(build_info, libraries=['floess'])
    c_sources = ['loess.c', 'loessc.c', 'misc.c', 'predict.c',]
    confgr.add_extension('_loess',
                         sources=[join('src','_loess.c')] + \
                                 [join('src', x) for x in c_sources],
                         depends = [join('src','*.h'),
                                    join('src','*.pyx'),
                                    join('src','*.pxd')
                                    ],
                         **build_info
                        )
    confgr.add_extension('_mloess',
                         sources=[join('src','_mloess.c')] + \
                                 [join('src', x) for x in c_sources],
                         depends = [join('src','*.h'),
                                    join('src','*.pyx'),
                                    join('src','*.pxd')
                                    ],
                         **build_info
                        )
    confgr.add_data_dir('tests')
    return confgr

if __name__ == "__main__":
    from numpy.distutils.core import setup
    config = configuration(top_path='').todict()
    setup(**config)
