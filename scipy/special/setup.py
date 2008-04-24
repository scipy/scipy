#!/usr/bin/env python

import os
import sys
from os.path import join
from distutils.sysconfig import get_python_inc

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('special', parent_package, top_path)

    define_macros = []
    if sys.platform=='win32':
#        define_macros.append(('NOINFINITIES',None))
#        define_macros.append(('NONANS',None))
        define_macros.append(('_USE_MATH_DEFINES',None))

    # C libraries
    config.add_library('c_misc',sources=[join('c_misc','*.c')])
    config.add_library('cephes',sources=[join('cephes','*.c')],
                       include_dirs=[get_python_inc()],
                       macros=define_macros)

    # Fortran libraries
    config.add_library('mach',sources=[join('mach','*.f')],
                       config_fc={'noopt':(__file__,1)})
    config.add_library('toms',sources=[join('amos','*.f')])
    config.add_library('amos',sources=[join('toms','*.f')])
    config.add_library('cdf',sources=[join('cdflib','*.f')])
    config.add_library('specfun',sources=[join('specfun','*.f')])

    # Extension _cephes
    sources = ['_cephesmodule.c', 'amos_wrappers.c', 'specfun_wrappers.c',
               'toms_wrappers.c','cdf_wrappers.c','ufunc_extras.c']
    config.add_extension('_cephes', sources=sources,
                         libraries=['amos','toms','c_misc','cephes','mach',
                                    'cdf', 'specfun'],
                         depends=["ufunc_extras.h", "cephes.h",
                                  "amos_wrappers.h", "toms_wrappers.h",
                                  "cdf_wrappers.h", "specfun_wrappers.h",
                                  "c_misc/misc.h", "cephes_doc.h",
                                  "cephes/mconf.h", "cephes/cephes_names.h"],
                         define_macros = define_macros
                         )

    # Extension specfun
    config.add_extension('specfun',
                         sources=['specfun.pyf'],
                         f2py_options=['--no-wrap-functions'],
                         define_macros=[],
                         libraries=['specfun'])

    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
