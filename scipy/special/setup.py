#!/usr/bin/env python

import os
import sys
from os.path import join
from distutils.sysconfig import get_python_inc
import numpy
from numpy.distutils.misc_util import get_numpy_include_dirs

try:
    from numpy.distutils.misc_util import get_info
except ImportError:
    raise ValueError("numpy >= 1.4 is required (detected %s from %s)" % \
                     (numpy.__version__, numpy.__file__))

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('special', parent_package, top_path)

    define_macros = []
    if sys.platform=='win32':
#        define_macros.append(('NOINFINITIES',None))
#        define_macros.append(('NONANS',None))
        define_macros.append(('_USE_MATH_DEFINES',None))

    # C libraries
    config.add_library('sc_c_misc',sources=[join('c_misc','*.c')],
                       include_dirs=[get_python_inc(), get_numpy_include_dirs()],
                       macros=define_macros)
    config.add_library('sc_cephes',sources=[join('cephes','*.c')],
                       include_dirs=[get_python_inc(), get_numpy_include_dirs()],
                       macros=define_macros)

    # Fortran libraries
    config.add_library('sc_mach',sources=[join('mach','*.f')],
                       config_fc={'noopt':(__file__,1)})
    config.add_library('sc_toms',sources=[join('amos','*.f')])
    config.add_library('sc_amos',sources=[join('toms','*.f')])
    config.add_library('sc_cdf',sources=[join('cdflib','*.f')])
    config.add_library('sc_specfun',sources=[join('specfun','*.f')])

    # Extension _cephes
    sources = ['_cephesmodule.c', 'amos_wrappers.c', 'specfun_wrappers.c',
               'toms_wrappers.c','cdf_wrappers.c','ufunc_extras.c']
    config.add_extension('_cephes', sources=sources,
                         libraries=['sc_amos','sc_toms','sc_c_misc','sc_cephes','sc_mach',
                                    'sc_cdf', 'sc_specfun'],
                         depends=["ufunc_extras.h", "cephes.h",
                                  "amos_wrappers.h", "toms_wrappers.h",
                                  "cdf_wrappers.h", "specfun_wrappers.h",
                                  "c_misc/misc.h", "cephes_doc.h",
                                  "cephes/mconf.h", "cephes/cephes_names.h"],
                         define_macros = define_macros,
                         extra_info=get_info("npymath")
                         )

    # Extension specfun
    config.add_extension('specfun',
                         sources=['specfun.pyf'],
                         f2py_options=['--no-wrap-functions'],
                         define_macros=[],
                         libraries=['sc_specfun'])

    # Extension orthogonal_eval
    config.add_extension('orthogonal_eval',
                         sources=['orthogonal_eval.c'],
                         define_macros=[],
                         extra_info=get_info("npymath"))

    # Extension lambertw
    config.add_extension('lambertw',
                         sources=['lambertw.c'],
                         define_macros=[],
                         extra_info=get_info("npymath"))

    config.add_data_files('tests/*.py')
    config.add_data_files('tests/data/README')
    config.add_data_files('tests/data/*.npz')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
