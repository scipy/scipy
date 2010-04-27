#!/usr/bin/env python
from os.path import join, dirname
import sys
import os

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info

    config = Configuration('dsolve',parent_package,top_path)
    config.add_data_dir('tests')

    lapack_opt = get_info('lapack_opt',notfound_action=2)
    if sys.platform=='win32':
        superlu_defs = [('NO_TIMER',1)]
    else:
        superlu_defs = []
    superlu_defs.append(('USE_VENDOR_BLAS',1))

    superlu_src = os.path.join(dirname(__file__), 'SuperLU', 'SRC')

    config.add_library('superlu_src',
                       sources = [join(superlu_src,'*.c')],
                       macros = superlu_defs,
                       include_dirs=[superlu_src],
                       )

    # Extension
    config.add_extension('_superlu',
                         sources = ['_superlumodule.c',
                                    '_superlu_utils.c',
                                    '_superluobject.c'],
                         libraries = ['superlu_src'],
                         extra_info = lapack_opt,
                         )

    config.add_subpackage('umfpack')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
