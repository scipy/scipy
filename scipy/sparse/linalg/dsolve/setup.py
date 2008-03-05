#!/usr/bin/env python
from os.path import join
import sys

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

    config.add_library('superlu_src',
                       sources = [join('SuperLU','SRC','*.c')],
                       macros = superlu_defs
                       )

    #SuperLU/SRC/util.h  has been modifed to use these by default
    #macs = [('USER_ABORT','superlu_python_module_abort'),
    #        ('USER_MALLOC','superlu_python_module_malloc'),
    #        ('USER_FREE','superlu_python_module_free')]

    # Extension
    config.add_extension('_zsuperlu',
                         sources = ['_zsuperlumodule.c','_superlu_utils.c',
                                    '_superluobject.c'],
                         libraries = ['superlu_src'],
                         extra_info = lapack_opt
                         )

    config.add_extension('_dsuperlu',
                         sources = ['_dsuperlumodule.c','_superlu_utils.c',
                                    '_superluobject.c'],
                         libraries = ['superlu_src'],
                         extra_info = lapack_opt
                         )

    config.add_extension('_csuperlu',
                         sources = ['_csuperlumodule.c','_superlu_utils.c',
                                    '_superluobject.c'],
                         libraries = ['superlu_src'],
                         extra_info = lapack_opt
                         )

    config.add_extension('_ssuperlu',
                         sources = ['_ssuperlumodule.c','_superlu_utils.c',
                                    '_superluobject.c'],
                         libraries = ['superlu_src'],
                         extra_info = lapack_opt
                         )

    config.add_subpackage('umfpack')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
