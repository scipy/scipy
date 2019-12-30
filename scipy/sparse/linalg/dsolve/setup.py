from __future__ import division, print_function, absolute_import

from os.path import join, dirname
import sys
import os
import glob


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils.system_info import get_info
    from scipy._build_utils import numpy_nodepr_api, combine_dict, uses_blas64

    config = Configuration('dsolve',parent_package,top_path)
    config.add_data_dir('tests')

    if sys.platform == 'win32':
        superlu_defs = [('NO_TIMER',1)]
    else:
        superlu_defs = []
    superlu_defs.append(('USE_VENDOR_BLAS',1))

    superlu_src = join(dirname(__file__), 'SuperLU', 'SRC')

    sources = sorted(glob.glob(join(superlu_src, '*.c')))
    headers = list(glob.glob(join(superlu_src, '*.h')))

    if uses_blas64():
        lapack_opt = get_info('lapack_ilp64_opt',notfound_action=2)
        lapack_opt = combine_dict(lapack_opt,
                                  define_macros=[('BLASLAPACK6432_SUFFIX', 'superlu_blas_')])
        sources.append(join(dirname(__file__), '..', '..', '..', '_build_utils',
                            'src', 'blaslapack6432.c'))
    else:
        lapack_opt = get_info('lapack_opt',notfound_action=2)

    cfg = combine_dict(lapack_opt, numpy_nodepr_api,
                       define_macros=superlu_defs,
                       include_dirs=[superlu_src,
                                     join(dirname(__file__), 'SuperLU')])

    config.add_library('superlu_src',
                       sources=sources,
                       macros=cfg['define_macros'],
                       include_dirs=cfg['include_dirs'])

    # Extension
    ext_sources = ['_superlumodule.c',
                   '_superlu_utils.c',
                   '_superluobject.c']

    cfg = combine_dict(cfg, libraries=['superlu_src'])
    config.add_extension('_superlu',
                         sources=ext_sources,
                         depends=(sources + headers),
                         **cfg
                         )

    # Add license files
    config.add_data_files('SuperLU/License.txt')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
