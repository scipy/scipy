from os.path import join
import pathlib

import numpy as np


def _is_32bit():
    return np.intp(0).itemsize < 8


def check_propack_submodule():
    if not (pathlib.Path(__file__).parent / 'PROPACK/README').exists():
        raise RuntimeError("Missing the `PROPACK` submodule! Run "
                           "`git submodule update --init` to fix this.")


def configuration(parent_package='', top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils import (gfortran_legacy_flag_hook,
                                    get_g77_abi_wrappers,
                                    needs_g77_abi_wrapper)

    lapack_opt = get_info('lapack_opt')
    pre_build_hook = gfortran_legacy_flag_hook
    f2py_options = None

    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    config = Configuration('_propack', parent_package, top_path)

    #  ------------------------------------------------------------
    #  Set up the libraries.
    #  We need a different python extension file for each, because
    #  names reuse between functions in the LAPACK extensions.  This
    #  could probably be remedied with some work.
    #  NOTES: this might not longer apply now that we build without
    #         LAPACK extensions
    type_dict = dict(s='single',
                     d='double',
                     c='complex8',
                     z='complex16')
    check_propack_submodule()

    for prefix, directory in type_dict.items():
        propack_lib = f'_{prefix}propack'
        # Use risc msg implementation for 64-bit machines, pentium for 32-bit
        src_dir = pathlib.Path(__file__).parent / 'PROPACK' / directory
        src = list(src_dir.glob('*.F'))
        if _is_32bit():
            # don't ask me why, 32-bit blows up without second.F
            src = [str(p) for p in src if 'risc' not in str(p)]
        else:
            src = [str(p) for p in src
                   if 'pentium' not in str(p) and 'second' not in str(p)]

        if not _is_32bit():
            # don't ask me why, 32-bit blows up with this wrapper
            src += get_g77_abi_wrappers(lapack_opt)

        cmacros = [('_OPENMP',)]
        if needs_g77_abi_wrapper(lapack_opt):
            cmacros += [('SCIPY_USE_G77_CDOTC_WRAP', 1)]

        config.add_library(propack_lib,
                           sources=src,
                           macros=cmacros,
                           include_dirs=src_dir,
                           depends=['setup.py'])
        ext = config.add_extension(f'_{prefix}propack',
                                   sources=f'{prefix}propack.pyf',
                                   libraries=[propack_lib],
                                   extra_info=lapack_opt,
                                   undef_macros=['_OPENMP'],
                                   f2py_options=f2py_options,
                                   depends=['setup.py'] + src)
        ext._pre_build_hook = pre_build_hook

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
