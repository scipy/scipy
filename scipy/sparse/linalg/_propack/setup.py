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
    config = Configuration('_propack', parent_package, top_path)
    lapack_opt = get_info('lapack_opt')
    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    #  ------------------------------------------------------------
    #  Set up the libraries.
    #  We need a different python extension file for each, because
    #  names resue between functions in the LAPACK extensions.  This
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

        # Need to use risc implementation for 32-bit machines
        if _is_32bit:
            src = list((pathlib.Path(
                __file__).parent / 'PROPACK' / directory).glob('*.F'))
            src = [str(p) for p in src if 'risc' not in str(p)]
        else:
            src = join('PROPACK', directory, '*.F')

        config.add_library(propack_lib,
                           sources=src,
                           macros=[('_OPENMP',)])
        config.add_extension(f'_{prefix}propack',
                             sources=f'{prefix}propack.pyf',
                             libraries=[propack_lib],
                             extra_info=lapack_opt,
                             undef_macros=['_OPENMP'])

        # add required data files to run example matrix tests
        path_list = ['PROPACK', directory, 'Examples']
        config.add_data_files('.', join(*path_list, '*.coord'))
        config.add_data_files('.', join(*path_list, '*.diag'))
        config.add_data_files('.', join(*path_list, '*.rra'))
        config.add_data_files('.', join(*path_list, '*.cua'))
        config.add_data_files('.', join(*path_list, 'Output', '*.ascii'))

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
