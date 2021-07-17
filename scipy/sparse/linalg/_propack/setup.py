from os.path import join
import pathlib


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

    #------------------------------------------------------------
    # Set up the libraries.
    #  We need a different python extension file for each, because
    #  names resue between functions in the LAPACK extensions.  This
    #  could probably be remedied with some work.
    type_dict = dict(s='single',
                     d='double',
                     c='complex8',
                     z='complex16')
    check_propack_submodule()

    for prefix, directory in type_dict.items():
        lapack_lib = f'sp_{prefix}lapack_util'
        propack_lib = f'_{prefix}propack'

        lapack_sources = join('PROPACK', directory, 'Lapack_Util', '*.f')
        propack_sources = join('PROPACK', directory, '*.F')

        config.add_library(lapack_lib,
                           sources=lapack_sources,
                           macros=[('_OPENMP',)])
        config.add_library(propack_lib,
                           sources=propack_sources,
                           macros=[('_OPENMP',)])
        config.add_extension(f'_{prefix}propack',
                             sources=f'{prefix}propack.pyf',
                             libraries=[lapack_lib, propack_lib],
                             extra_info=lapack_opt,
                             undef_macros=['_OPENMP'])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
