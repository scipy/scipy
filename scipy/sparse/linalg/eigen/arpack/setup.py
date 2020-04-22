from os.path import join


def configuration(parent_package='',top_path=None):
    from scipy._build_utils.system_info import get_info
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils import (get_g77_abi_wrappers,
                                    gfortran_legacy_flag_hook)

    lapack_opt = get_info('lapack_opt')

    config = Configuration('arpack', parent_package, top_path)

    arpack_sources = [join('ARPACK','SRC', '*.f')]
    arpack_sources.extend([join('ARPACK','UTIL', '*.f')])

    arpack_sources += get_g77_abi_wrappers(lapack_opt)

    config.add_library('arpack_scipy', sources=arpack_sources,
                       include_dirs=[join('ARPACK', 'SRC')],
                       _pre_build_hook=gfortran_legacy_flag_hook)

    ext_sources = ['arpack.pyf.src']
    ext = config.add_extension('_arpack',
                               sources=ext_sources,
                               libraries=['arpack_scipy'],
                               extra_info=lapack_opt,
                               depends=arpack_sources,
                               )
    ext._pre_build_hook = gfortran_legacy_flag_hook

    config.add_data_dir('tests')

    # Add license files
    config.add_data_files('ARPACK/COPYING')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
