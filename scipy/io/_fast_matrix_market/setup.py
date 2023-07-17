import os

import pybind11.setup_helpers
from pybind11.setup_helpers import Pybind11Extension


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('_fast_matrix_market', parent_package, top_path)

    define_macros = [
        ('FMM_SCIPY_PRUNE', 1),
        ('FMM_FROM_CHARS_INT_SUPPORTED', 1),
        ('FMM_TO_CHARS_INT_SUPPORTED', 1)
    ]

    if pybind11.setup_helpers.WIN:
        # MSVC has complete <charconv> support so dependencies can be skipped.
        define_macros.append(('FMM_FROM_CHARS_DOUBLE_SUPPORTED', 1))
        define_macros.append(('FMM_FROM_CHARS_LONG_DOUBLE_SUPPORTED', 1))
        define_macros.append(('FMM_TO_CHARS_DOUBLE_SUPPORTED', 1))
        define_macros.append(('FMM_TO_CHARS_LONG_DOUBLE_SUPPORTED', 1))
        # GCC/Clang will use fallbacks as floating-point <charconv> methods are not available.
        # Performance will suffer. Use the meson build instead.

    ext = Pybind11Extension(config.name + "." + '_fmm_core',
                            sources=[os.path.join(config.package_path, "src", s) for s in (
                                '_fmm_core.cpp',
                                '_fmm_core_read_array.cpp',
                                '_fmm_core_read_coo.cpp',
                                '_fmm_core_write_array.cpp',
                                '_fmm_core_write_coo_32.cpp',
                                '_fmm_core_write_coo_64.cpp',
                            )],
                            language='c++',
                            cxx_std=17,
                            include_dirs=[os.path.join(config.package_path, 'fast_matrix_market', 'include')],
                            define_macros=define_macros,
                            )

    config.ext_modules.append(ext)

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
