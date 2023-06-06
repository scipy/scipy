def configuration(parent_package='', top_path=None):
    import sys
    from numpy.distutils.misc_util import Configuration
    config = Configuration('_fast_matrix_market', parent_package, top_path)

    define_macros = [
        ('FMM_FROM_CHARS_INT_SUPPORTED', 1),
        ('FMM_TO_CHARS_INT_SUPPORTED', 1)
    ]
    extra_compile_args = []

    if sys.platform.startswith('win'):
        # MSVC
        extra_compile_args.append('/std:c++17')

        # MSVC has complete <charconv> support so dependencies can be skipped.
        define_macros.append(('FMM_FROM_CHARS_DOUBLE_SUPPORTED', 1))
        define_macros.append(('FMM_FROM_CHARS_LONG_DOUBLE_SUPPORTED', 1))
        define_macros.append(('FMM_TO_CHARS_DOUBLE_SUPPORTED', 1))
        define_macros.append(('FMM_TO_CHARS_LONG_DOUBLE_SUPPORTED', 1))
    else:
        # GCC/Clang
        extra_compile_args.append('-std=c++17')
        # Will use fallbacks as floating-point <charconv> methods are not available.
        # Performance will suffer. Use the meson build instead.

    config.add_extension('_core',
                         sources=['src/_core.cpp'],
                         include_dirs=['fast_matrix_market/include'],
                         define_macros=define_macros,
                         extra_compile_args=extra_compile_args,
                         )

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
