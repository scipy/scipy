import os
from os.path import join


def pre_build_hook(build_ext, ext):
    # Copy & paste from scipy.spatial
    from scipy._build_utils.compiler_helper import (set_cxx_flags_hook,
                                                    try_add_flag)
    cc = build_ext._cxx_compiler
    args = ext.extra_compile_args

    set_cxx_flags_hook(build_ext, ext)

    if cc.compiler_type == 'msvc':
        # Ignore "structured exceptions" which are non-standard MSVC extensions
        args.append('/EHsc')
    else:
        # Don't export library symbols
        try_add_flag(args, cc, '-fvisibility=hidden')


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    from scipy._build_utils import (get_f2py_int64_options,
                                    ilp64_pre_build_hook,
                                    uses_blas64, numpy_nodepr_api)
    import pybind11

    if uses_blas64():
        # TODO: Note that fitpack does not use BLAS/LAPACK.
        # The reason why we use 64-bit ints only in this case
        # is because scipy._build_utils knows the 64-bit int
        # flags for too few Fortran compilers, so we cannot turn
        # this on by default.
        pre_build_hook = ilp64_pre_build_hook
        f2py_options = get_f2py_int64_options()
        define_macros = [("HAVE_ILP64", None)]
    else:
        pre_build_hook = None
        f2py_options = None
        define_macros = []

    config = Configuration('interpolate', parent_package, top_path)

    fitpack_src = [join('fitpack', '*.f')]
    config.add_library('fitpack', sources=fitpack_src,
                       _pre_build_hook=pre_build_hook)

    config.add_extension('interpnd',
                         sources=['interpnd.c'])

    config.add_extension('_rgi_cython',
                         sources=['_rgi_cython.c'])

    config.add_extension('_ppoly',
                         sources=['_ppoly.c'])

    config.add_extension('_bspl',
                         sources=['_bspl.c'],
                         depends=['src/__fitpack.h'])

    config.add_extension('_fitpack',
                         sources=['src/_fitpackmodule.c'],
                         libraries=['fitpack'],
                         define_macros=define_macros + numpy_nodepr_api['define_macros'],
                         depends=(['src/__fitpack.h']
                                  + fitpack_src)
                         )

    config.add_extension('dfitpack',
                         sources=['src/fitpack.pyf'],
                         libraries=['fitpack'],
                         define_macros=define_macros,
                         depends=fitpack_src,
                         f2py_options=f2py_options
                         )

    if int(os.environ.get('SCIPY_USE_PYTHRAN', 1)):
        from pythran.dist import PythranExtension
        ext = PythranExtension(
            'scipy.interpolate._rbfinterp_pythran',
            sources=['scipy/interpolate/_rbfinterp_pythran.py'],
            config=['compiler.blas=none']
            )
        config.ext_modules.append(ext)
    
    pava_pybind_includes = [
        pybind11.get_include(True),
        pybind11.get_include(False),
        get_numpy_include_dirs()]
    ext = config.add_extension('_pava_pybind',
                               sources=[join('src', 'pava_pybind.cpp')],
                               depends=[],
                               include_dirs=pava_pybind_includes,
                               language='c++',
                               **numpy_nodepr_api)
    ext._pre_build_hook = pre_build_hook

    config.add_data_dir('tests')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
