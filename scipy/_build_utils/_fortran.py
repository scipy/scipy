import re
import os


__all__ = ['needs_g77_abi_wrapper', 'get_g77_abi_wrappers',
           'gfortran_legacy_flag_hook']


def uses_mkl(info):
    r_mkl = re.compile("mkl")
    libraries = info.get('libraries', '')
    for library in libraries:
        if r_mkl.search(library):
            return True

    return False


def needs_g77_abi_wrapper(info):
    """Returns True if g77 ABI wrapper must be used."""
    try:
        needs_wrapper = int(os.environ["SCIPY_USE_G77_ABI_WRAPPER"]) != 0
    except KeyError:
        needs_wrapper = uses_mkl(info)
    return needs_wrapper


def get_g77_abi_wrappers(info):
    """
    Returns file names of source files containing Fortran ABI wrapper
    routines.
    """
    wrapper_sources = []

    path = os.path.abspath(os.path.dirname(__file__))
    if needs_g77_abi_wrapper(info):
        wrapper_sources += [
            os.path.join(path, 'src', 'wrap_g77_abi_f.f'),
            os.path.join(path, 'src', 'wrap_g77_abi_c.c'),
        ]
    else:
        wrapper_sources += [
            os.path.join(path, 'src', 'wrap_dummy_g77_abi.f'),
        ]
    return wrapper_sources


def gfortran_legacy_flag_hook(cmd, ext):
    """
    Pre-build hook to add dd gfortran legacy flag -fallow-argument-mismatch
    """
    from .compiler_helper import try_add_flag
    from distutils.version import LooseVersion

    if isinstance(ext, dict):
        # build_clib
        compilers = ((cmd._f_compiler, ext.setdefault('extra_f77_compile_args', [])),
                      (cmd._f_compiler, ext.setdefault('extra_f90_compile_args', [])))
    else:
        # build_ext
        compilers = ((cmd._f77_compiler, ext.extra_f77_compile_args),
                     (cmd._f90_compiler, ext.extra_f90_compile_args))

    for compiler, args in compilers:
        if compiler is None:
            continue

        if compiler.compiler_type == "gnu95" and compiler.version >= LooseVersion("10"):
            try_add_flag(args, compiler, "-fallow-argument-mismatch")
