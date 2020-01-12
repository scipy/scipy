import re
import os
import glob
from distutils.dep_util import newer


__all__ = ['needs_g77_abi_wrapper', 'get_g77_abi_wrappers']


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
