import re
import sys

def _uses_veclib(info):
    r_accelerate = re.compile("Accelerate|vecLib")

    extra_link_args = info.get('extra_link_args', '')
    for arg in extra_link_args:
        if r_accelerate.search(arg):
            return True

    return False

def _uses_mkl(info):
    r_mkl = re.compile("mkl_core")

    libraries = info.get('libraries', '')
    for library in libraries:
        if r_mkl.search(library):
            return True

    return False

def needs_g77_abi_wrapper(info):
    """Returns true if g77 ABI wrapper must be used."""
    if _uses_veclib(info):
        return True
    # XXX: is this really true only on Mac OS X ?
    elif _uses_mkl(info) and sys.platform == "darwin":
        return True
    else:
        return False

