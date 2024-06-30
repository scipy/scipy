"""
Helper functions and variables for generation of BLAS/LAPACK wrappers.
"""

import os
from stat import ST_MTIME

# Used to convert from types in signature files to C types
C_TYPES = {'int': 'int',
           'c': 'npy_complex64',
           'd': 'double',
           's': 'float',
           'z': 'npy_complex128',
           'char': 'char',
           'bint': 'int',
           'void': 'void',
           'cselect1': '_cselect1',
           'cselect2': '_cselect2',
           'dselect2': '_dselect2',
           'dselect3': '_dselect3',
           'sselect2': '_sselect2',
           'sselect3': '_sselect3',
           'zselect1': '_zselect1',
           'zselect2': '_zselect2'}

# Used to convert complex types in signature files to Numpy complex types
NPY_TYPES = {'c': 'npy_complex64', 'z': 'npy_complex128',
             'cselect1': '_cselect1', 'cselect2': '_cselect2',
             'dselect2': '_dselect2', 'dselect3': '_dselect3',
             'sselect2': '_sselect2', 'sselect3': '_sselect3',
             'zselect1': '_zselect1', 'zselect2': '_zselect2'}

# BLAS/LAPACK functions with complex return values (use 'wrp'-suffixed
# wrappers from G77 ABI wrapper)
WRAPPED_FUNCS = ['cdotc', 'cdotu', 'zdotc', 'zdotu', 'cladiv', 'zladiv']

# Missing from new Accelerate so use standard old Accelerate symbols
USE_OLD_ACCELERATE = ['lsame', 'dcabs1']

C_PREAMBLE = """
#include "npy_cblas.h"
#include "fortran_defs.h"
"""

LAPACK_DECLS = """
typedef int (*_cselect1)(npy_complex64*);
typedef int (*_cselect2)(npy_complex64*, npy_complex64*);
typedef int (*_dselect2)(double*, double*);
typedef int (*_dselect3)(double*, double*, double*);
typedef int (*_sselect2)(float*, float*);
typedef int (*_sselect3)(float*, float*, float*);
typedef int (*_zselect1)(npy_complex128*);
typedef int (*_zselect2)(npy_complex128*, npy_complex128*);
"""

CPP_GUARD_BEGIN = """
#ifdef __cplusplus
extern "C" {
#endif

"""

CPP_GUARD_END = """
#ifdef __cplusplus
}
#endif
"""


def read_signatures(lines):
    """
    Read BLAS/LAPACK signatures and split into name, return type, argument
    names, and argument types.
    """
    sigs = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        line = line[:-1].split('(')
        args = line[1]
        name_and_type = line[0].split(' ')
        ret_type = name_and_type[0]
        name = name_and_type[1]
        argtypes, argnames = zip(*[arg.split(' *') for arg in args.split(', ')])
        # Argname cannot be same as abbreviated return type
        if ret_type in argnames:
            argnames = [n if n != ret_type else n + '_' for n in argnames]
        # Argname should not be Python keyword
        argnames = [n if n not in ['lambda', 'in'] else n + '_' for n in argnames]
        sigs.append({
            'name': name,
            'return_type': ret_type,
            'argnames': argnames,
            'argtypes': list(argtypes)
        })
    return sigs

def newer(dst, src):
    """
    Return true if 'dst' exists and is more recently modified than
    'src', or if 'dst' exists and 'src' doesn't.  Return false if
    both exist and 'dst' is the same age or younger than 'src'.
    """
    if not os.path.exists(dst):
        raise ValueError("file '%s' does not exist" % os.path.abspath(dst))
    if not os.path.exists(src):
        return 1

    mtime1 = os.stat(dst)[ST_MTIME]
    mtime2 = os.stat(src)[ST_MTIME]

    return mtime1 > mtime2


def all_newer(dst_files, src_files):
    """True only if all dst_files exist and are newer than all src_files."""
    return all(os.path.exists(dst) and newer(dst, src)
               for dst in dst_files for src in src_files)


def get_blas_macro_and_name(name, accelerate):
    """Complex-valued and some Accelerate functions have special symbols."""
    if accelerate:
        if name in USE_OLD_ACCELERATE:
            return '', f'{name}_'
        # Not in new Accelerate but old symbol has double underscore suffix
        elif name == 'xerbla_array':
            return '', name + '__'
    if name in WRAPPED_FUNCS:
        name = name + 'wrp'
        return 'F_FUNC', f'{name},{name.upper()}'
    return 'BLAS_FUNC', name


def write_files(file_dict):
    """
    Takes a mapping of full filepath to file contents to write at that path.
    """
    for file_path, content in file_dict.items():
        with open(file_path, 'w') as f:
            f.write(content)
