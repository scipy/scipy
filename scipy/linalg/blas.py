"""
Low-level BLAS functions
========================

This module contains low-level functions from the BLAS library.

.. versionadded:: 0.12.0

.. warning::

   These functions do little to no error checking.
   It is possible to cause crashes by mis-using them,
   so prefer using the higher-level routines in `scipy.linalg`.

Finding functions
=================

.. autosummary::
   :toctree: generated/

   get_blas_funcs
   find_best_blas_type

BLAS Level 1 functions
======================

.. autosummary::
   :toctree: generated/

    caxpy
    ccopy
    cdotc
    cdotu
    crotg
    cscal
    csrot
    csscal
    cswap
    dasum
    daxpy
    dcopy
    ddot
    dnrm2
    drot
    drotg
    drotm
    drotmg
    dscal
    dswap
    dzasum
    dznrm2
    icamax
    idamax
    isamax
    izamax
    sasum
    saxpy
    scasum
    scnrm2
    scopy
    sdot
    snrm2
    srot
    srotg
    srotm
    srotmg
    sscal
    sswap
    zaxpy
    zcopy
    zdotc
    zdotu
    zdrot
    zdscal
    zrotg
    zscal
    zswap

BLAS Level 2 functions
======================

.. autosummary::
   :toctree: generated/

    cgemv
    cgerc
    cgeru
    chemv
    ctrmv
    csyr
    cher
    cher2
    dgemv
    dger
    dsymv
    dtrmv
    dsyr
    dsyr2
    sgemv
    sger
    ssymv
    strmv
    ssyr
    ssyr2
    zgemv
    zgerc
    zgeru
    zhemv
    ztrmv
    zsyr
    zher
    zher2

BLAS Level 3 functions
======================

.. autosummary::
   :toctree: generated/

    cgemm
    chemm
    cherk
    cher2k
    csymm
    csyrk
    csyr2k
    dgemm
    dsymm
    dsyrk
    dsyr2k
    sgemm
    ssymm
    ssyrk
    ssyr2k
    zgemm
    zhemm
    zherk
    zher2k
    zsymm
    zsyrk
    zsyr2k

"""
#
# Author: Pearu Peterson, March 2002
#         refactoring by Fabian Pedregosa, March 2010
#

from __future__ import division, print_function, absolute_import

__all__ = ['get_blas_funcs', 'find_best_blas_type']

import numpy as _np

from scipy.linalg import _fblas
try:
    from scipy.linalg import _cblas
except ImportError:
    _cblas = None

# Expose all functions (only fblas --- cblas is an implementation detail)
empty_module = None
from scipy.linalg._fblas import *
del empty_module

# Backward compatibility
from scipy.lib._util import DeprecatedImport as _DeprecatedImport
cblas = _DeprecatedImport("scipy.linalg.blas.cblas", "scipy.linalg.blas")
fblas = _DeprecatedImport("scipy.linalg.blas.fblas", "scipy.linalg.blas")

# 'd' will be default for 'i',..
_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z', 'G':'z'}

# some convenience alias for complex functions
_blas_alias = {'cnrm2': 'scnrm2', 'znrm2': 'dznrm2',
               'cdot': 'cdotc', 'zdot': 'zdotc',
               'cger': 'cgerc', 'zger': 'zgerc',
               'sdotc': 'sdot', 'sdotu': 'sdot',
               'ddotc': 'ddot', 'ddotu': 'ddot'}


def find_best_blas_type(arrays=(), dtype=None):
    """Find best-matching BLAS/LAPACK type.

    Arrays are used to determine the optimal prefix of BLAS routines.

    Parameters
    ----------
    arrays : sequence of ndarrays, optional
        Arrays can be given to determine optimal prefix of BLAS
        routines. If not given, double-precision routines will be
        used, otherwise the most generic type in arrays will be used.
    dtype : str or dtype, optional
        Data-type specifier. Not used if `arrays` is non-empty.

    Returns
    -------
    prefix : str
        BLAS/LAPACK prefix character.
    dtype : dtype
        Inferred Numpy data type.
    prefer_fortran : bool
        Whether to prefer Fortran order routines over C order.

    """
    dtype = _np.dtype(dtype)
    prefer_fortran = False

    if arrays:
        # use the most generic type in arrays
        dtypes = [ar.dtype for ar in arrays]
        dtype = _np.find_common_type(dtypes, ())
        try:
            index = dtypes.index(dtype)
        except ValueError:
            index = 0
        if arrays[index].flags['FORTRAN']:
            # prefer Fortran for leading array with column major order
            prefer_fortran = True

    prefix = _type_conv.get(dtype.char, 'd')

    return prefix, dtype, prefer_fortran


def _get_funcs(names, arrays, dtype,
               lib_name, fmodule, cmodule,
               fmodule_name, cmodule_name, alias):
    """
    Return available BLAS/LAPACK functions.

    Used also in lapack.py. See get_blas_funcs for docstring.
    """

    funcs = []
    unpack = False
    dtype = _np.dtype(dtype)
    module1 = (cmodule, cmodule_name)
    module2 = (fmodule, fmodule_name)

    if isinstance(names, str):
        names = (names,)
        unpack = True

    prefix, dtype, prefer_fortran = find_best_blas_type(arrays, dtype)

    if prefer_fortran:
        module1, module2 = module2, module1

    for i, name in enumerate(names):
        func_name = prefix + name
        func_name = alias.get(func_name, func_name)
        func = getattr(module1[0], func_name, None)
        module_name = module1[1]
        if func is None:
            func = getattr(module2[0], func_name, None)
            module_name = module2[1]
        if func is None:
            raise ValueError(
                '%s function %s could not be found' % (lib_name, func_name))
        func.module_name, func.typecode = module_name, prefix
        func.dtype = dtype
        func.prefix = prefix  # Backward compatibility
        funcs.append(func)

    if unpack:
        return funcs[0]
    else:
        return funcs


def get_blas_funcs(names, arrays=(), dtype=None):
    """Return available BLAS function objects from names.

    Arrays are used to determine the optimal prefix of BLAS routines.

    Parameters
    ----------
    names : str or sequence of str
        Name(s) of BLAS functions without type prefix.

    arrays : sequence of ndarrays, optional
        Arrays can be given to determine optimal prefix of BLAS
        routines. If not given, double-precision routines will be
        used, otherwise the most generic type in arrays will be used.

    dtype : str or dtype, optional
        Data-type specifier. Not used if `arrays` is non-empty.


    Returns
    -------
    funcs : list
        List containing the found function(s).


    Notes
    -----
    This routine automatically chooses between Fortran/C
    interfaces. Fortran code is used whenever possible for arrays with
    column major order. In all other cases, C code is preferred.

    In BLAS, the naming convention is that all functions start with a
    type prefix, which depends on the type of the principal
    matrix. These can be one of {'s', 'd', 'c', 'z'} for the numpy
    types {float32, float64, complex64, complex128} respectively.
    The code and the dtype are stored in attributes `typecode` and `dtype`
    of the returned functions.
    """
    return _get_funcs(names, arrays, dtype,
                      "BLAS", _fblas, _cblas, "fblas", "cblas",
                      _blas_alias)
