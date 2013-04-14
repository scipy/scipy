"""
Low-level LAPACK functions
==========================

This module contains low-level functions from the LAPACK library.

.. versionadded:: 0.12.0

.. warning::

   These functions do little to no error checking.
   It is possible to cause crashes by mis-using them,
   so prefer using the higher-level routines in `scipy.linalg`.

Finding functions
=================

.. autosummary::

   get_lapack_funcs

All functions
=============

.. autosummary::
   :toctree: generated/

   cgbsv
   cgbtrf
   cgbtrs
   cgebal
   cgees
   cgeev
   cgegv
   cgehrd
   cgelss
   cgeqp3
   cgeqrf
   cgerqf
   cgesdd
   cgesv
   cgetrf
   cgetri
   cgetrs
   cgges
   cggev
   chbevd
   chbevx
   cheev
   cheevd
   cheevr
   chegv
   chegvd
   chegvx
   claswp
   clauum
   cpbsv
   cpbtrf
   cpbtrs
   cposv
   cpotrf
   cpotri
   cpotrs
   ctrsyl
   ctrtri
   ctrtrs
   cungqr
   cungrq
   cunmqr
   dgbsv
   dgbtrf
   dgbtrs
   dgebal
   dgees
   dgeev
   dgegv
   dgehrd
   dgelss
   dgeqp3
   dgeqrf
   dgerqf
   dgesdd
   dgesv
   dgetrf
   dgetri
   dgetrs
   dgges
   dggev
   dlamch
   dlaswp
   dlauum
   dorgqr
   dorgrq
   dormqr
   dpbsv
   dpbtrf
   dpbtrs
   dposv
   dpotrf
   dpotri
   dpotrs
   dsbev
   dsbevd
   dsbevx
   dsyev
   dsyevd
   dsyevr
   dsygv
   dsygvd
   dsygvx
   dtrsyl
   dtrtri
   dtrtrs
   sgbsv
   sgbtrf
   sgbtrs
   sgebal
   sgees
   sgeev
   sgegv
   sgehrd
   sgelss
   sgeqp3
   sgeqrf
   sgerqf
   sgesdd
   sgesv
   sgetrf
   sgetri
   sgetrs
   sgges
   sggev
   slamch
   slaswp
   slauum
   sorgqr
   sorgrq
   sormqr
   spbsv
   spbtrf
   spbtrs
   sposv
   spotrf
   spotri
   spotrs
   ssbev
   ssbevd
   ssbevx
   ssyev
   ssyevd
   ssyevr
   ssygv
   ssygvd
   ssygvx
   strsyl
   strtri
   strtrs
   zgbsv
   zgbtrf
   zgbtrs
   zgebal
   zgees
   zgeev
   zgegv
   zgehrd
   zgelss
   zgeqp3
   zgeqrf
   zgerqf
   zgesdd
   zgesv
   zgetrf
   zgetri
   zgetrs
   zgges
   zggev
   zhbevd
   zhbevx
   zheev
   zheevd
   zheevr
   zhegv
   zhegvd
   zhegvx
   zlaswp
   zlauum
   zpbsv
   zpbtrf
   zpbtrs
   zposv
   zpotrf
   zpotri
   zpotrs
   ztrsyl
   ztrtri
   ztrtrs
   zungqr
   zungrq
   zunmqr

"""
#
# Author: Pearu Peterson, March 2002
#

from __future__ import division, print_function, absolute_import

__all__ = ['get_lapack_funcs']

from .blas import _get_funcs

# Backward compatibility:
from .blas import find_best_blas_type as find_best_lapack_type

from scipy.linalg import _flapack
try:
    from scipy.linalg import _clapack
except ImportError:
    _clapack = None

# Backward compatibility
from scipy.lib._util import DeprecatedImport as _DeprecatedImport
clapack = _DeprecatedImport("scipy.linalg.blas.clapack", "scipy.linalg.lapack")
flapack = _DeprecatedImport("scipy.linalg.blas.flapack", "scipy.linalg.lapack")

# Expose all functions (only flapack --- clapack is an implementation detail)
empty_module = None
from scipy.linalg._flapack import *
del empty_module

# some convenience alias for complex functions
_lapack_alias = {
    'corgqr': 'cungqr', 'zorgqr': 'zungqr',
    'cormqr': 'cunmqr', 'zormqr': 'zunmqr',
    'corgrq': 'cungrq', 'zorgrq': 'zungrq',
}


def get_lapack_funcs(names, arrays=(), dtype=None):
    """Return available LAPACK function objects from names.

    Arrays are used to determine the optimal prefix of LAPACK routines.

    Parameters
    ----------
    names : str or sequence of str
        Name(s) of LAPACK functions withouth type prefix.

    arrays : sequency of ndarrays, optional
        Arrays can be given to determine optiomal prefix of LAPACK
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
    This routines automatically chooses between Fortran/C
    interfaces. Fortran code is used whenever possible for arrays with
    column major order. In all other cases, C code is preferred.

    In LAPACK, the naming convention is that all functions start with a
    type prefix, which depends on the type of the principal
    matrix. These can be one of {'s', 'd', 'c', 'z'} for the numpy
    types {float32, float64, complex64, complex128} respectevely, and
    are stored in attribute `typecode` of the returned functions.
    """
    return _get_funcs(names, arrays, dtype,
                      "LAPACK", _flapack, _clapack,
                      "flapack", "clapack", _lapack_alias)
