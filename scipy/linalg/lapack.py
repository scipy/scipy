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


   sgbsv
   dgbsv
   cgbsv
   zgbsv

   sgbtrf
   dgbtrf
   cgbtrf
   zgbtrf

   sgbtrs
   dgbtrs
   cgbtrs
   zgbtrs

   sgebal
   dgebal
   cgebal
   zgebal

   sgees
   dgees
   cgees
   zgees

   sgeev
   dgeev
   cgeev
   zgeev

   sgeev_lwork
   dgeev_lwork
   cgeev_lwork
   zgeev_lwork

   sgegv
   dgegv
   cgegv
   zgegv

   sgehrd
   dgehrd
   cgehrd
   zgehrd

   sgehrd_lwork
   dgehrd_lwork
   cgehrd_lwork
   zgehrd_lwork

   sgelss
   dgelss
   cgelss
   zgelss

   sgeqp3
   dgeqp3
   cgeqp3
   zgeqp3

   sgeqrf
   dgeqrf
   cgeqrf
   zgeqrf

   sgerqf
   dgerqf
   cgerqf
   zgerqf

   sgesdd
   dgesdd
   cgesdd
   zgesdd

   sgesdd_lwork
   dgesdd_lwork
   cgesdd_lwork
   zgesdd_lwork

   sgesv
   dgesv
   cgesv
   zgesv

   sgetrf
   dgetrf
   cgetrf
   zgetrf

   sgetri
   dgetri
   cgetri
   zgetri

   sgetri_lwork
   dgetri_lwork
   cgetri_lwork
   zgetri_lwork

   sgetrs
   dgetrs
   cgetrs
   zgetrs

   sgges
   dgges
   cgges
   zgges

   sggev
   dggev
   cggev
   zggev

   chbevd
   zhbevd

   chbevx
   zhbevx

   cheev
   zheev

   cheevd
   zheevd

   cheevr
   zheevr

   chegv
   zhegv

   chegvd
   zhegvd

   chegvx
   zhegvx

   slaswp
   dlaswp
   claswp
   zlaswp

   slauum
   dlauum
   clauum
   zlauum

   spbsv
   dpbsv
   cpbsv
   zpbsv

   spbtrf
   dpbtrf
   cpbtrf
   zpbtrf

   spbtrs
   dpbtrs
   cpbtrs
   zpbtrs

   sposv
   dposv
   cposv
   zposv

   spotrf
   dpotrf
   cpotrf
   zpotrf

   spotri
   dpotri
   cpotri
   zpotri

   spotrs
   dpotrs
   cpotrs
   zpotrs

   strsyl
   dtrsyl
   ctrsyl
   ztrsyl

   strtri
   dtrtri
   ctrtri
   ztrtri

   strtrs
   dtrtrs
   ctrtrs
   ztrtrs

   cunghr
   zunghr

   cungqr
   zungqr

   cungrq
   zungrq

   cunmqr
   zunmqr

   sgtsv
   dgtsv
   cgtsv
   zgtsv

   sptsv
   dptsv
   cptsv
   zptsv

   slamch
   dlamch

   sorghr
   dorghr

   sorgqr
   dorgqr

   sorgrq
   dorgrq

   sormqr
   dormqr

   ssbev
   dsbev

   ssbevd
   dsbevd

   ssbevx
   dsbevx

   ssyev
   dsyev

   ssyevd
   dsyevd

   ssyevr
   dsyevr

   ssygv
   dsygv

   ssygvd
   dsygvd

   ssygvx
   dsygvx

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
from scipy._lib._util import DeprecatedImport as _DeprecatedImport
clapack = _DeprecatedImport("scipy.linalg.blas.clapack", "scipy.linalg.lapack")
flapack = _DeprecatedImport("scipy.linalg.blas.flapack", "scipy.linalg.lapack")

# Expose all functions (only flapack --- clapack is an implementation detail)
empty_module = None
from scipy.linalg._flapack import *
del empty_module

# some convenience alias for complex functions
_lapack_alias = {
    'corghr': 'cunghr', 'zorghr': 'zunghr',
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
        Name(s) of LAPACK functions without type prefix.

    arrays : sequence of ndarrays, optional
        Arrays can be given to determine optimal prefix of LAPACK
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

    In LAPACK, the naming convention is that all functions start with a
    type prefix, which depends on the type of the principal
    matrix. These can be one of {'s', 'd', 'c', 'z'} for the numpy
    types {float32, float64, complex64, complex128} respectevely, and
    are stored in attribute `typecode` of the returned functions.
    """
    return _get_funcs(names, arrays, dtype,
                      "LAPACK", _flapack, _clapack,
                      "flapack", "clapack", _lapack_alias)
