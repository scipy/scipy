"""
Linear Solvers
==============

The default solver is SuperLU (included in the scipy distribution),
which can solve real or complex linear systems in both single and
double precisions.  It is automatically replaced by UMFPACK, if
available.  Note that UMFPACK works in double precision only, so
switch it off by using::

    >>> from scipy.sparse.linalg import spsolve
    >>> spsolve(..., use_umfpack=False)

to solve in the single precision.

Example session::

    >>> from scipy.sparse import csc_array, dia_array
    >>> from numpy import array
    >>>
    >>> print("Inverting a sparse linear system:")
    >>> print("The sparse matrix (constructed from diagonals):")
    >>> a = dia_array(([[1, 2, 3, 4, 5], [6, 5, 8, 9, 10]], [0, 1]), shape=(5, 5))
    >>> b = array([1, 2, 3, 4, 5])
    >>> print("Solve: single precision complex:")
    >>> a = a.astype('F')
    >>> x = spsolve(a, b, use_umfpack=False)
    >>> print(x)
    >>> print("Error: ", a@x-b)
    >>>
    >>> print("Solve: double precision complex:")
    >>> a = a.astype('D')
    >>> x = spsolve(a, b, use_umfpack=True)
    >>> print(x)
    >>> print("Error: ", a@x-b)
    >>>
    >>> print("Solve: double precision:")
    >>> a = a.astype('d')
    >>> x = spsolve(a, b)
    >>> print(x)
    >>> print("Error: ", a@x-b)
    >>>
    >>> print("Solve: single precision:")
    >>> a = a.astype('f')
    >>> x = spsolve(a, b.astype('f'), use_umfpack=False)
    >>> print(x)
    >>> print("Error: ", a@x-b)

"""

#import umfpack
#__doc__ = '\n\n'.join( (__doc__,  umfpack.__doc__) )
#del umfpack

from .linsolve import *
from ._superlu import SuperLU
from . import _add_newdocs
from . import linsolve

__all__ = [
    'SuperLU', 'factorized',
    'spilu', 'splu', 'spsolve', 'is_sptriangular',
    'spsolve_triangular', 'spbandwidth',
]

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
