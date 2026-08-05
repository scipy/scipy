"""
Sparse linear algebra (:mod:`scipy.sparse.linalg`)
==================================================

.. currentmodule:: scipy.sparse.linalg

Abstract linear operators
-------------------------

.. autosummary::
   :toctree: generated/

   LinearOperator -- abstract representation of a linear operator
   aslinearoperator -- convert an object to an abstract linear operator

Matrix Operations
-----------------

.. autosummary::
   :toctree: generated/

   inv -- compute the sparse matrix inverse
   expm -- compute the sparse matrix exponential
   expm_multiply -- compute the product of a matrix exponential and a matrix
   funm_multiply_krylov -- use a Krylov method to compute f(A)b for a general f
   matrix_power -- compute the matrix power by raising a matrix to an exponent

Matrix norms
------------

.. autosummary::
   :toctree: generated/

   norm -- Norm of a sparse matrix
   onenormest -- Estimate the 1-norm of a sparse matrix

Solving linear problems
-----------------------

Direct methods for linear equation systems:

.. autosummary::
   :toctree: generated/

   spsolve -- Solve the sparse linear system Ax=b
   spsolve_triangular -- Solve sparse linear system Ax=b for a triangular A.
   is_sptriangular -- Check if sparse A is triangular.
   spbandwidth -- Find the bandwidth of a sparse matrix.
   factorized -- Pre-factorize matrix to a function solving a linear system
   MatrixRankWarning -- Warning on exactly singular matrices
   use_solver -- Select direct solver to use

Iterative methods for linear equation systems:

.. autosummary::
   :toctree: generated/

   bicg -- Use BIConjugate Gradient iteration to solve Ax = b
   bicgstab -- Use BIConjugate Gradient STABilized iteration to solve Ax = b
   cg -- Use Conjugate Gradient iteration to solve Ax = b
   cgs -- Use Conjugate Gradient Squared iteration to solve Ax = b
   gmres -- Use Generalized Minimal RESidual iteration to solve Ax = b
   lgmres -- Solve a matrix equation using the LGMRES algorithm
   minres -- Use MINimum RESidual iteration to solve Ax = b
   qmr -- Use Quasi-Minimal Residual iteration to solve Ax = b
   gcrotmk -- Solve a matrix equation using the GCROT(m,k) algorithm
   tfqmr -- Use Transpose-Free Quasi-Minimal Residual iteration to solve Ax = b

Iterative methods for least-squares problems:

.. autosummary::
   :toctree: generated/

   lsqr -- Find the least-squares solution to a sparse linear equation system
   lsmr -- Find the least-squares solution to a sparse linear equation system

Matrix factorizations
---------------------

Eigenvalue problems:

.. autosummary::
   :toctree: generated/

   eigs -- Find k eigenvalues and eigenvectors of the square matrix A
   eigsh -- Find k eigenvalues and eigenvectors of a symmetric matrix
   lobpcg -- Solve symmetric partial eigenproblems with optional preconditioning

Singular values problems:

.. autosummary::
   :toctree: generated/

   svds -- Compute k singular values/vectors for a sparse matrix

The `svds` function supports the following solvers:

.. toctree::

    sparse.linalg.svds-arpack
    sparse.linalg.svds-lobpcg
    sparse.linalg.svds-propack

Complete or incomplete LU factorizations

.. autosummary::
   :toctree: generated/

   splu -- Compute a LU decomposition for a sparse matrix
   spilu -- Compute an incomplete LU decomposition for a sparse matrix
   SuperLU -- Object representing an LU factorization

Sparse arrays with structure
----------------------------

.. autosummary::
   :toctree: generated/

   LaplacianNd -- Laplacian on a uniform rectangular grid in ``N`` dimensions

Exceptions
----------

.. autosummary::
   :toctree: generated/

   ArpackNoConvergence
   ArpackError

"""

from . import _isolve
from ._isolve import *
from . import _dsolve
from ._dsolve import *
from . import _interface
from ._interface import *
from . import _eigen
from ._eigen import *
from . import _matfuncs
from ._matfuncs import *
from . import _onenormest
from ._onenormest import *
from . import _norm
from ._norm import *
from . import _expm_multiply
from ._expm_multiply import *
from . import _funm_multiply_krylov
from ._funm_multiply_krylov import *
from . import _special_sparse_arrays
from ._special_sparse_arrays import *

# Deprecated namespaces, to be removed in v2.0.0
from . import isolve, dsolve, interface, eigen, matfuncs

__all__ = []
__all__ += _isolve.__all__
__all__ += _dsolve.__all__
__all__ += _interface.__all__
__all__ += _eigen.__all__
__all__ += _matfuncs.__all__
__all__ += _onenormest.__all__
__all__ += _norm.__all__
__all__ += _expm_multiply.__all__
__all__ += _funm_multiply_krylov.__all__
__all__ += _special_sparse_arrays.__all__
__all__ += ['isolve', 'dsolve', 'interface', 'eigen', 'matfuncs']

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
