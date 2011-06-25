"""
==================================================
Sparse linear algebra (:mod:`scipy.sparse.linalg`)
==================================================

.. currentmodule:: scipy.sparse.linalg

Abstract linear operators
-------------------------

.. autosummary::
   :toctree: generated/

   LinearOperator -- abstract representation of a linear operator
   aslinearoperator -- convert an object to an abstract linear operator

Solving linear problems
-----------------------

Direct methods for linear equation systems:

.. autosummary::
   :toctree: generated/

   spsolve -- Solve the sparse linear system Ax=b
   factorized -- Pre-factorize matrix to a function solving a linear system

Iterative methods for linear equation systems:

.. autosummary::
   :toctree: generated/

   bicg -- Use BIConjugate Gradient iteration to solve A x = b
   bicgstab -- Use BIConjugate Gradient STABilized iteration to solve A x = b
   cg -- Use Conjugate Gradient iteration to solve A x = b
   cgs -- Use Conjugate Gradient Squared iteration to solve A x = b
   gmres -- Use Generalized Minimal RESidual iteration to solve A x = b
   lgmres -- Solve a matrix equation using the LGMRES algorithm
   minres -- Use MINimum RESidual iteration to solve Ax = b
   qmr -- Use Quasi-Minimal Residual iteration to solve A x = b

Iterative methods for least-squares problems:

.. autosummary::
   :toctree: generated/

   lsqr -- Find the least-squares solution to a sparse linear equation system

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

Complete or incomplete LU factorizations

.. autosummary::
   :toctree: generated/

   splu -- Compute a LU decomposition for a sparse matrix
   spilu -- Compute an incomplete LU decomposition for a sparse matrix


Exceptions
----------

.. autosummary::
   :toctree: generated/

   ArpackNoConvergence
   ArpackError

"""

__docformat__ = "restructuredtext en"

postpone_import = 1
