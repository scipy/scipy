"""
Linear Algebra
==============

Basics
------

.. autosummary::
   :toctree: generated/

   inv - Find the inverse of a square matrix
   solve - Solve a linear system of equations
   solve_banded - Solve a banded linear system
   solveh_banded - Solve a Hermitian or symmetric banded system
   det - Find the determinant of a square matrix
   norm - Matrix and vector norm
   lstsq - Solve a linear least-squares problem
   pinv - Pseudo-inverse (Moore-Penrose) using lstsq
   pinv2 - Pseudo-inverse using svd

Eigenvalue Problems
-------------------

.. autosummary::
   :toctree: generated/

   eig - Find the eigenvalues and eigenvectors of a square matrix
   eigvals - Find just the eigenvalues of a square matrix
   eigh - Find the e-vals and e-vectors of a Hermitian or symmetric matrix
   eigvalsh - Find just the eigenvalues of a Hermitian or symmetric matrix
   eig_banded - Find the eigenvalues and eigenvectors of a banded matrix
   eigvals_banded - Find just the eigenvalues of a banded matrix

Decompositions
--------------

.. autosummary::
   :toctree: generated/

   lu - LU decomposition of a matrix
   lu_factor - LU decomposition returning unordered matrix and pivots
   lu_solve - Solve Ax=b using back substitution with output of lu_factor
   svd - Singular value decomposition of a matrix
   svdvals - Singular values of a matrix
   diagsvd - Construct matrix of singular values from output of svd
   orth - Construct orthonormal basis for the range of A using svd
   cholesky - Cholesky decomposition of a matrix
   cholesky_banded - Cholesky decomp. of a sym. or Hermitian banded matrix
   cho_factor - Cholesky decomposition for use in solving a linear system
   cho_solve - Solve previously factored linear system
   cho_solve_banded - Solve previously factored banded linear system
   qr - QR decomposition of a matrix
   schur - Schur decomposition of a matrix
   rsf2csf - Real to complex Schur form
   hessenberg - Hessenberg form of a matrix

Matrix Functions
----------------

.. autosummary::
   :toctree: generated/

   expm - Matrix exponential using Pade approximation
   expm2 - Matrix exponential using eigenvalue decomposition
   expm3 - Matrix exponential using Taylor-series expansion
   logm - Matrix logarithm
   cosm - Matrix cosine
   sinm - Matrix sine
   tanm - Matrix tangent
   coshm - Matrix hyperbolic cosine
   sinhm - Matrix hyperbolic sine
   tanhm - Matrix hyperbolic tangent
   signm - Matrix sign
   sqrtm - Matrix square root
   funm - Evaluating an arbitrary matrix function

Special Matrices
----------------

.. autosummary::
   :toctree: generated/

   block_diag - Construct a block diagonal matrix from submatrices
   circulant - Circulant matrix
   companion - Companion matrix
   hadamard - Hadamard matrix of order 2**n
   hankel - Hankel matrix
   kron - Kronecker product of two arrays
   leslie - Leslie matrix
   toeplitz - Toeplitz matrix
   tri - Construct a matrix filled with ones at and below a given diagonal
   tril - Construct a lower-triangular matrix from a given matrix
   triu - Construct an upper-triangular matrix from a given matrix

"""

postpone_import = 1
depends = ['misc','lib.lapack']
