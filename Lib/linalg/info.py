"""
Linear algebra routines
=======================

 Linear Algebra Basics:

   inv        --- Find the inverse of a square matrix
   solve      --- Solve a linear system of equations
   solve_banded --- Solve a linear system of equations with a banded matrix
   det        --- Find the determinant of a square matrix
   norm       --- matrix and vector norm
   lstsq      --- Solve linear least-squares problem
   pinv       --- Pseudo-inverse (Moore-Penrose) using lstsq
   pinv2      --- Pseudo-inverse using svd

 Eigenvalues and Decompositions:

   eig        --- Find the eigenvalues and vectors of a square matrix
   eigvals    --- Find the eigenvalues of a square matrix
   lu         --- LU decomposition of a matrix
   lu_factor  --- LU decomposition returning unordered matrix and pivots
   lu_solve   --- solve Ax=b using back substitution with output of lu_factor
   svd        --- Singular value decomposition of a matrix
   svdvals    --- Singular values of a matrix
   diagsvd    --- construct matrix of singular values from output of svd
   orth       --- construct orthonormal basis for range of A using svd
   cholesky   --- Cholesky decomposition of a matrix
   cho_factor --- Cholesky decomposition for use in solving linear system
   cho_solve  --- Solve previously factored linear system
   qr         --- QR decomposition of a matrix
   schur      --- Schur decomposition of a matrix
   rsf2csf    --- Real to complex schur form
   hessenberg --- Hessenberg form of a matrix

 matrix Functions:
 
   expm       --- matrix exponential using Pade approx.
   expm2      --- matrix exponential using Eigenvalue decomp.
   expm3      --- matrix exponential using Taylor-series expansion
   logm       --- matrix logarithm
   cosm       --- matrix cosine
   sinm       --- matrix sine
   tanm       --- matrix tangent
   coshm      --- matrix hyperbolic cosine
   sinhm      --- matrix hyperbolic sine
   tanhm      --- matrix hyperbolic tangent
   signm      --- matrix sign
   sqrtm      --- matrix square root
   funm       --- Evaluating an arbitrary matrix function.

 Iterative linear systems solutions
   
   cg         --- Conjugate gradient (symmetric systems only)
   cgs        --- Conjugate gradient squared
   qmr        --- Quasi-minimal residual
   gmres      --- Generalized minimal residual
   bicg       --- Bi-conjugate gradient
   bicgstab   --- Bi-conjugate gradient stabilized
 
"""

postpone_import = 1
