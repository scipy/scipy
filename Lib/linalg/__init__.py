""" Linear algebra routines.

 Solving Linear Systems:

   inv        --- Find the inverse of a square matrix
   solve      --- Solve a linear system of equations
   det        --- determinant of a matrix
   pinv       --- Moore-Penrose pseudo inverse (using least-squares)
   pinv2      --- Moore-Penrose pseudo inverse (using SVD)
   lstsq      --- Least-squares solve

 matrix Factorizations:
 
   lu         --- LU decomposition
   cholesky   --- Cholesky factorization
   qr         --- QR factorization
   schur      --- Schur decomposition
   rsf2csf    --- Real to complex schur form.
   norm       --- vector and matrix norm
   eig        --- eigenvectors and eigenvalues
   eigvals    --- only eigenvalues
   svd        --- singular value decomposition

 matrix Functions
 
   expm       --- exponential (using Pade approximation)
   cosm       --- cosine
   sinm       --- sine
   tanm       --- tangent
   coshm      --- hyperbolic cosine
   sinhm      --- hyperbolic sine
   tanhm      --- hyperbolic tangent
   funm       --- arbitrary function 
"""

import fblas, flapack, cblas, clapack
from linear_algebra import *
