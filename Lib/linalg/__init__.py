""" Linear algebra routines.

 Linear Algebra Basics:

   inv        --- Find the inverse of a square matrix
   solve      --- Solve a linear system of equations
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
   cholesky   --- Cholesky decomposition of a matrix
   qr         --- QR decomposition of a matrix
   schur      --- Schur decomposition of a matrix
   rsf2csf    --- Real to complex schur form

 matrix Functions:
 
   expm       --- matrix exponential using Pade approx.
   expm2      --- matrix exponential using Eigenvalue decomp.
   expm3      --- matrix exponential using Taylor-series expansion
   cosm       --- matrix cosine
   sinm       --- matrix sine
   tanm       --- matrix tangent
   coshm      --- matrix hyperbolic cosine
   sinhm      --- matrix hyperbolic sine
   tanhm      --- matrix hyperbolic tangent
   funm       --- Evaluating an arbitrary matrix function.
 
"""

from linalg_version import linalg_version as __version__

from basic import *
from decomp import *
from matfuncs import *
from blas import *

################## test functions #########################

def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level))
    return runner

def test_suite(level=1):
    import scipy_test.testing
    exec 'import %s as this_mod' % (__name__)
    return scipy_test.testing.harvest_test_suites(this_mod,level=level)
