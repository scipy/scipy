/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*! \mainpage SuperLU Documentation
 
SuperLU is a general purpose library for the direct solution of large,
sparse, nonsymmetric systems of linear equations on high performance
machines. The library is written in C and is callable from either C or
Fortran. The library routines perform an LU decomposition with
partial pivoting and triangular system solves through forward and back
substitution.  The library also provides threshold-based ILU factorization
preconditioners.

The factorization routines can handle non-square
matrices but the triangular solves are performed only for square
matrices. The matrix columns may be preordered (before factorization)
either through library or user supplied routines. This preordering for
sparsity is completely separate from the factorization. Working
precision iterative refinement subroutines are provided for improved
backward stability. Routines are also provided to equilibrate the
system, estimate the condition number, calculate the relative backward
error, and estimate error bounds for the refined solutions. 
 
 */
