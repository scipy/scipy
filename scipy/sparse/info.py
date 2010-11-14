"""
Sparse Matrices
===============

SciPy 2-D sparse matrix package.

Original code by Travis Oliphant.
Modified and extended by Ed Schofield, Robert Cimrman, and Nathan Bell.

There are seven available sparse matrix types:
    1. csc_matrix: Compressed Sparse Column format
    2. csr_matrix: Compressed Sparse Row format
    3. bsr_matrix: Block Sparse Row format
    4. lil_matrix: List of Lists format
    5. dok_matrix: Dictionary of Keys format
    6. coo_matrix: COOrdinate format (aka IJV, triplet format)
    7. dia_matrix: DIAgonal format

To construct a matrix efficiently, use either lil_matrix (recommended) or
dok_matrix. The lil_matrix class supports basic slicing and fancy
indexing with a similar syntax to NumPy arrays.  As illustrated below,
the COO format may also be used to efficiently construct matrices.

To perform manipulations such as multiplication or inversion, first
convert the matrix to either CSC or CSR format. The lil_matrix format is
row-based, so conversion to CSR is efficient, whereas conversion to CSC
is less so.

All conversions among the CSR, CSC, and COO formats are efficient,
linear-time operations.

Example 1
---------
Construct a 1000x1000 lil_matrix and add some values to it:

>>> from scipy.sparse import lil_matrix
>>> from scipy.sparse.linalg import spsolve
>>> from numpy.linalg import solve, norm
>>> from numpy.random import rand

>>> A = lil_matrix((1000, 1000))
>>> A[0, :100] = rand(100)
>>> A[1, 100:200] = A[0, :100]
>>> A.setdiag(rand(1000))

Now convert it to CSR format and solve A x = b for x:

>>> A = A.tocsr()
>>> b = rand(1000)
>>> x = spsolve(A, b)

Convert it to a dense matrix and solve, and check that the result
is the same:

>>> x_ = solve(A.todense(), b)

Now we can compute norm of the error with:

>>> err = norm(x-x_)
>>> err < 1e-10
True

It should be small :)


Example 2
---------

Construct a matrix in COO format:

>>> from scipy import sparse
>>> from numpy import array
>>> I = array([0,3,1,0])
>>> J = array([0,3,1,2])
>>> V = array([4,5,7,9])
>>> A = sparse.coo_matrix((V,(I,J)),shape=(4,4))

Notice that the indices do not need to be sorted.

Duplicate (i,j) entries are summed when converting to CSR or CSC.

>>> I = array([0,0,1,3,1,0,0])
>>> J = array([0,2,1,3,1,0,0])
>>> V = array([1,1,1,1,1,1,1])
>>> B = sparse.coo_matrix((V,(I,J)),shape=(4,4)).tocsr()

This is useful for constructing finite-element stiffness and mass matrices.

Further Details
---------------

CSR column indices are not necessarily sorted.  Likewise for CSC row
indices.  Use the .sorted_indices() and .sort_indices() methods when
sorted indices are required (e.g. when passing data to other libraries).

Package Contents
================

Modules
-------

.. autosummary::
   :toctree: generated/

   base - Base class for sparse matrices
   bsr - Compressed Block Sparse Row matrix format
   compressed - Sparse matrix base class using compressed storage
   construct - Functions to construct sparse matrices
   coo - A sparse matrix in COOrdinate or 'triplet' format
   csc - Compressed Sparse Column matrix format
   csgraph - Compressed Sparse graph algorithms
   csr - Compressed Sparse Row matrix format
   data - Base class for sparse matrice with a .data attribute
   dia - Sparse DIAgonal format
   dok - Dictionary Of Keys based matrix
   extract - Functions to extract parts of sparse matrices
   lil - LInked List sparse matrix class
   linalg -
   sparsetools - A collection of routines for sparse matrix operations
   spfuncs - Functions that operate on sparse matrices
   sputils - Utility functions for sparse matrix module

Classes
-------

.. autosummary::
   :toctree: generated/

   SparseEfficiencyWarning -
   SparseWarning -
   bsr_matrix - Block Sparse Row matrix
   coo_matrix - A sparse matrix in COOrdinate format
   csc_matrix - Compressed Sparse Column matrix
   csr_matrix - Compressed Sparse Row matrix
   dia_matrix - Sparse matrix with DIAgonal storage
   dok_matrix - Dictionary Of Keys based sparse matrix
   lil_matrix - Row-based linked list sparse matrix

Functions
---------

.. autosummary::
   :toctree: generated/

   bmat - Build a sparse matrix from sparse sub-blocks
   cs_graph_components -
   eye - Sparse MxN matrix whose k-th diagonal is all ones
   find -
   hstack - Stack sparse matrices horizontally (column wise)
   identity - Identity matrix in sparse format
   issparse -
   isspmatrix -
   isspmatrix_bsr -
   isspmatrix_coo -
   isspmatrix_csc -
   isspmatrix_csr -
   isspmatrix_dia -
   isspmatrix_dok -
   isspmatrix_lil -
   kron - kronecker product of two sparse matrices
   kronsum - kronecker sum of sparse matrices
   lil_diags - Generate a lil_matrix with the given diagonals
   lil_eye - RxC lil_matrix whose k-th diagonal set to one
   rand - Random values in a given shape
   spdiags - Return a sparse matrix from diagonals
   tril - Lower triangular portion of a matrix in sparse format
   triu - Upper triangular portion of a matrix in sparse format
   vstack - Stack sparse matrices vertically (row wise)

"""

__docformat__ = "restructuredtext en"

postpone_import = 1
