"""
Sparse matrix
=============

Scipy 2D sparse matrix module.

Original code by Travis Oliphant.
Modified and extended by Ed Schofield and Robert Cimrman.

There are four available sparse matrix types:
    (1) csc_matrix: Compressed Sparse Column format
    (2) csr_matrix: Compressed Sparse Row format
    (3) lil_matrix: List of Lists format
    (4) dok_matrix: Dictionary of Keys format

To construct a matrix efficiently, use either lil_matrix (recommended) or
dok_matrix. The lil_matrix class supports basic slicing and fancy
indexing with a similar syntax to NumPy arrays.

To perform manipulations such as multiplication or inversion, first
convert the matrix to either CSC or CSR format. The lil_matrix format is
row-based, so conversion to CSR is efficient, whereas conversion to CSC
is less so.

Example:
    Construct a 10x1000 lil_matrix and add some values to it:
    >>> from scipy import sparse, linsolve
    >>> from numpy import linalg
    >>> from numpy.random import rand
    >>> A = sparse.lil_matrix((1000, 1000))
    >>> A[0, :100] = rand(100)
    >>> A[1, 100:200] = A[0, :100]
    >>> A.setdiag(rand(1000))

    Now convert it to CSR format and solve (A A^T) x = b for x:
    >>> A = A.tocsr()
    >>> b = rand(1000)
    >>> x = linsolve.spsolve(A * A.T, b)

    Convert it to a dense matrix and solve, and check that the result
    is the same:
    >>> A_ = A.todense()
    >>> x_ = linalg.solve(A_ * A_.T, b)
    >>> err = linalg.norm(x-x_)

    Now we can print the error norm with:
        print "Norm error =", err
    It should be small :)

"""

postpone_import = 1
