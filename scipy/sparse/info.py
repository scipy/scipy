"""
Sparse matrix
=============

Scipy 2D sparse matrix module.

Original code by Travis Oliphant.
Modified and extended by Ed Schofield, Robert Cimrman, and Nathan Bell.

There are five available sparse matrix types:
    (1) csc_matrix: Compressed Sparse Column format
    (2) csr_matrix: Compressed Sparse Row format
    (3) lil_matrix: List of Lists format
    (4) dok_matrix: Dictionary of Keys format
    (5) coo_matrix: COOrdinate format (IJV triplets)

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

Example:
    Construct a 1000x1000 lil_matrix and add some values to it:
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

Example:
    Construct a matrix in COO format:
    >>> from scipy import sparse
    >>> from numpy import array
    >>> I = array([0,3,1,0])
    >>> J = array([0,3,1,2])
    >>> V = array([4,5,7,9])
    >>> A = sparse.coo_matrix((V,(I,J)),dims=(4,4))

    Notice that the indices do not need to be sorted.

    Duplicate (i,j) entries are summed when converting to CSR or CSC. 
    >>> I = array([0,0,1,3,1,0,0])
    >>> J = array([0,2,1,3,1,0,0])
    >>> V = array([1,1,1,1,1,1,1])
    >>> B = sparse.coo_matrix((V,(I,J)),dims=(4,4)).tocsr()
    
    This is useful for constructing finite-element stiffness and 
    mass matrices.
  
Further Details:
    CSR column indices are not necessarily sorted.  Likewise for CSC row
    indices.  Use the .ensure_sorted_indices() method when sorted indices 
    are necessary.  Note that there is no expectation for sorted indices 
    in the sparsetools module.  Furthermore some sparsetools functions 
    produce matrices with unsorted indices even when sorted input is given.
"""

postpone_import = 1
