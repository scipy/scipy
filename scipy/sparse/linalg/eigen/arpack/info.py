"""
Eigenvalue solver using iterative methods.

Find k eigenvectors and eigenvalues of a matrix A using the
Arnoldi/Lanczos iterative methods from ARPACK.

These methods are most useful for large sparse matrices.

  - eigen(A,k)
  - eigen_symmetric(A,k)

Reference
---------
 - http://www.caam.rice.edu/
software/ARPACK/
 - http://www.caam.rice.edu/software/ARPACK/UG/ug.html
 - http://books.google.com/books?hl=en&id=4E9PY7NT8a0C&dq=arpack+users+guide

"""
global_symbols = []
postpone_import = 1
