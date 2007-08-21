"""
Eigenvalue solver using iterative methods.

Find k eigenvectors and eigenvalues of a matrix A using the
Arnoldi/Lanczos iterative methods from ARPACK.

These methods are most useful for large sparse matrices.

 - eigen(A,k) 
 - eigen_symmetric(A,k)

"""
global_symbols = []
postpone_import = 1
