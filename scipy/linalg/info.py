"""
Linear Algebra
==============

Linear Algebra Basics:

    inv:
        Find the inverse of a square matrix
    solve:
        Solve a linear system of equations
    solve_banded:
        Solve a linear system of equations with a banded matrix
    solveh_banded:
        Solve a linear system of equations with a Hermitian or symmetric
        banded matrix
    det:
        Find the determinant of a square matrix
    norm:
        matrix and vector norm
    lstsq:
        Solve linear least-squares problem
    pinv:
        Pseudo-inverse (Moore-Penrose) using lstsq
    pinv2:
        Pseudo-inverse using svd

Eigenvalue Problem:

    eig:
        Find the eigenvalues and vectors of a square matrix
    eigvals:
        Find the eigenvalues of a square matrix
    eigh:
        Find the eigenvalues and eigenvectors of a complex Hermitian or
        real symmetric matrix.
    eigvalsh:
        Find the eigenvalues of a complex Hermitian or real symmetric
        matrix.
    eig_banded:
        Find the eigenvalues and vectors of a band matrix
    eigvals_banded:
        Find the eigenvalues of a band matrix

Decompositions:

    lu:
        LU decomposition of a matrix
    lu_factor:
        LU decomposition returning unordered matrix and pivots
    lu_solve:
        solve Ax=b using back substitution with output of lu_factor
    svd:
        Singular value decomposition of a matrix
    svdvals:
        Singular values of a matrix
    diagsvd:
        construct matrix of singular values from output of svd
    orth:
        construct orthonormal basis for range of A using svd
    cholesky:
        Cholesky decomposition of a matrix
    cholesky_banded:
        Cholesky decomposition of a banded symmetric or Hermitian matrix
    cho_factor:
        Cholesky decomposition for use in solving linear system
    cho_solve:
        Solve previously factored linear system
    cho_solve_banded:
        Solve previously factored banded linear system.
    qr:
        QR decomposition of a matrix
    schur:
        Schur decomposition of a matrix
    rsf2csf:
        Real to complex schur form
    hessenberg:
        Hessenberg form of a matrix

Matrix Functions:

    expm:
        matrix exponential using Pade approx.
    expm2:
        matrix exponential using Eigenvalue decomp.
    expm3:
        matrix exponential using Taylor-series expansion
    logm:
        matrix logarithm
    cosm:
        matrix cosine
    sinm:
        matrix sine
    tanm:
        matrix tangent
    coshm:
        matrix hyperbolic cosine
    sinhm:
        matrix hyperbolic sine
    tanhm:
        matrix hyperbolic tangent
    signm:
        matrix sign
    sqrtm:
        matrix square root
    funm:
        Evaluating an arbitrary matrix function.

Special Matrices:

    block_diag:
        Construct a block diagonal matrix from submatrices.
    circulant:
        Circulant matrix
    companion:
        Companion matrix
    hadamard:
        Hadamard matrix of order 2^n
    hankel:
        Hankel matrix
    kron:
        Kronecker product of two arrays.
    leslie:
        Leslie matrix
    toeplitz:
        Toeplitz matrix
    tri:
        Construct a matrix filled with ones at and below a given diagonal. 
    tril:
        Construct a lower-triangular matrix from a given matrix.
    triu:
        Construct an upper-triangular matrix from a given matrix.
"""

postpone_import = 1
depends = ['misc','lib.lapack']
