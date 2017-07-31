"""Basic linear factorizations needed by the solver."""

from __future__ import division, print_function, absolute_import
from scipy.sparse import (linalg, bmat, csc_matrix, eye, issparse)
from scipy.sparse.linalg import LinearOperator
import scipy.linalg
try:
    from sksparse.cholmod import cholesky_AAt
except ImportError as e:
    pass
import numpy as np

__all__ = [
    'orthogonality',
    'projections',
]


def orthogonality(A, g):
    """Measure orthogonality between a vector and the null space of a matrix.

    Compute a measure of orthogonality between the null space
    of the (possibly sparse) matrix ``A`` and a given vector ``g``.

    The formula is a simplified (and cheaper) version of formula (3.13)
    from [1]_.
    ``orth =  norm(A g, ord=2)/(norm(A, ord='fro')*norm(g, ord=2))``.

    References
    ----------
    .. [1] Gould, Nicholas IM, Mary E. Hribar, and Jorge Nocedal.
           "On the solution of equality constrained quadratic
            programming problems arising in optimization."
            SIAM Journal on Scientific Computing 23.4 (2001): 1376-1395.
    """
    # Compute vector norms
    norm_g = np.linalg.norm(g)
    # Compute Frobenius norm of the matrix A
    if issparse(A):
        norm_A = linalg.norm(A, ord='fro')
    else:
        norm_A = np.linalg.norm(A, ord='fro')

    # Check if norms are zero
    if norm_g == 0 or norm_A == 0:
        return 0

    norm_A_g = np.linalg.norm(A.dot(g))
    # Orthogonality measure
    orth = norm_A_g / (norm_A*norm_g)
    return orth


def projections(A, method=None, orth_tol=1e-12, max_refin=3):
    """Return three linear operators related with a given matrix A.

    Parameters
    ----------
    A : sparse matrix (or ndarray), shape (m, n)
        Matrix ``A`` used in the projection.
    method : string, optional
        Method used for compute the given linear
        operators. Should be one of:

            - 'NormalEquation': The operators
               will be computed using the
               so-called normal equation approach
               explained in [1]_. In order to do
               so the Cholesky factorization of
               ``(A A.T)`` is computed. Exclusive
               for sparse matrices.
            - 'AugmentedSystem': The operators
               will be computed using the
               so-called augmented system approach
               explained in [1]_. Exclusive
               for sparse matrices.
            - 'QRFactorization': Compute projections
               using QR factorization. Exclusive for
               dense matrices.

    orth_tol : float, optional
        Tolerance for iterative refinements.
    max_refin : int, optional
        Maximum number of iterative refinements

    Returns
    -------
    Z : LinearOperator, shape (n, n)
        Null-space operator. For a given vector ``x``,
        the null space operator is equivalent to apply
        a projection matrix ``P = I - A.T inv(A A.T) A``
        to the vector. It can be shown that this is
        equivalent to project ``x`` into the null space
        of A.
    LS : LinearOperator, shape (m, n)
        Least-Square operator. For a given vector ``x``,
        the least-square operator is equivalent to apply a
        pseudoinverse matrix ``pinv(A.T) = inv(A A.T) A``
        to the vector. It can be shown that this vector
        ``pinv(A.T) x`` is the least_square solution to
        ``A.T y = x``.
    Y : LinearOperator, shape (n, m)
        Row-space operator. For a given vector ``x``,
        the row-space operator is equivalent to apply a
        projection matrix ``Q = A.T inv(A A.T)``
        to the vector.  It can be shown that this
        vector ``y = Q x``  the minimum norm solution
        of ``A y = x``.

    Notes
    -----
    Uses iterative refinements described in [1]
    during the computation of ``Z`` in order to
    cope with the possibility of large roundoff errors.

    References
    ----------
    .. [1] Gould, Nicholas IM, Mary E. Hribar, and Jorge Nocedal.
        "On the solution of equality constrained quadratic
        programming problems arising in optimization."
        SIAM Journal on Scientific Computing 23.4 (2001): 1376-1395.
    """
    m, n = np.shape(A)

    # Check Argument
    if issparse(A):
        if method is None:
            method = "AugmentedSystem"

        if method not in ("NormalEquation", "AugmentedSystem"):
            raise ValueError("Method not allowed for the given matrix.")

    else:
        if method is None:
            method = "QRFactorization"

        if method != "QRFactorization":
            raise ValueError("Method not allowed for the given matrix.")

    if method == 'NormalEquation':
        # Cholesky factorization
        factor = cholesky_AAt(A)

        # z = x - A.T inv(A A.T) A x
        def null_space(x):
            v = factor(A.dot(x))
            z = x - A.T.dot(v)

            # Iterative refinement to improve roundoff
            # errors described in [2]_, algorithm 5.1.
            k = 0
            while orthogonality(A, z) > orth_tol:
                if k >= max_refin:
                    break
                # z_next = z - A.T inv(A A.T) A z
                v = factor(A.dot(z))
                z = z - A.T.dot(v)
                k += 1

            return z

        # z = inv(A A.T) A x
        def least_squares(x):
            return factor(A.dot(x))

        # z = A.T inv(A A.T) x
        def row_space(x):
            return A.T.dot(factor(x))

    elif method == 'AugmentedSystem':
        # Form augmented system
        K = csc_matrix(bmat([[eye(n), A.T], [A, None]]))
        # LU factorization
        # TODO: Use a symmetric indefinite factorization
        #       to solve the system twice as fast (because
        #       of the symmetry).
        factor = linalg.splu(K)

        # z = x - A.T inv(A A.T) A x
        # is computed solving the extended system:
        # [I A.T] * [ z ] = [x]
        # [A  O ]   [aux]   [0]
        def null_space(x):
            # v = [x]
            #     [0]
            v = np.hstack([x, np.zeros(m)])
            # lu_sol = [ z ]
            #          [aux]
            lu_sol = factor.solve(v)
            z = lu_sol[:n]

            # Iterative refinement to improve roundoff
            # errors described in [2]_, algorithm 5.2.
            k = 0
            while orthogonality(A, z) > orth_tol:
                if k >= max_refin:
                    break
                # new_v = [x] - [I A.T] * [ z ]
                #         [0]   [A  O ]   [aux]
                new_v = v - K.dot(lu_sol)
                # [I A.T] * [delta  z ] = new_v
                # [A  O ]   [delta aux]
                lu_update = factor.solve(new_v)
                #  [ z ] += [delta  z ]
                #  [aux]    [delta aux]
                lu_sol += lu_update
                z = lu_sol[:n]
                k += 1

            # return z = x - A.T inv(A A.T) A x
            return z

        # z = inv(A A.T) A x
        # is computed solving the extended system:
        # [I A.T] * [aux] = [x]
        # [A  O ]   [ z ]   [0]
        def least_squares(x):
            # v = [x]
            #     [0]
            v = np.hstack([x, np.zeros(m)])
            # lu_sol = [aux]
            #          [ z ]
            lu_sol = factor.solve(v)
            # return z = inv(A A.T) A x
            return lu_sol[n:m+n]

        # z = A.T inv(A A.T) x
        # is computed solving the extended system:
        # [I A.T] * [ z ] = [0]
        # [A  O ]   [aux]   [x]
        def row_space(x):
            # v = [0]
            #     [x]
            v = np.hstack([np.zeros(n), x])
            # lu_sol = [ z ]
            #          [aux]
            lu_sol = factor.solve(v)
            # return z = A.T inv(A A.T) x
            return lu_sol[:n]

    elif method == "QRFactorization":
        # QRFactorization
        Q, R, P = scipy.linalg.qr(A.T, pivoting=True,  mode='economic')

        # z = x - A.T inv(A A.T) A x
        def null_space(x):
            # v = P inv(R) Q.T x
            aux1 = Q.T.dot(x)
            aux2 = scipy.linalg.solve_triangular(R, aux1, lower=False)
            v = np.zeros(m)
            v[P] = aux2
            z = x - A.T.dot(v)

            # Iterative refinement to improve roundoff
            # errors described in [2]_, algorithm 5.1.
            k = 0
            while orthogonality(A, z) > orth_tol:
                if k >= max_refin:
                    break
                # v = P inv(R) Q.T x
                aux1 = Q.T.dot(z)
                aux2 = scipy.linalg.solve_triangular(R, aux1, lower=False)
                v[P] = aux2
                # z_next = z - A.T v
                z = z - A.T.dot(v)
                k += 1

            return z

        # z = inv(A A.T) A x
        def least_squares(x):
            # z = P inv(R) Q.T x
            aux1 = Q.T.dot(x)
            aux2 = scipy.linalg.solve_triangular(R, aux1, lower=False)
            z = np.zeros(m)
            z[P] = aux2
            return z

        # z = A.T inv(A A.T) x
        def row_space(x):
            # z = Q inv(R.T) P.T x
            aux1 = x[P]
            aux2 = scipy.linalg.solve_triangular(R, aux1,
                                                 lower=False,
                                                 trans='T')
            z = Q.dot(aux2)
            return z

    Z = LinearOperator((n, n), null_space)
    LS = LinearOperator((m, n), least_squares)
    Y = LinearOperator((n, m), row_space)

    return Z, LS, Y
