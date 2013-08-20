
import numpy as np
from numpy.testing import assert_allclose, assert_equal
import scipy.linalg


def _levinson_durbin_linear_memory(c, r=None, y=None):
    """
    c is the first column
    r is the first row
    y is the rhs of the matrix equation
    """
    if r is None:
        r = np.conjugate(c)

    N = c.shape[0]

    # Key relating entries of the toeplitz matrix to entries of c, r,
    # assuming n is a positive integer less than N:
    # M[0, 0] == c[0]
    # M[n, :n] == c[n:0:-1]
    # M[0, 1:n+1] == r[1:n+1]

    # Initialize the forward, backward, and solution vectors.
    f = np.zeros(N)
    b = np.zeros(N)
    x = np.zeros(N)
    f_prev = np.zeros(N)
    b_prev = np.zeros(N)
    x_prev = np.zeros(N)
    f[0] = 1 / c[0]
    b[0] = 1 / c[0]
    x[0] = y[0] / c[0]

    # Compute forward, backward, and solution vectors recursively.
    for n in range(1, N):
        f, f_prev = f_prev, f
        b, b_prev = b_prev, b
        x, x_prev = x_prev, x
        f.fill(0)
        b.fill(0)
        x.fill(0)
        eps_f = np.dot(c[n:0:-1], f_prev[:n])
        eps_x = np.dot(c[n:0:-1], x_prev[:n])
        eps_b = np.dot(r[1:n+1], b_prev[:n])
        coeff = 1 / (1 - eps_f * eps_b)
        f[:n] += coeff * f_prev[:n]
        f[1:n+1] -= coeff * eps_f * b_prev[:n]
        b[1:n+1] += coeff * b_prev[:n]
        b[:n] -= coeff * eps_b * f_prev[:n]
        x[:n+1] = x_prev[:n+1] + (y[n] - eps_x) * b[:n+1]

    return x


def _levinson_durbin_linear_extra_memory(M, y):
    """
    Implement the Levinson-Durbin solver using O(n) extra memory.
    Eventually this can be done without the explicit toeplitz matrix.

    """
    N = M.shape[0]
    assert_equal(M.shape, (N, N))
    assert_equal(y.shape, (N,))

    # Initialize the forward, backward, and solution vectors.
    f = np.zeros(N)
    b = np.zeros(N)
    x = np.zeros(N)
    f_prev = np.zeros(N)
    b_prev = np.zeros(N)
    x_prev = np.zeros(N)
    f[0] = 1 / M[0, 0]
    b[0] = 1 / M[0, 0]
    x[0] = y[0] / M[0, 0]

    # Compute forward, backward, and solution vectors recursively.
    for n in range(1, N):
        f, f_prev = f_prev, f
        b, b_prev = b_prev, b
        x, x_prev = x_prev, x
        f.fill(0)
        b.fill(0)
        x.fill(0)
        eps_f = np.dot(M[n, :n], f_prev[:n])
        eps_x = np.dot(M[n, :n], x_prev[:n])
        eps_b = np.dot(M[0, 1:n+1], b_prev[:n])
        c = 1 / (1 - eps_f * eps_b)
        f[:n] += c * f_prev[:n]
        f[1:n+1] -= c * eps_f * b_prev[:n]
        b[1:n+1] += c * b_prev[:n]
        b[:n] -= c * eps_b * f_prev[:n]
        x[:n+1] = x_prev[:n+1] + (y[n] - eps_x) * b[:n+1]

    return x


def _levinson_durbin_quadratic_memory_interleaved(M, y):
    """
    Implement the Levinson-Durbin solver using O(n^2) memory,
    but interleaving the calculations of the forward, backward,
    and solution vectors as in the hypothetical O(n) memory solution.

    """
    ntotal = M.shape[0]

    # Initialize the forward and backward vectors.
    F = np.zeros_like(M)
    B = np.zeros_like(M)
    X = np.zeros_like(M)
    F[0, 0] = 1 / M[0, 0]
    B[0, 0] = 1 / M[0, 0]
    X[0, 0] = y[0] / M[0, 0]

    # Compute forward and backward vectors recursively.
    f_epsilons = []
    b_epsilons = []
    x_epsilons = []
    for n in range(1, ntotal):
        eps_f = np.dot(M[n, :n], F[n-1, :n])
        eps_x = np.dot(M[n, :n], X[n-1, :n])
        eps_b = np.dot(M[0, 1:n+1], B[n-1, :n])
        c = 1 / (1 - eps_f * eps_b)
        F[n, :n] += c * F[n-1, :n]
        F[n, 1:n+1] -= c * eps_f * B[n-1, :n]
        B[n, 1:n+1] += c * B[n-1, :n]
        B[n, :n] -= c * eps_b * F[n-1, :n]
        X[n, :n+1] = X[n-1, :n+1] + (y[n] - eps_x) * B[n, :n+1]
        f_epsilons.append(eps_f)
        b_epsilons.append(eps_b)
        x_epsilons.append(eps_x)

    return F, B, X, f_epsilons, b_epsilons, x_epsilons


def _levinson_durbin_vectors(M):
    """
    First step of Levinson-Durbin algorithm.

    Input is a square toeplitz matrix as an ndarray.
    The diagonal must be nonzero.
    Outputs are forward and backward vectors.

    """
    ntotal = M.shape[0]

    # Initialize the forward and backward vectors.
    F = np.zeros_like(M)
    B = np.zeros_like(M)
    F[0, 0] = 1 / M[0, 0]
    B[0, 0] = 1 / M[0, 0]

    # Compute forward and backward vectors recursively.
    f_epsilons = []
    b_epsilons = []
    for n in range(1, ntotal):
        eps_f = np.dot(M[n, :n], F[n-1, :n])
        eps_b = np.dot(M[0, 1:n+1], B[n-1, :n])
        c = 1 / (1 - eps_f * eps_b)
        F[n, :n] += c * F[n-1, :n]
        F[n, 1:n+1] -= c * eps_f * B[n-1, :n]
        B[n, 1:n+1] += c * B[n-1, :n]
        B[n, :n] -= c * eps_b * F[n-1, :n]
        f_epsilons.append(eps_f)
        b_epsilons.append(eps_b)

    # Return the forward and backward vectors.
    return F, B, f_epsilons, b_epsilons

def _levinson_durbin_solve(M, B, y):
    ntotal = M.shape[0]
    X = np.zeros_like(M)
    X[0, 0] = y[0] / M[0, 0]
    x_epsilons = []
    for n in range(1, ntotal):
        eps_x = np.dot(M[n, :n], X[n-1, :n])
        X[n, :n+1] = X[n-1, :n+1] + (y[n] - eps_x) * B[n, :n+1]
        x_epsilons.append(eps_x)
    return X, x_epsilons

def solve_levinson_durbin(M, y):
    F, B, f_epsilons, b_epsilons = _levinson_durbin_vectors(M)
    X, x_epsilons = _levinson_durbin_solve(M, B, y)
    return X[-1]

def check_levinson_durbin_solve():
    n = 4

    # Construct a random toeplitz matrix.
    first_column = np.random.randn(n)
    first_row = np.random.randn(n)
    first_row[0] = first_column[0]
    M = scipy.linalg.toeplitz(first_column, r=first_row)

    y = np.random.randn(n)
    x_generic = scipy.linalg.solve(M, y)
    F, B, f_epsilons, b_epsilons = _levinson_durbin_vectors(M)
    X, x_epsilons = _levinson_durbin_solve(M, B, y)
    x_levinson_durbin = X[-1]
    print M
    print y
    print 'generic solution        :', x_generic
    print 'levinson-durbin solution:', x_levinson_durbin
    print

    print 'checking forward vectors...'
    print 'forward epsilons:', f_epsilons
    for i in range(1, n):
        u = np.array(F[i-1, :i].tolist() + [0])
        T = M[:i+1, :i+1]
        print 'iteration', i
        print 'T:', T
        print 'u:', u
        print 'dot(T, u):', np.dot(T, u)
        print
    print

    print 'checking backward vectors...'
    print 'backward epsilons:', b_epsilons
    for i in range(1, n):
        u = np.array([0] + B[i-1, :i].tolist())
        T = M[:i+1, :i+1]
        print 'iteration', i
        print 'T:', T
        print 'u:', u
        print 'dot(T, u):', np.dot(T, u)
        print
    print

    print 'checking solution vectors...'
    print 'x:', x_levinson_durbin
    print 'y:', y
    print 'x epsilons:', x_epsilons
    for i in range(1, n):
        u = np.array(X[i-1, :i].tolist() + [0])
        T = M[:i+1, :i+1]
        print 'iteration', i
        print 'T:', T
        print 'u:', u
        print 'dot(T, u):', np.dot(T, u)
        print
    print

    # check that the interleaved implementation gives the same result
    info = _levinson_durbin_quadratic_memory_interleaved(M, y)
    q_F, q_B, q_X, q_f_epsilons, q_b_epsilons, q_x_epsilons = info
    assert_allclose(F, q_F)
    assert_allclose(B, q_B)
    assert_allclose(X, q_X)
    assert_allclose(f_epsilons, q_f_epsilons)
    assert_allclose(b_epsilons, q_b_epsilons)
    assert_allclose(x_epsilons, q_x_epsilons)

    x_linear_extra =  _levinson_durbin_linear_extra_memory(M, y)
    x_linear_total = _levinson_durbin_linear_memory(first_column, first_row, y)
    print 'revisit solutions...'
    print 'generic         :', x_generic
    print 'quadratic extra :', x_levinson_durbin
    print 'linear extra    :', x_linear_extra
    print 'linear total    :', x_linear_total
    print


def main():
    check_levinson_durbin_solve()

if __name__ == '__main__':
    main()
