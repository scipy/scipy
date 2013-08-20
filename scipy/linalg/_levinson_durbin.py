
import numpy as np
import scipy.linalg



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


def main():
    check_levinson_durbin_solve()

if __name__ == '__main__':
    main()
