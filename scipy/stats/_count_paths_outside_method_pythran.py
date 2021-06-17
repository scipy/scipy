import scipy.special as special
import numpy as np


#pythran export _count_paths_outside_method(int, int, int, int)
def _count_paths_outside_method(m, n, g, h):
    if m < n:
        m, n = n, m
    mg = m // g
    ng = n // g

    # Not every x needs to be considered.
    # xj holds the list of x values to be checked.
    # Wherever n*x/m + ng*h crosses an integer
    lxj = n + (mg-h)//mg
    xj = [(h + mg * j + ng-1)//ng for j in range(lxj)]
    # B is an array just holding a few values of B(x,y), the ones needed.
    # B[j] == B(x_j, j)
    if lxj == 0:
        return np.round(special.binom(m + n, n))
    B = np.zeros(lxj)
    B[0] = 1
    # Compute the B(x, y) terms
    # The binomial coefficient is an integer, but special.binom()
    # may return a float. Round it to the nearest integer.
    for j in range(1, lxj):
        Bj = np.round(special.binom(xj[j] + j, j))
        if not np.isfinite(Bj):
            raise FloatingPointError()
        for i in range(j):
            bin = np.round(special.binom(xj[j] - xj[i] + j - i, j-i))
            Bj -= bin * B[i]
        B[j] = Bj
        if not np.isfinite(Bj):
            raise FloatingPointError()
    # Compute the number of path extensions...
    num_paths = 0
    for j in range(lxj):
        bin = np.round(special.binom((m-xj[j]) + (n - j), n-j))
        term = B[j] * bin
        if not np.isfinite(term):
            raise FloatingPointError()
        num_paths += term
    return np.round(num_paths)