# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
# cython: cpow=True
# cython: language_level=3
# distutils : language = c++

__all__ = ['CythonLinewalker']

"""
Cython Linewalker
"""
import sys
import numpy as np
import cython

class CythonLinewalker:
    """ Compute banded matrix (with resource management)
    """
    def __init__(self, n):
        # Define scratch memory to avoid repeated allocations
        self.__n = n
        self.__a = np.zeros(self.__n)
        self.__b = np.zeros(self.__n)
        self.__c = np.zeros(self.__n)
        self.__d = np.zeros(self.__n)
        self.__x = np.zeros(self.__n)

    def compute_banded_matrix(self, n, e, a, d, c, f, b):
        """
        Compute banded matrix

        Source: https://www.math.uakron.edu/~kreider/anpde/_penta.f

        RESULTS: matrix has 5 bands, eadcf, with d being the main diagonal,
        e and a are the lower diagonals, and c and f are the upper diagonals.

        e is defined for rows i = 2:n-1, but is defined as e[0] to e[n-3]
        a is defined for rows i = 1:n-1, but is defined as a[0] to a[n-2]
        d is defined for rows i = 0:n-1
        c is defined for rows i = 0:n-2, but the last element isn't used
        f is defined for rows i = 0:n-3, but the last 2 elements aren't used
        b is defined for rows i = 0:n-1
        x is defined for rows i = 0:n-1, This is were the fit is stored
        """

        # Function signatures is cosnsitent with non-Cython version but 
        # we check that n is consistent across calls
        assert n == self.__n

        # Vectors a, d, c, b are temporary storage and changed in the algorithm below,
        # so we need to "deepcopy" them to avoid errors in subsequent calls to penta.
        # Vectors e and f are never changed.
        cdef double[::1] mv_a = a
        cdef double[::1] mv_b = b
        cdef double[::1] mv_c = c
        cdef double[::1] mv_d = d
        cdef double[::1] __mv_a = self.__a
        cdef double[::1] __mv_b = self.__b
        cdef double[::1] __mv_c = self.__c
        cdef double[::1] __mv_d = self.__d
        cdef double[::1] __mv_x = self.__x
        cdef double[::1] mv_e = e
        cdef double[::1] mv_f = f

        cdef size_t cn = n
        cdef size_t i

        for i in range(cn):
            __mv_a[i] = mv_a[i]
            __mv_b[i] = mv_b[i]
            __mv_c[i] = mv_c[i]
            __mv_d[i] = mv_d[i]
            __mv_x[i] = 0

        # Forward elimination
        cdef double _xmult = 0.0
        cdef size_t im1
        cdef size_t ip1
        for i in range (1, cn - 1):
            im1 = i - 1
            ip1 = i + 1
            _xmult = __mv_a[im1] / __mv_d[im1]
            __mv_d[i] -= _xmult * __mv_c[im1]
            __mv_c[i] -= _xmult * mv_f[im1]
            __mv_b[i] -= _xmult * __mv_b[im1]
            _xmult = mv_e[im1] / __mv_d[im1]
            __mv_a[i] -= _xmult * __mv_c[im1]
            __mv_d[ip1] -= _xmult * mv_f[im1]
            __mv_b[ip1] -= _xmult * __mv_b[im1]

        # Backward substitution
        cdef size_t cnm1 = cn - 1
        cdef size_t cnm2 = cn - 2
        _xmult = __mv_a[cnm2] / __mv_d[cnm2]
        __mv_d[cnm1] -= _xmult * __mv_c[cnm2]
        __mv_x[cnm1] = (__mv_b[cnm1] - _xmult * __mv_b[cnm2]) / __mv_d[cnm1]
        __mv_x[cnm2] = (__mv_b[cnm2] - __mv_c[cnm2] * __mv_x[cnm1]) / __mv_d[cnm2]
        for i in range(cn - 3, -1, -1):
            __mv_x[i] = (__mv_b[i] - mv_f[i] * __mv_x[i + 2] - __mv_c[i] * __mv_x[i + 1]) / __mv_d[i]

        return self.__x

    def compute_fit_vectors5(self, n, alpha, mu):
        """
        Compute the fit matrix using banded matrix
        n ==> Number of grid(line segment) points
        alpha ==> first - derivative smoothing parameter
        mu ==> second - derivative smoothing parameter

        Source : https://www.math.uakron.edu/~kreider/anpde/penta.f

        RESULTS:  matrix has 5 bands, eadcf, with d being the main diagonal,
        e and a are the lower diagonals, and c and f are the upper diagonals.

        e is defined for rows i = 2 : n - 1, but is defined as e[0] to e[n - 3]
        a is defined for rows i = 1 : n - 1, but is defined as a[0] to a[n - 2]
        d is defined for rows i = 0 : n - 1
        c is defined for rows i = 0 : n - 2, but the last element isn't used
        f is defined for rows i = 0 : n - 3, but the last 2 elements aren't used
        """
        # Function signatures is cosnsitent with non-Cython version but 
        # we check that n is consistent across calls
        assert n == self.__n

        e = np.zeros(n)
        a = np.zeros(n)
        d = np.zeros(n)
        c = np.zeros(n)
        f = np.zeros(n)

        cdef double[::1] mv_e = e
        cdef double[::1] mv_a = a
        cdef double[::1] mv_d = d
        cdef double[::1] mv_c = c
        cdef double[::1] mv_f = f

        cdef double calpha = alpha
        cdef double cmu = mu

        cdef size_t cn = n
        cdef size_t i

        for i in range(n):
            mv_e[i] = cmu
            mv_f[i] = cmu

        cdef double k = -(calpha + 4 * cmu)
        for i in range(n):
            mv_a[i] = k
            mv_c[i] = k

        k = -(calpha + 2 * cmu)
        mv_a[0] = k
        mv_a[cn - 2] = mv_a[0]
        mv_c[0] = k
        mv_c[cn - 2] = mv_c[0]

        k = (2 * calpha + 6 * cmu)
        for i in range(cn):
            mv_d[i] = k
        mv_d[0] = cmu
        mv_d[1] = 2 * calpha + 5 * cmu
        mv_d[cn - 2] = mv_d[1]
        mv_d[cn - 1] = mv_d[0]

        return e, a, d, c, f

    def clear(self):
        self.__n = 0
        self.__a = None
        self.__b = None
        self.__c = None
        self.__d = None
        self.__x = None

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.clear()
