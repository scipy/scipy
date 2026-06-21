# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

cimport numpy as np
import numpy as np

np.import_array()


def _solve_WH_order2(const double[::1] y, double lamb):
    """Solve Whittaker-Henderson smoothing of order 2 according to Weinert."""
    cdef Py_ssize_t n = y.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] x_arr = np.empty(n, dtype=np.float64)
    cdef double[::1] x = x_arr
    # b, e and f hold the right-hand side of the forward solve and the two
    # subdiagonals of the unit lower-triangular factor L.
    cdef double[::1] b = np.empty(n, dtype=np.float64)
    cdef double[::1] e = np.empty(n, dtype=np.float64)
    cdef double[::1] f = np.empty(n, dtype=np.float64)
    cdef double d, mu, mu_old
    cdef Py_ssize_t i

    # Dividing A @ x = y by lamb puts the system in Weinert's convention,
    # (1/lamb * I + D.T @ D) @ x = 1/lamb * y, so we work with 1/lamb below.
    lamb = 1.0 / lamb

    with nogil:
        # Forward solve LD @ b = lamb * y, building the factor on the fly.
        # Indices follow the reference, shifted from Weinert's 1-based Eq. 2.2-2.6
        # to 0-based. The first and last two rows have fewer penalty neighbours
        # and so use special coefficients.
        d = 1 + lamb
        f[0] = 1 / d
        mu = 2
        e[0] = mu * f[0]
        b[0] = f[0] * lamb * y[0]
        mu_old = mu

        if n == 3:
            d = 4 + lamb - mu_old * e[0]
            mu = 2 - e[0]
        else:
            d = 5 + lamb - mu_old * e[0]
            mu = 4 - e[0]
        f[1] = 1 / d
        e[1] = mu * f[1]
        b[1] = f[1] * (lamb * y[1] + mu_old * b[0])
        mu_old = mu

        for i in range(2, n - 2):
            d = 6 + lamb - mu_old * e[i-1] - f[i-2]
            f[i] = 1 / d
            mu = 4 - e[i-1]
            e[i] = mu * f[i]
            b[i] = f[i] * (lamb * y[i] + mu_old * b[i-1] - b[i-2])
            mu_old = mu

        if n >= 4:
            i = n - 2
            d = 5 + lamb - mu_old * e[i-1] - f[i-2]
            f[i] = 1 / d
            mu = 2 - e[i-1]
            e[i] = mu * f[i]
            b[i] = f[i] * (lamb * y[i] + mu_old * b[i-1] - b[i-2])
            mu_old = mu

        i = n - 1
        d = 1 + lamb - mu_old * e[i-1] - f[i-2]
        f[i] = 1 / d
        b[i] = f[i] * (lamb * y[i] + mu_old * b[i-1] - b[i-2])

        # Back substitution L.T @ x = b.
        x[n-1] = b[n-1]
        x[n-2] = b[n-2] + e[n-2] * x[n-1]
        for i in range(n - 3, -1, -1):
            x[i] = b[i] + e[i] * x[i+1] - f[i] * x[i+2]

    return x_arr
