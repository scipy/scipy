from libc cimport math
cimport cython
cimport numpy as np
from numpy.math cimport PI
from numpy cimport ndarray, int64_t, intp_t

import numpy as np
import scipy.stats, scipy.special


cdef double von_mises_cdf_series(double k, double x, unsigned int p):
    cdef double s, c, sn, cn, R, V
    cdef unsigned int n
    s = math.sin(x)
    c = math.cos(x)
    sn = math.sin(p * x)
    cn = math.cos(p * x)
    R = 0
    V = 0
    for n in range(p - 1, 0, -1):
        sn, cn = sn * c - cn * s, cn * c + sn * s
        R = 1. / (2 * n / k + R)
        V = R * (sn / n + V)

    with cython.cdivision(True):
        return 0.5 + x / (2 * PI) + V / PI


cdef von_mises_cdf_normalapprox(k, x):
    b = math.sqrt(2 / PI) / scipy.special.i0e(k) # Check for negative k
    z = b * np.sin(x / 2.)
    return scipy.stats.norm.cdf(z)


@cython.boundscheck(False)
def von_mises_cdf(k_obj, x_obj):
    cdef double[:] temp, temp_xs, temp_ks
    cdef unsigned int i, p
    cdef double a1, a2, a3, a4, CK
    cdef np.ndarray k = np.asarray(k_obj)
    cdef np.ndarray x = np.asarray(x_obj)
    cdef bint zerodim = k.ndim == 0 and x.ndim == 0

    k = np.atleast_1d(k)
    x = np.atleast_1d(x)
    ix = np.round(x / (2 * PI))
    x = x - ix * (2 * PI)

    # These values should give 12 decimal digits
    CK = 50
    a1, a2, a3, a4 = 28., 0.5, 100., 5.

    bx, bk = np.broadcast_arrays(x, k)
    result = np.empty_like(bx, float)

    c_small_k = bk < CK
    temp = result[c_small_k]
    temp_xs = bx[c_small_k].astype(float)
    temp_ks = bk[c_small_k].astype(float)
    for i in range(len(temp)):
        p = <int>(1 + a1 + a2 * temp_ks[i] - a3 / (temp_ks[i] + a4))
        temp[i] = von_mises_cdf_series(temp_ks[i], temp_xs[i], p)
        temp[i] = 0 if temp[i] < 0 else 1 if temp[i] > 1 else temp[i]
    result[c_small_k] = temp
    result[~c_small_k] = von_mises_cdf_normalapprox(bk[~c_small_k], bx[~c_small_k])

    if not zerodim:
        return result + ix
    else:
        return (result + ix)[0]

@cython.wraparound(False)
@cython.boundscheck(False)
def _kendall_dis(intp_t[:] x, intp_t[:] y):
    cdef:
        intp_t sup = 1 + np.max(y)
        intp_t[::1] arr = np.zeros(sup, dtype=np.intp)
        intp_t i = 0, k = 0, size = x.size, idx
        int64_t dis = 0

    with nogil:
        while i < size:
            while k < size and x[i] == x[k]:
                dis += i
                idx = y[k]
                while idx != 0:
                    dis -= arr[idx]
                    idx = idx & (idx - 1)

                k += 1

            while i < k:
                idx = y[i]
                while idx < sup:
                    arr[idx] += 1
                    idx += idx & -idx
                i += 1

    return dis
