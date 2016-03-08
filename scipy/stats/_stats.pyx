from libc cimport math
cimport cython
cimport numpy as np

import numpy as np
import scipy.stats
from scipy.special import i0, i0e
import numpy.testing

cimport cython
from libc.math cimport cos, sin, sqrt
cimport numpy as np

cdef extern from "numpy/npy_math.h":
    double NPY_PI
    double NPY_2_PI


ctypedef fused dtype:
    np.uint8_t
    np.uint16_t
    np.uint32_t
    np.uint64_t
    np.int8_t
    np.int16_t
    np.int32_t
    np.int64_t
    np.float32_t
    np.float64_t
    np.longdouble_t


@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef double ks_2samp(dtype[:] data1, dtype[:] data2):
    cdef:
        size_t i = 0, j = 0, n1 = data1.shape[0], n2 = data2.shape[0]
        dtype d1i, d2j
        double d = 0, maxd = 0, inv_n1 = 1. / n1, inv_n2 = 1. / n2
    while i < n1 and j < n2:
        d1i = data1[i]
        d2j = data2[j]
        if d1i <= d2j:
            d += inv_n1
            i += 1
        if d1i >= d2j:
            d -= inv_n2
            j += 1
        maxd = max(maxd, math.fabs(d))
    return maxd


cdef double von_mises_cdf_series(double k,double x,unsigned int p):
    cdef double s, c, sn, cn, R, V
    cdef unsigned int n
    s = sin(x)
    c = cos(x)
    sn = sin(p*x)
    cn = cos(p*x)
    R = 0
    V = 0
    for n in range(p-1,0,-1):
        sn, cn = sn*c - cn*s, cn*c + sn*s
        R = 1./(2*n/k + R)
        V = R*(sn/n+V)

    with cython.cdivision(True):
        return 0.5 + x / (2 * NPY_PI) + V / NPY_PI

cdef von_mises_cdf_normalapprox(k, x):
    b = sqrt(NPY_2_PI) / i0e(k) # Check for negative k
    z = b*np.sin(x/2.)
    return scipy.stats.norm.cdf(z)

@cython.boundscheck(False)
def von_mises_cdf(k_obj, x_obj):
    cdef np.ndarray[double, ndim=1] temp, temp_xs, temp_ks
    cdef unsigned int i, p
    cdef double a1, a2, a3, a4, CK
    #k,x = np.broadcast_arrays(np.asarray(k),np.asarray(x))
    cdef np.ndarray k = np.asarray(k_obj)
    cdef np.ndarray x = np.asarray(x_obj)
    cdef bint zerodim = k.ndim==0 and x.ndim==0

    k = np.atleast_1d(k)
    x = np.atleast_1d(x)
    ix = np.round(x / (2 * NPY_PI))
    x = x - ix * (2 * NPY_PI)

    # These values should give 12 decimal digits
    CK=50
    a1, a2, a3, a4 = [28., 0.5, 100., 5.0]

    bx, bk = np.broadcast_arrays(x,k)
    result = np.empty(bx.shape,dtype=np.float)

    c_small_k = bk<CK
    temp = result[c_small_k]
    temp_xs = bx[c_small_k].astype(np.float)
    temp_ks = bk[c_small_k].astype(np.float)
    for i in range(len(temp)):
        p = <int>(1+a1+a2*temp_ks[i]-a3/(temp_ks[i]+a4))
        temp[i] = von_mises_cdf_series(temp_ks[i],temp_xs[i],p)
        if temp[i]<0:
            temp[i]=0
        elif temp[i]>1:
            temp[i]=1
    result[c_small_k] = temp
    result[~c_small_k] = von_mises_cdf_normalapprox(bk[~c_small_k],bx[~c_small_k])

    if not zerodim:
        return result+ix
    else:
        return (result+ix)[0]
