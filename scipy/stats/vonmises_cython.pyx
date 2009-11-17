import numpy as np
import scipy.stats
from scipy.special import i0
import numpy.testing
cimport numpy as np

cdef extern from "math.h":
    double cos(double theta)
    double sin(double theta)


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

    return 0.5+x/(2*np.pi) + V/np.pi

def von_mises_cdf_normalapprox(k,x,C1):
    b = np.sqrt(2/np.pi)*np.exp(k)/i0(k)
    z = b*np.sin(x/2.)
    C = 24*k
    chi = z - z**3/((C-2*z**2-16)/3.-(z**4+7/4.*z**2+167./2)/(C+C1-z**2+3))**2
    return scipy.stats.norm.cdf(z)

cimport cython
@cython.boundscheck(False)
def von_mises_cdf(k,x):
    cdef np.ndarray[double, ndim=1] temp, temp_xs, temp_ks
    cdef unsigned int i, p
    cdef double a1, a2, a3, a4, C1, CK
    #k,x = np.broadcast_arrays(np.asarray(k),np.asarray(x))
    k = np.asarray(k)
    x = np.asarray(x)
    zerodim = k.ndim==0 and x.ndim==0

    k = np.atleast_1d(k)
    x = np.atleast_1d(x)
    ix = np.round(x/(2*np.pi))
    x = x-ix*2*np.pi

    # These values should give 12 decimal digits
    CK=50
    a1, a2, a3, a4 = [28., 0.5, 100., 5.0]
    C1 = 50.1

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
    result[~c_small_k] = von_mises_cdf_normalapprox(bk[~c_small_k],bx[~c_small_k],C1)

    if not zerodim:
        return result+ix
    else:
        return (result+ix)[0]
