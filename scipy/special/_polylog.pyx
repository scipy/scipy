from __future__ import division
import cython
import warnings
from special import bernoulli,gamma,digamma
import numpy as np
from math import pi

cdef extern from "math.h":
    double exp(double x) nogil
    double log(double x) nogil
    double fabs(double x) nogil

cdef extern from "numpy/npy_math.h":
    ctypedef struct npy_cdouble:
        double real
        double imag

    double npy_cabs(npy_cdouble z) nogil
    npy_cdouble npy_clog(npy_cdouble z) nogil
    npy_cdouble npy_cexp(npy_cdouble z) nogil
    int npy_isnan(double x) nogil
    double NPY_INFINITY
    double NPY_PI

cdef extern from "cephes.h":
     #double zeta(double x, double q)
     double zetac(double x)

EPS=1.e-16
def fac(n) : return gamma(n+1)

def zetan(s): return 1+zetac(s)

def altzeta(s):
    if s == 1:
        return np.log(2) + 0*s
    return (1-pow(2, 1-s)) * zetan(s)


def harmonic(x):
    if x == 0 or x == 1:
        return x
    return digamma(x+1) -digamma(1)

cpdef double complex ln(double complex x):
    cdef npy_cdouble r
    r = npy_clog((<npy_cdouble*>&x)[0])
    return (<double complex*>&r)[0]

cpdef double complex pyln(double complex z):
    if np.isreal(z) and z.real>0 : return np.log(z)
    if np.isreal(z) and z.real<0 : return complex(np.log(-z),np.pi)
    return np.log(complex(z))

def bernpoly(n, z):
    # Slow implementation:
    #return sum(binomial(n,k)*bernoulli(k)*z**(n-k) for k in xrange(0,n+1))
    n = int(n)
    if n < 0:
        raise ValueError("Bernoulli polynomials only defined for n >= 0")
    if z == 0 or (z == 1 and n > 1):
        return bernoulli(n)[-1]
    if z == 0.5:
        return (np.ldexp(1,1-n)-1)*bernoulli(n)[-1]
    if n <= 3:
        if n == 0: return z ** 0
        if n == 1: return z - 0.5
        if n == 2: return (6*z*(z-1)+1)/6
        if n == 3: return z*(z*(z-1.5)+0.5)
    if np.isinf(z):
        return z ** n
    if np.isnan(z):
        return z
    if abs(z) > 2:
        def terms():
            t = 1
            yield t
            r = 1/z
            k = 1
            while k <= n:
                t = t*(n+1-k)/k*r
                if not (k > 2 and k & 1):
                    yield t*bernoulli(k)[-1]
                k += 1
        return np.sum(terms) * z**n
    else:
        def terms():
            yield bernoulli(n)[-1]
            t = 1
            k = 1
            while k <= n:
                t = t*(n+1-k)/k * z
                m = n-k
                if not (m > 2 and m & 1):
                    yield t*bernoulli(m)[-1]
                k += 1
        return np.sum(terms)

def polylog_series(s, z,tol=EPS):
    l = 0.
    k = 1
    zk = z
    while 1:
        term = zk / k**s
        l += term
        if abs(term) < tol:
            break
        zk *= z
        k += 1
    return l

def polylog_continuation(n, z):
    if n < 0:
        return z*0
    twopij = 2j * pi
    a = -twopij**n/fac(n) * bernpoly(n, ln(z)/twopij)
    if np.isreal(z) and z.real < 0:
        a = a.real
    if z.imag < 0 or (z.imag == 0 and z.real >= 1):
        a -= twopij*ln(z)**(n-1)/fac(n-1)
    return a

@cython.cdivision(True)
cdef double complex polylog_unitcircle(int n, double complex z, double tol=EPS):
    if n > 1:
        l = 0.
        logz = ln(z)
        logmz = 1
        m = 0
        while 1:
            if (n-m) != 1:
                term = zetan(n-m) * logmz / fac(m)
                if np.isnan(term):break
                if term and abs(term) < EPS:
                    break
                l += term
            logmz *= logz
            m += 1
        l += ln(z)**(n-1)/fac(n-1)*(harmonic(n-1)-ln(-ln(z)))
    elif n < 1:  # else
        logz = ln(z)
        l = fac(-n)*(-logz)**(n-1)
        logkz = 1
        k = 0
        while 1:
            b = bernoulli(k-n+1)[-1]
            if b:
                term = b*logkz/(fac(k)*(k-n+1))
                if abs(term) < tol:
                    break
                l -= term
            logkz *= logz
            k += 1
    else:
        raise ValueError
    if np.isreal(z) and z.real < 0:
        l = l.real
    return l

def polylog_general(s, z,tol=EPS):
    v = 0
    u = ln(z)
    if not abs(u) < 5: # theoretically |u| < 2*pi
        raise NotImplementedError("polylog for arbitrary s and z")
    t = 1
    k = 0
    while 1:
        term = zetan(s-k) * t
        if abs(term) < tol:
            break
        v += term
        k += 1
        t *= u
        t /= k
    return gamma(1-s)*(-u)**(s-1) + v


@cython.cdivision(True)
def _polylog(s, z, tol=EPS):
    if z == 1:
        if isinstance(s,complex):
            raise NotImplementedError
        return zetan(s)
    if z == -1:
        return -altzeta(s)
    if s == 0:
        return z/(1-z)
    if s == 1:
        value = -ln(1-z)
        if np.isnan(value):
            value = -np.log(1-np.complex(z))
        return value
    if s == -1:
        return z/(1-z)**2
    if abs(z) <= 0.75 or (not isinstance(s,int) and abs(z) < 0.9):
        return polylog_series(s, z,tol)
    if abs(z) >= 1.4 and isinstance(s,int):
        return (-1)**(s+1)*polylog_series(s, 1/z,tol) + polylog_continuation(s, z)
    if isinstance(s,int):
        return polylog_unitcircle(int(s), z, tol)
    return polylog_general(s, z)

    #raise NotImplementedError("polylog for arbitrary s and z")
    # This could perhaps be used in some cases
    #from quadrature import quad
    #return quad(lambda t: t**(s-1)/(exp(t)/z-1),[0,inf])/gamma(s)


def polylog(s, z, tol=EPS):
    return _polylog(s, z, tol)