from __future__ import division
import cython
import warnings
from special import bernoulli,gamma,digamma
import numpy as np
from math import pi

cdef extern from "math.h":
    double exp(double x) nogil
    double log(double x) nogil
    double fabs(double complex x) nogil
    double M_LN2

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

cdef int fac(int n) : return gamma(n+1)

cdef double complex zetan(double s): return 1+zetac(s)

#TO DO : DEAL WITH POW
cdef double complex altzeta(double s):
    if s == 1:
        return M_LN2 + 0*s
    return (1-pow(2, 1-s)) * zetan(s)

cdef double complex harmonic(int x):
    if x == 0 or x == 1:
        return x
    return digamma(x+1) - digamma(1)

cpdef double complex ln(double complex x):
    cdef npy_cdouble r
    r = npy_clog((<npy_cdouble*>&x)[0])
    return (<double complex*>&r)[0]

@cython.cdivision(True)
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
    if fabs(z) > 2:
        def terms():
            cdef double complex t = 1
            yield t
            cdef double complex r = 1/z
            cdef unsigned int k = 1
            while k <= n:
                t = t*(n+1-k)/k*r
                if not (k > 2 and k & 1):
                    yield t*bernoulli(k)[-1]
                k += 1
        return np.sum(terms) * z**n
    else:
        def terms():
            yield bernoulli(n)[-1]
            cdef double complex t = 1
            cdef unsigned int k = 1
            cdef int m=0
            while k <= n:
                t = t*(n+1-k)/k * z
                m = n-k
                if not (m > 2 and m & 1):
                    yield t*bernoulli(m)[-1]
                k += 1
        return np.sum(terms)


@cython.cdivision(True)
def polylog_series(s, z, tol=EPS):
    print s,z
    l = 0.
    k = 1
    zk = z
    term = 0.
    while 1:
        print term
        term = zk / k**s
        l += term
        if abs(term) < tol:
            break
        zk *= z
        k += 1
    return l

@cython.cdivision(True)
cdef double complex polylog_continuation(int n, double complex z):
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
                if term and fabs(term) < EPS:
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
                if fabs(term) < tol:
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
    if not fabs(u) < 5: # theoretically |u| < 2*pi
        raise NotImplementedError("polylog for arbitrary s and z")
    t = 1
    k = 0
    while 1:
        term = zetan(s-k) * t
        if fabs(term) < tol:
            break
        v += term
        k += 1
        t *= u
        t /= k
    return gamma(1-s)*(-u)**(s-1) + v


@cython.cdivision(True)
def _polylog_scalar(s, z, tol=EPS):
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
        if value.imag==0: return value.real
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
    r"""
Computes the polylogarithm, defined by the sum

.. math ::

    \mathrm{Li}_s(z) = \sum_{k=1}^{\infty} \frac{z^k}{k^s}.

This series is convergent only for `|z| < 1`, so elsewhere
the analytic continuation is implied.

The polylogarithm should not be confused with the logarithmic
integral (also denoted by Li or li), which is implemented
as :func:`~mpmath.li` in the mpmath package.

**Examples**

The polylogarithm satisfies a huge number of functional identities.
A sample of polylogarithm evaluations is shown below::

    >>> polylog(1,0.5), log(2)
    (0.693147180559945, 0.693147180559945)
    >>> polylog(2,0.5), (pi**2-6*log(2)**2)/12
    (0.582240526465012, 0.582240526465012)
    >>> polylog(2,-phi), -log(phi)**2-pi**2/10
    (-1.21852526068613, -1.21852526068613)
    >>> polylog(3,0.5), 7*zeta(3)/8-pi**2*log(2)/12+log(2)**3/6
    (0.53721319360804, 0.53721319360804)

:func:`polylog` can evaluate the analytic continuation of the
polylogarithm when `s` is an integer::

    >>> polylog(2, 10)
    (0.536301287357863 - 7.23378441241546j)
    >>> polylog(2, -10)
    -4.1982778868581
    >>> polylog(2, 10j)
    (-3.05968879432873 + 3.71678149306807j)
    >>> polylog(-2, 10)
    -0.150891632373114
    >>> polylog(-2, -10)
    0.067618332081142
    >>> polylog(-2, 10j)
    (0.0384353698579347 + 0.0912451798066779j)

Some more examples, with arguments on the unit circle (note that
the series definition cannot be used for computation here)::

    >>> polylog(2,j)
    (-0.205616758356028 + 0.915965594177219j)
    >>> j*catalan-pi**2/48
    (-0.205616758356028 + 0.915965594177219j)
    >>> polylog(3,exp(2*pi*j/3))
    (-0.534247512515375 + 0.765587078525922j)
    >>> -4*zeta(3)/9 + 2*j*pi**3/81
    (-0.534247512515375 + 0.765587078525921j)

Polylogarithms of different order are related by integration
and differentiation::

    >>> s, z = 3, 0.5
    >>> polylog(s+1, z)
    0.517479061673899
    >>> quad(lambda t: polylog(s,t)/t, [0, z])
    0.517479061673899
    >>> z*diff(lambda t: polylog(s+2,t), z)
    0.517479061673899

Taylor series expansions around `z = 0` are::

    >>> for n in range(-3, 4):
    ...     nprint(taylor(lambda x: polylog(n,x), 0, 5))
    ...
    [0.0, 1.0, 8.0, 27.0, 64.0, 125.0]
    [0.0, 1.0, 4.0, 9.0, 16.0, 25.0]
    [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    [0.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    [0.0, 1.0, 0.5, 0.333333, 0.25, 0.2]
    [0.0, 1.0, 0.25, 0.111111, 0.0625, 0.04]
    [0.0, 1.0, 0.125, 0.037037, 0.015625, 0.008]

The series defining the polylogarithm is simultaneously
a Taylor series and an L-series. For certain values of `z`, the
polylogarithm reduces to a pure zeta function::

    >>> polylog(pi, 1), zeta(pi)
    (1.17624173838258, 1.17624173838258)
    >>> polylog(pi, -1), -altzeta(pi)
    (-0.909670702980385, -0.909670702980385)

Evaluation for arbitrary, nonintegral `s` is supported
for `z` within the unit circle:

    >>> polylog(3+4j, 0.25)
    (0.24258605789446 - 0.00222938275488344j)
    >>> nsum(lambda k: 0.25**k / k**(3+4j), [1,inf])
    (0.24258605789446 - 0.00222938275488344j)

It is also currently supported outside of the unit circle for `z`
not too large in magnitude::

    >>> polylog(1+j, 20+40j)
    (-7.1421172179728 - 3.92726697721369j)
    >>> polylog(1+j, 200+400j)
    Traceback (most recent call last):
      ...
    NotImplementedError: polylog for arbitrary s and z

**References**

1. Richard Crandall, "Note on fast polylogarithm computation"
   http://people.reed.edu/~crandall/papers/Polylog.pdf
2. http://en.wikipedia.org/wiki/Polylogarithm
3. http://mathworld.wolfram.com/Polylogarithm.html

   """
    return _polylog_scalar(s, z, tol)