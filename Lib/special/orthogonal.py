#!/usr/bin/env python
#
# Author:  Travis Oliphant 2000
"""
A collection of functions to find the weights and abscissas for
Gaussian Quadrature.

These calculations are done by finding the eigenvalues of a
tridiagonal matrix whose entries are dependent on the coefficients
in the recursion formula for the orthogonal polynomials with the
corresponding weighting function over the interval.

Many recursion relations for orthogonal polynomials are given:

a1n f_n+1 (x) = (a2n + a3n x ) f_n (x) - a4n f_n-1 (x)

The recursion relation of interest is

P_n+1 (x) = (x - A_n) P_n (x) - B_n P_n-1 (x)

where P has a different normalization than f.

The coefficients can be found as:

A_n = -a2n / a3n

B_n = ( a4n / a3n sqrt(h_n-1 / h_n))**2

     where
             h_n = int_a^b w(x) f_n(x)^2
assume:
P_0(x) = 1
P_-1(x) == 0
             
See Numerical Recipies in C, page 156 and
Abramowitz and Stegun p. 774, 782

Functions:

  gen_roots_and_weights  -- Generic roots and weights.
  J_roots                -- Jacobi 
  Js_roots               -- Shifted Jacobi
  La_roots               -- Generalized Laguerre
  H_roots                -- Hermite
  He_roots               -- Hermite (unit-variance)
  Cg_roots               -- Ultraspherical (Gegenbauer)
  T_roots                -- Chebyshev of the first kind
  U_roots                -- Chebyshev of the second kind
  C_roots                -- Chebyshev of the first kind ([-2,2] interval)
  S_roots                -- Chebyshev of the second kind ([-2,2] interval)
  Ts_roots               -- Shifted Chebyshev of the first kind.
  Us_roots               -- Shifted Chebyshev of the second kind.
  P_roots                -- Legendre
  Ps_roots               -- Shifted Legendre
  L_roots                -- Laguerre
"""

from __future__ import nested_scopes
from scipy.linalg import eig
from scipy_base import *
import cephes
_gam = cephes.gamma


class orthopoly1d(poly1d):
    def __init__(self, roots, weights=None, hn=1.0, An=1.0, wfunc=None, limits=None, monic=0):
        poly1d.__init__(self, roots, r=1)
        self.__dict__['weights'] = weights
        self.__dict__['weight_func'] = wfunc
        self.__dict__['limits'] = limits
        mu = sqrt(hn)
        if monic:
            mu = mu / abs(An)
            An = 1.0
        self.__dict__['normcoef'] = mu
        self.__dict__['coeffs'] *= An
    
def gen_roots_and_weights(n,an_func,sqrt_bn_func,mu):
    """[x,w] = gen_roots_and_weights(n,an_func,sqrt_bn_func,mu)

    Returns the roots (x) of an nth order orthogonal polynomail,
    and weights (w) to use in appropriate Gaussian quadrature with that
    orthogonal polynomial.

    The polynomials have the recurrence relation
          P_n+1(x) = (x - A_n) P_n(x) - B_n P_n-1(x)

    an_func(n)          should return A_n
    sqrt_bn_func(n)     should return sqrt(B_n)
    mu ( = h_0 )        is the integral of the weight over the orthogonal interval
    """
    nn = arange(1.0,n)
    sqrt_bn = sqrt_bn_func(nn)
    an = an_func(concatenate(([0],nn)))
    [x,v] = eig((diag(an)+diag(sqrt_bn,1)+diag(sqrt_bn,-1)))
    answer = []
    sortind = argsort(x)
    answer.append(take(x,sortind))
    answer.append(take(mu*v[:,0]**2,sortind))
    return answer    

# Jacobi Polynomials 1               P^(alpha,beta)_n(x)
def J_roots(n,alpha,beta,mu=0):
    """[x,w] = J_roots(n,alpha,beta)

    Returns the roots (x) of the nth order Jacobi polynomial, P^(alpha,beta)_n(x)
    and weights (w) to use in Gaussian Quadrature over [-1,1] with weighting
    function (1-x)**alpha (1+x)**beta with alpha,beta > -1.
    """
    if (alpha <= -1) or (beta <= -1):
        raise ValueError, "alpha and beta must be greater than -1."
    assert(n>0), "n must be positive."

    (p,q) = (alpha,beta)
    # from recurrence relation
    sbn_J = lambda k: 2.0/(2*k+p+q)*sqrt((k+p)*(k+q)*k*(k+p+q)/((2*k+p+q+1)*(2*k+p+q-1)))
    if (p == q):
        an_J = lambda k: 0*k
    else:
        an_J = lambda k: (q*q - p*p)/((2.0*k+p+q)*(2.0*k+p+q+2))
    g = cephes.gamma
    mu0 = 2**(p+q+1)*g(p+1)*g(q+1)/(g(p+q+2))
    val = gen_roots_and_weights(n,an_J,sbn_J,mu0)
    if mu:
        return val + [mu0]
    else:
        return val

def jacobi(n,alpha,beta,monic=0):
    """Returns the nth order Jacobi polynomial, P^(alpha,beta)_n(x)
    orthogonal over [-1,1] with weighting function
    (1-x)**alpha (1+x)**beta with alpha,beta > -1.
    """
    assert(n>=0), "n must be nonnegative"
    if n==0: n1 = n+1
    else: n1 = n
    x,w,mu = J_roots(n1,alpha,beta,mu=1)
    if n==0: x,w = [],[]
    wfunc = lambda x: (1-x)**alpha * (1+x)**beta
    ab1 = alpha+beta+1.0
    hn = 2**ab1/(2*n+ab1)*_gam(n+alpha+1)
    hn *= _gam(n+beta+1.0) / _gam(n+1) / _gam(n+ab1)
    An = _gam(2*n+ab1)/2.0**n / _gam(n+1) / _gam(n+ab1)
    p = orthopoly1d(x,w,hn,An,wfunc,(-1,1),monic)
    return p

# Jacobi Polynomials shifted         G_n(p,q,x)
def Js_roots(n,p1,q1,mu=0):
    """[x,w] = Js_roots(n,p,q)

    Returns the roots (x) of the nth order shifted Jacobi polynomial, G_n(p,q,x),
    and weights (w) to use in Gaussian Quadrature over [0,1] with weighting
    function (1-x)**(p-q) x**(q-1) with p-q > -1 and q > 0.
    """
    # from recurrence relation
    if not ( ( (p1 - q1) > -1 ) and ( q1 > 0 ) ):
        raise ValueError, "(p - q) > -1 and q > 0 please."
    assert(n>0), "n must be positive."
    
    (p,q) = (p1,q1)
    if (p == 0):
        sbn_Js = lambda k: sqrt((k+q-1)*(k-q))/((2*k-1)*2)
    else:
        sbn_Js = lambda k: sqrt((k*(k+q-1)*(k+p-1)*(k+p-q))/((2*k+p)*(2*k+p-2)))/(2*k+p-1)
    if (p == 1):
        an_Js = lambda k: 0*k + 0.5
    else:
        an_Js = lambda k: (2*k*(k+p)+q*(p-1))/((2*k+p+1)*(2*k+p-1))

    # integral of weight over interval
    g = cephes.gamma
    mu0 =  g(q)*g(p-q+1)/g(p+1)
    val = gen_roots_and_weights(n,an_Js,sbn_Js,mu0)
    if mu:
        return val + [mu0]
    else:
        return val

def sh_jacobi(n, p, q, monic=0):
    """Returns the nth order Jacobi polynomial, G_n(p,q,x)
    orthogonal over [0,1] with weighting function
    (1-x)**(p-q) (1+x)**(q-1) with p>q-1 and q > 0.
    """
    assert(n>=0), "n must be nonnegative"
    if n==0: n1 = n+1
    else: n1 = n
    x,w,mu0 = Js_roots(n1,p,q,mu=1)
    if n==0: x,w = [],[]
    wfunc = lambda x: (1.0-x)**(p-q) * (1+x)**(q-1.)
    hn = _gam(n+1)*_gam(n+q)*_gam(n+p)*_gam(n+p-q+1)
    hn /= (2*n+p)*(_gam(2*n+p)**2)
    An = 1.0
    p = orthopoly1d(x,w,hn,An,wfunc=wfunc,limits=(0,1),monic=monic)
    return p

# Generalized Laguerre               L^(alpha)_n(x)
def La_roots(n,alpha,mu=0):
    """[x,w] = La_roots(n,alpha)

    Returns the roots (x) of the nth order generalized (associated) Laguerre
    polynomial, L^(alpha)_n(x), and weights (w) to use in Gaussian quadrature over
    [0,inf] with weighting function exp(-x) x**alpha with alpha > -1.
    """
    if not (alpha > -1):
        raise ValueError, "alpha > -1"
    assert(n>0), "n must be positive."
    (p,q) = (alpha,0.0)
    sbn_La = lambda k: -sqrt(k*(k + p))  # from recurrence relation
    an_La = lambda k: 2*k + p + 1                 
    mu0 = cephes.gamma(alpha+1)           # integral of weight over interval 
    val = gen_roots_and_weights(n,an_La,sbn_La,mu0)
    if mu:
        return val + [mu0]
    else:
        return val

def genlaguerre(n,alpha,monic=0):
    """Returns the nth order generalized (associated) Laguerre polynomial,
    L^(alpha)_n(x), orthogonal over [0,inf) with weighting function
    exp(-x) x**alpha with alpha > -1
    """
    assert(n>=0), "n must be nonnegative"
    if n==0: n1 = n+1
    else: n1 = n
    x,w,mu0 = La_roots(n1,alpha,mu=1)
    wfunc = lambda x: exp(-x) * x**alpha
    if n==0: x,w = [],[]
    hn = _gam(n+alpha+1)/_gam(n+1)
    An = (-1)**n / _gam(n+1)
    p = orthopoly1d(x,w,hn,An,wfunc,(0,inf),monic)
    return p

# Hermite  1                         H_n(x)
def H_roots(n,mu=0):
    """[x,w] = H_roots(n)

    Returns the roots (x) of the nth order Hermite polynomial,
    H_n(x), and weights (w) to use in Gaussian Quadrature over
    [-inf,inf] with weighting function exp(-x**2).
    """
    assert(n>0), "n must be positive."
    sbn_H = lambda k: sqrt(k/2)  # from recurrence relation
    an_H = lambda k: 0*k                    
    mu0 = sqrt(pi)               # integral of weight over interval 
    val = gen_roots_and_weights(n,an_H,sbn_H,mu0)
    if mu:
        return val + [mu0]
    else:
        return val

def hermite(n,monic=0):
    """Return the nth order Hermite polynomial, H_n(x), orthogonal over
    (-inf,inf) with weighting function exp(-x**2)
    """
    assert(n>=0), "n must be nonnegative"
    if n==0: n1 = n+1
    else: n1 = n
    x,w,mu0 = H_roots(n1,mu=1)
    wfunc = lambda x: exp(-x*x)
    if n==0: x,w = [],[]
    hn = 2**n * _gam(n+1)*sqrt(pi)
    An = 2**n
    p = orthopoly1d(x,w,hn,An,wfunc,(-inf,inf),monic)
    return p
    
# Hermite  2                         He_n(x)
def He_roots(n,mu=0):
    """[x,w] = He_roots(n)

    Returns the roots (x) of the nth order Hermite polynomial,
    He_n(x), and weights (w) to use in Gaussian Quadrature over
    [-inf,inf] with weighting function exp(-(x/2)**2).
    """
    assert(n>0), "n must be positive."
    sbn_He = lambda k: sqrt(k)   # from recurrence relation
    an_He  = lambda k: 0*k                
    mu0 = sqrt(2*pi)             # integral of weight over interval 
    val = gen_roots_and_weights(n,an_He,sbn_He,mu0)
    if mu:
        return val + [mu0]
    else:
        return val

def hermitenorm(n,monic=0):
    """Return the nth order normalized Hermite polynomial, He_n(x), orthogonal
    over (-inf,inf) with weighting function exp(-(x/2)**2)
    """
    assert(n>=0), "n must be nonnegative"
    if n==0: n1 = n+1
    else: n1 = n
    x,w,mu0 = He_roots(n1,mu=1)
    wfunc = lambda x: exp(-x*x/4.0)
    if n==0: x,w = [],[]
    hn = sqrt(2*pi)*_gam(n+1)
    An = 1.0
    p = orthopoly1d(x,w,hn,An,wfunc=wfunc,limits=(-inf,inf),monic=monic)
    return p

## The remainder of the polynomials can be derived from the ones above.

# Ultraspherical (Gegenbauer)        C^(alpha)_n(x)
def Cg_roots(n,alpha,mu=0):
    """[x,w] = Cg_roots(n,alpha)

    Returns the roots (x) of the nth order Ultraspherical (Gegenbauer)
    polynomial, C^(alpha)_n(x), and weights (w) to use in Gaussian Quadrature
    over [-1,1] with weighting function (1-x**2)**(alpha-1/2) with alpha>-1/2.
    """
    return J_roots(n,alpha-0.5,alpha-0.5,mu=mu)

def gegenbauer(n,alpha,monic=0):
    """Return the nth order Gegenbauer (ultraspherical) polynomial,
    C^(alpha)_n(x), orthogonal over [-1,1] with weighting function
    (1-x**2)**(alpha-1/2) with alpha > -1/2
    """
    return jacobi(n,alpha-0.5,alpha-0.5,monic=monic)

# Chebyshev of the first kind        T_n(x)
def T_roots(n,mu=0):
    """[x,w] = T_roots(n)

    Returns the roots (x) of the nth order Chebyshev (of the first kind)
    polynomial, T_n(x), and weights (w) to use in Gaussian Quadrature
    over [-1,1] with weighting function (1-x**2)**(-1/2).
    """
    return J_roots(n,-0.5,-0.5,mu=mu)

def chebyt(n,monic=0):
    """Return nth order Chebyshev polynomial of first kind, Tn(x).  Orthogonal
    over [-1,1] with weight function (1-x**2)**(-1/2).
    """
    return jacobi(n,-0.5,-0.5,monic=monic)

# Chebyshev of the second kind       U_n(x)
def U_roots(n,mu=0):
    """[x,w] = U_roots(n)

    Returns the roots (x) of the nth order Chebyshev (of the second kind)
    polynomial, U_n(x), and weights (w) to use in Gaussian Quadrature
    over [-1,1] with weighting function (1-x**2)**1/2.
    """
    return J_roots(n,0.5,0.5,mu=mu)

def chebyu(n):
    """Return nth order Chebyshev polynomial of second kind, Un(x).  Orthogonal
    over [-1,1] with weight function (1-x**2)**(1/2).
    """
    return jacobi(n,0.5,0.5,monic=monic)

# Chebyshev of the first kind        C_n(x)
def C_roots(n,mu=0):
    """[x,w] = C_roots(n)

    Returns the roots (x) of the nth order Chebyshev (of the first kind)
    polynomial, C_n(x), and weights (w) to use in Gaussian Quadrature
    over [-2,2] with weighting function (1-(x/2)**2)**(-1/2).
    """
    if mu:
        [x,w,mu0] = J_roots(n,-0.5,-0.5,mu=1)
        return [x*2,w,mu0]
    else:
        [x,w] = J_roots(n,-0.5,-0.5,mu=0)
        return [x*2,w]

def chebyc(n,monic=0):
    """Return nth order Chebyshev polynomial of first kind, Cn(x).  Orthogonal
    over [-2,2] with weight function (1-(x/2)**2)**(-1/2).
    """
    assert(n>=0), "n must be nonnegative"
    if n==0: n1 = n+1
    else: n1 = n
    x,w,mu0 = C_roots(n1,mu=1)
    if n==0: x,w = [],[]
    hn = 4*pi * ((n==0)+1)
    An = 1.0
    p = orthopoly1d(x,w,hn,An,wfunc=lambda x: 1.0/sqrt(1-x*x/4.0),limits=(-2,2),monic=monic)
    if not monic:
        p = p * 2.0/p(2)
    return p

# Chebyshev of the second kind       S_n(x)
def S_roots(n,mu=0):
    """[x,w] = S_roots(n)

    Returns the roots (x) of the nth order Chebyshev (of the second kind)
    polynomial, S_n(x), and weights (w) to use in Gaussian Quadrature
    over [-2,2] with weighting function (1-(x/2)**2)**1/2.
    """
    if mu:
        [x,w,mu0] = J_roots(n,0.5,0.5,mu=1)
        return [x*2,w,mu0]
    else:
        [x,w] = J_roots(n,0.5,0.5,mu=0)
        return [x*2,w]

def chebys(n,monic=0):
    """Return nth order Chebyshev polynomial of second kind, Sn(x).  Orthogonal
    over [-2,2] with weight function (1-(x/)**2)**(1/2).
    """
    assert(n>=0), "n must be nonnegative"
    if n==0: n1 = n+1
    else: n1 = n
    x,w,mu0 = S_roots(n1,mu=1)
    if n==0: x,w = [],[]
    hn = pi
    An = 1.0
    p = orthopoly1d(x,w,hn,An,wfunc=lambda x: sqrt(1-x*x/4.0),limits=(-2,2),monic=monic)
    if not monic:
        p = p * (n+1.0)/p(2)
    return p

# Shifted Chebyshev of the first kind     T^*_n(x)
def Ts_roots(n,mu=0):
    """[x,w] = Ts_roots(n)

    Returns the roots (x) of the nth order shifted Chebyshev (of the first kind)
    polynomial, T^*_n(x), and weights (w) to use in Gaussian Quadrature
    over [0,1] with weighting function (x-x**2)**(-1/2).
    """
    return Js_roots(n,0.0,0.5,mu=mu)

def sh_chebyt(n,monic=0):
    """Return nth order shifted Chebyshev polynomial of first kind, Tn(x).
    Orthogonal over [0,1] with weight function (x-x**2)**(-1/2).
    """
    return sh_jacobi(n,0.0,0.5,monic=monic)


# Shifted Chebyshev of the second kind    U^*_n(x)
def Us_roots(n,mu=0):
    """[x,w] = Us_roots(n)

    Returns the roots (x) of the nth order shifted Chebyshev (of the second kind)
    polynomial, U^*_n(x), and weights (w) to use in Gaussian Quadrature
    over [0,1] with weighting function (x-x**2)**1/2.
    """
    return Js_roots(n,2.0,1.5,mu=mu)

def sh_chebyu(n,monic=0):
    """Return nth order shifted Chebyshev polynomial of second kind, Un(x).
    Orthogonal over [0,1] with weight function (x-x**2)**(1/2).
    """
    return sh_jacobi(n,2.0,1.5,monic=monic)

# Legendre 
def P_roots(n,mu=0):
    """[x,w] = P_roots(n)

    Returns the roots (x) of the nth order Legendre polynomial, P_n(x),
    and weights (w) to use in Gaussian Quadrature over [-1,1] with weighting
    function 1.
    """
    return J_roots(n,0.0,0.0,mu=mu)

def legendre(n,monic=0):
    """Returns the nth order Legendre polynomial, P_n(x), orthogonal over
    [-1,1] with weight function 1.
    """
    assert(n>=0), "n must be nonnegative"
    if n==0: n1 = n+1
    else: n1 = n
    x,w,mu0 = P_roots(n1,mu=1)
    if n==0: x,w = [],[]
    hn = 2.0/(2*n+1)
    An = _gam(2*n+1)/_gam(n+1)**2 / 2.0**n
    p = orthopoly1d(x,w,hn,An,wfunc=lambda x: 1.0,limits=(-1,1),monic=monic)
    return p

# Shifted Legendre              P^*_n(x)
def Ps_roots(n,mu=0):
    """[x,w] = Ps_roots(n)

    Returns the roots (x) of the nth order shifted Legendre polynomial, P^*_n(x),
    and weights (w) to use in Gaussian Quadrature over [0,1] with weighting
    function 1.
    """
    return Js_roots(n,1.0,1.0,mu=mu)

def sh_legendre(n,monic=0):
    """Returns the nth order shifted Legendre polynomial, P^*_n(x), orthogonal
    over [0,1] with weighting function 1.
    """
    assert(n>=0), "n must be nonnegative"
    if n==0: n1 = n+1
    else: n1 = n
    x,w,mu0 = Ps_roots(n1,mu=1)
    if n==0: x,w = [],[]
    hn = 1.0/(2*n+1.0)
    An = 1.0
    p = orthopoly1d(x,w,hn,An,wfunc=lambda x: 1.0,limits=(0,1),monic=monic)
    return p 

# Laguerre                      L_n(x)
def L_roots(n,mu=0):
    """[x,w] = L_roots(n)

    Returns the roots (x) of the nth order Laguerre polynomial, L_n(x),
    and weights (w) to use in Gaussian Quadrature over [0,inf] with weighting
    function exp(-x).
    """
    return La_roots(n,0.0,mu=mu)

def laguerre(n,monic=0):
    """Return the nth order Laguerre polynoimal, L_n(x), orthogonal over
    [0,inf) with weighting function exp(-x)
    """
    assert(n>=0), "n must be nonnegative"
    if n==0: n1 = n+1
    else: n1 = n
    x,w,mu0 = L_roots(n1,mu=1)
    if n==0: x,w = [],[]
    hn = 1.0
    An = (-1)**n / _gam(n+1)
    p = orthopoly1d(x,w,hn,An,lambda x: exp(-x),(0,inf),monic)
    return p 
    







