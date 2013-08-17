# Docstrings for generated ufuncs
#
# The syntax is designed to look like the function add_newdoc is being
# called from numpy.lib, but in this file add_newdoc puts the
# docstrings in a dictionary. This dictionary is used in
# generate_ufuncs.py to generate the docstrings for the ufuncs in
# scipy.special at the C level when the ufuncs are created at compile
# time.

from __future__ import division, print_function, absolute_import

docdict = {}


def get(name):
    return docdict.get(name)


def add_newdoc(place, name, doc):
    docdict['.'.join((place, name))] = doc


add_newdoc("scipy.special", "_lambertw",
    """
    Internal function, use `lambertw` instead.
    """)

add_newdoc("scipy.special", "airy",
    """
    (Ai,Aip,Bi,Bip)=airy(z) calculates the Airy functions and their derivatives
    evaluated at real or complex number z.  The Airy functions Ai and Bi
    are two independent solutions of y''(x)=xy.  Aip and Bip are the first derivatives
    evaluated at x of Ai and Bi respectively.
    """)

add_newdoc("scipy.special", "airye",
    """
    (Aie,Aipe,Bie,Bipe)=airye(z) calculates the exponentially scaled Airy functions and
    their derivatives evaluated at real or complex number z.
    airye(z)[0:1] = airy(z)[0:1] * exp(2.0/3.0*z*sqrt(z))
    airye(z)[2:3] = airy(z)[2:3] * exp(-abs((2.0/3.0*z*sqrt(z)).real))
    """)

add_newdoc("scipy.special", "bdtr",
    """
    y=bdtr(k,n,p) returns the sum of the terms 0 through k of the
    Binomial probability density:  sum(nCj p**j (1-p)**(n-j),j=0..k)
    """)

add_newdoc("scipy.special", "bdtrc",
    """
    y=bdtrc(k,n,p) returns the sum of the terms k+1 through n of the
    Binomial probability density: sum(nCj p**j (1-p)**(n-j), j=k+1..n)
    """)

add_newdoc("scipy.special", "bdtri",
    """
    p=bdtri(k,n,y) finds the probability p such that the sum of the
    terms 0 through k of the Binomial probability density is equal to the
    given cumulative probability y.
    """)

add_newdoc("scipy.special", "bdtrik",
    """
    """)

add_newdoc("scipy.special", "bdtrin",
    """
    """)

add_newdoc("scipy.special", "binom",
    """
    binom(n, k)

    Binomial coefficient
    """)

add_newdoc("scipy.special", "btdtria",
    """
    """)

add_newdoc("scipy.special", "btdtrib",
    """
    """)

add_newdoc("scipy.special", "bei",
    """
    y=bei(x) returns the Kelvin function bei x
    """)

add_newdoc("scipy.special", "beip",
    """
    y=beip(x) returns the derivative of the Kelvin function bei x
    """)

add_newdoc("scipy.special", "ber",
    """
    y=ber(x) returns the Kelvin function ber x
    """)

add_newdoc("scipy.special", "berp",
    """
    y=berp(x) returns the derivative of the Kelvin function ber x
    """)

add_newdoc("scipy.special", "besselpoly",
    """
    y=besselpoly(a,lam,nu) returns the value of the integral:
    integral(x**lam * jv(nu,2*a*x),x=0..1).
    """)

add_newdoc("scipy.special", "beta",
    """
    y=beta(a,b) returns gamma(a) * gamma(b) / gamma(a+b)
    """)

add_newdoc("scipy.special", "betainc",
    """
    betainc(a, b, x)

    Compute the incomplete beta integral of the arguments, evaluated
    from zero to x::

        gamma(a+b) / (gamma(a)*gamma(b)) * integral(t**(a-1) (1-t)**(b-1), t=0..x).

    Notes
    -----
    The incomplete beta is also sometimes defined without the terms
    in gamma, in which case the above definition is the so-called regularized
    incomplete beta. Under this definition, you can get the incomplete beta by
    multiplying the result of the scipy function by beta(a, b).

    """)

add_newdoc("scipy.special", "betaincinv",
    """
    betaincinv(a,b,y)

    Compute x such that betainc(a,b,x) = y.
    """)

add_newdoc("scipy.special", "betaln",
    """
    y=betaln(a,b) returns the natural logarithm of the absolute value of
    beta: ln(abs(beta(x))).
    """)

add_newdoc("scipy.special", "btdtr",
    """
    y=btdtr(a,b,x) returns the area from zero to x under the beta
    density function: gamma(a+b)/(gamma(a)*gamma(b)))*integral(t**(a-1)
    (1-t)**(b-1), t=0..x).  SEE ALSO betainc
    """)

add_newdoc("scipy.special", "btdtri",
    """
    x=btdtri(a,b,p) returns the pth quantile of the beta distribution.  It is
    effectively the inverse of btdtr returning the value of x for which
    btdtr(a,b,x) = p.   SEE ALSO betaincinv
    """)

add_newdoc("scipy.special", "cbrt",
    """
    y=cbrt(x) returns the real cube root of x.
    """)

add_newdoc("scipy.special", "chdtr",
    """
    p=chdtr(v,x) Returns the area under the left hand tail (from 0 to x) of the Chi
    square probability density function with v degrees of freedom:
    1/(2**(v/2) * gamma(v/2)) * integral(t**(v/2-1) * exp(-t/2), t=0..x)
    """)

add_newdoc("scipy.special", "chdtrc",
    """
    p=chdtrc(v,x) returns the area under the right hand tail (from x to
    infinity) of the Chi square probability density function with v
    degrees of freedom:
    1/(2**(v/2) * gamma(v/2)) * integral(t**(v/2-1) * exp(-t/2), t=x..inf)
    """)

add_newdoc("scipy.special", "chdtri",
    """
    x=chdtri(v,p) returns the argument x such that chdtrc(v,x) is equal
    to p.
    """)

add_newdoc("scipy.special", "chdtriv",
    """
    """)

add_newdoc("scipy.special", "chndtr",
    """
    """)

add_newdoc("scipy.special", "chndtrix",
    """
    """)

add_newdoc("scipy.special", "chndtridf",
    """
    """)

add_newdoc("scipy.special", "chndtrinc",
    """
    """)

add_newdoc("scipy.special", "cosdg",
    """
    y=cosdg(x) calculates the cosine of the angle x given in degrees.
    """)

add_newdoc("scipy.special", "cosm1",
    """
    y=calculates cos(x) - 1 for use when x is near zero.
    """)

add_newdoc("scipy.special", "cotdg",
    """
    y=cotdg(x) calculates the cotangent of the angle x given in degrees.
    """)

add_newdoc("scipy.special", "dawsn",
    """
    y=dawsn(x) returns dawson's integral: exp(-x**2) *
    integral(exp(t**2),t=0..x).

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva
    """)

add_newdoc("scipy.special", "ellipe",
    """
    y=ellipe(m) returns the complete integral of the second kind:
    integral(sqrt(1-m*sin(t)**2),t=0..pi/2)
    """)

add_newdoc("scipy.special", "ellipeinc",
    """
    y=ellipeinc(phi,m) returns the incomplete elliptic integral of the
    second kind: integral(sqrt(1-m*sin(t)**2),t=0..phi)
    """)

add_newdoc("scipy.special", "ellipj",
    """
    (sn,cn,dn,ph)=ellipj(u,m) calculates the Jacobian elliptic functions of
    parameter m between 0 and 1, and real u.  The returned functions are
    often written sn(u|m), cn(u|m), and dn(u|m).  The value of ph is such
    that if u = ellik(ph,m), then sn(u|m) = sin(ph) and cn(u|m) = cos(ph).
    """)

add_newdoc("scipy.special", "ellipkm1",
    """
    ellipkm1(p)

    The complete elliptic integral of the first kind around m=1.

    This function is defined as

    .. math:: K(p) = \\int_0^{\\pi/2} [1 - m \\sin(t)^2]^{-1/2} dt

    where `m = 1 - p`.

    Parameters
    ----------
    p : array_like
        Defines the parameter of the elliptic integral as m = 1 - p.

    Returns
    -------
    K : array_like
        Value of the elliptic integral.

    See Also
    --------
    ellipk

    """)

add_newdoc("scipy.special", "ellipkinc",
    """
    y=ellipkinc(phi,m) returns the incomplete elliptic integral of the first
    kind: integral(1/sqrt(1-m*sin(t)**2),t=0..phi)
    """)

add_newdoc("scipy.special", "erf",
    """
    erf(z)

    Returns the error function of complex argument.

    It is defined as ``2/sqrt(pi)*integral(exp(-t**2), t=0..z)``.

    Parameters
    ----------
    x : ndarray
        Input array.

    Returns
    -------
    res : ndarray
        The values of the error function at the given points x.

    See Also
    --------
    erfc, erfinv, erfcinv

    Notes
    -----
    The cumulative of the unit normal distribution is given by
    ``Phi(z) = 1/2[1 + erf(z/sqrt(2))]``.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Error_function
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover,
        1972. http://www.math.sfu.ca/~cbm/aands/page_297.htm
    .. [3] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    """)

add_newdoc("scipy.special", "erfc",
    """
    y=erfc(x) returns 1 - erf(x).

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    """)

add_newdoc("scipy.special", "erfi",
    """
    Imaginary error function, -i erf(i z)

    .. versionadded:: 0.12.0

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    """)

add_newdoc("scipy.special", "erfcx",
    """
    Scaled complementary error function, exp(x^2) erfc(x)

    .. versionadded:: 0.12.0

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    """)

add_newdoc("scipy.special", "eval_jacobi",
    """
    eval_jacobi(n, alpha, beta, x, out=None)

    Evaluate Jacobi polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_sh_jacobi",
    """
    eval_sh_jacobi(n, p, q, x, out=None)

    Evaluate shifted Jacobi polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_gegenbauer",
    """
    eval_gegenbauer(n, alpha, x, out=None)

    Evaluate Gegenbauer polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_chebyt",
    """
    eval_chebyt(n, x, out=None)

    Evaluate Chebyshev T polynomial at a point.

    This routine is numerically stable for `x` in ``[-1, 1]`` at least
    up to order ``10000``.
    """)

add_newdoc("scipy.special", "eval_chebyu",
    """
    eval_chebyu(n, x, out=None)

    Evaluate Chebyshev U polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_chebys",
    """
    eval_chebys(n, x, out=None)

    Evaluate Chebyshev S polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_chebyc",
    """
    eval_chebyc(n, x, out=None)

    Evaluate Chebyshev C polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_sh_chebyt",
    """
    eval_sh_chebyt(n, x, out=None)

    Evaluate shifted Chebyshev T polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_sh_chebyu",
    """
    eval_sh_chebyu(n, x, out=None)

    Evaluate shifted Chebyshev U polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_legendre",
    """
    eval_legendre(n, x, out=None)

    Evaluate Legendre polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_sh_legendre",
    """
    eval_sh_legendre(n, x, out=None)

    Evaluate shifted Legendre polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_genlaguerre",
    """
    eval_genlaguerre(n, alpha, x, out=None)

    Evaluate generalized Laguerre polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_laguerre",
     """
    eval_laguerre(n, x, out=None)

    Evaluate Laguerre polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_hermite",
    """
    eval_hermite(n, x, out=None)

    Evaluate Hermite polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_hermitenorm",
    """
    eval_hermitenorm(n, x, out=None)

    Evaluate normalized Hermite polynomial at a point.
    """)

add_newdoc("scipy.special", "exp1",
    """
    y=exp1(z) returns the exponential integral (n=1) of complex argument
    z: integral(exp(-z*t)/t,t=1..inf).
    """)

add_newdoc("scipy.special", "exp10",
    """
    y=exp10(x) returns 10 raised to the x power.
    """)

add_newdoc("scipy.special", "exp2",
    """
    y=exp2(x) returns 2 raised to the x power.
    """)

add_newdoc("scipy.special", "expi",
    """
    y=expi(x) returns an exponential integral of argument x defined as
    integral(exp(t)/t,t=-inf..x).  See expn for a different exponential
    integral.
    """)

add_newdoc('scipy.special', 'expit',
    """
    Expit ufunc for ndarrays.

    The expit function is defined as expit(x) = 1/(1+exp(-x)).
    Note that expit is the inverse logit function.

    Parameters
    ----------
    x : ndarray
        The ndarray to apply expit to element-wise.

    Returns
    -------
    out : ndarray
        An ndarray of the same shape as x. Its entries
        are expit of the corresponding entry of x.

    Notes
    -----
    As a ufunc logit takes a number of optional
    keywork arguments. For more information
    see `ufuncs <http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_
    """)

add_newdoc("scipy.special", "expm1",
    """
    y=expm1(x) calculates exp(x) - 1 for use when x is near zero.
    """)

add_newdoc("scipy.special", "expn",
    """
    y=expn(n,x) returns the exponential integral for integer n and
    non-negative x and n: integral(exp(-x*t) / t**n, t=1..inf).
    """)

add_newdoc("scipy.special", "fdtr",
    """
    y=fdtr(dfn,dfd,x) returns the area from zero to x under the F density
    function (also known as Snedcor's density or the variance ratio
    density).  This is the density of X = (unum/dfn)/(uden/dfd), where unum and
    uden are random variables having Chi square distributions with dfn and
    dfd degrees of freedom, respectively.
    """)

add_newdoc("scipy.special", "fdtrc",
    """
    y=fdtrc(dfn,dfd,x) returns the complemented F distribution function.
    """)

add_newdoc("scipy.special", "fdtri",
    """
    x=fdtri(dfn,dfd,p) finds the F density argument x such that
    fdtr(dfn,dfd,x)=p.
    """)

add_newdoc("scipy.special", "fdtridfd",
    """
    x=fdtridfd(dfn,p,x) finds the F density argument dfd such that
    fdtr(dfn,dfd,x)=p.
    """)

add_newdoc("scipy.special", "fdtridfn",
    """
    x=fdtridfn(p,dfd,x) finds the F density argument dfn such that
    fdtr(dfn,dfd,x)=p.
    """)

add_newdoc("scipy.special", "fresnel",
    """
    (ssa,cca)=fresnel(z) returns the Fresnel sin and cos integrals: integral(sin(pi/2
    * t**2),t=0..z) and integral(cos(pi/2 * t**2),t=0..z) for real or
    complex z.
    """)

add_newdoc("scipy.special", "gamma",
    """
    y=gamma(z) returns the gamma function of the argument.  The gamma
    function is often referred to as the generalized factorial since
    z*gamma(z) = gamma(z+1) and gamma(n+1) = n! for natural number n.
    """)

add_newdoc("scipy.special", "gammainc",
    """
    y=gammainc(a,x) returns the incomplete gamma integral defined as
    1 / gamma(a) * integral(exp(-t) * t**(a-1), t=0..x).  a must be
    positive and x must be >= 0.
    """)

add_newdoc("scipy.special", "gammaincc",
    """
    y=gammaincc(a,x) returns the complemented incomplete gamma integral
    defined as 1 / gamma(a) * integral(exp(-t) * t**(a-1), t=x..inf) = 1 -
    gammainc(a,x).  a must be positive and x must be >= 0.
    """)

add_newdoc("scipy.special", "gammainccinv",
    """
    x=gammainccinv(a,y) returns x such that gammaincc(a,x) = y.
    """)

add_newdoc("scipy.special", "gammaincinv",
    """
    gammaincinv(a, y) returns x such that gammainc(a, x) = y.
    """)

add_newdoc("scipy.special", "gammaln",
    """
    y=gammaln(z) returns the base e logarithm of the absolute value of the
    gamma function of z: ln(abs(gamma(z)))

    See Also
    --------
    gammasgn
    """)

add_newdoc("scipy.special", "gammasgn",
    """
    y=gammasgn(x) returns the sign of the gamma function.

    See Also
    --------
    gammaln
    """)

add_newdoc("scipy.special", "gdtr",
    """
    y=gdtr(a,b,x) returns the integral from zero to x of the gamma
    probability density function: a**b / gamma(b) * integral(t**(b-1) exp(-at),t=0..x).
    The arguments a and b are used differently here than in other definitions.
    """)

add_newdoc("scipy.special", "gdtrc",
    """
    y=gdtrc(a,b,x) returns the integral from x to infinity of the gamma
    probability density function.  SEE gdtr, gdtri
    """)

add_newdoc("scipy.special", "gdtria",
    """
    gdtria(p, b, x, out=None)

    Inverse with respect to `a` of `gdtr(a, b, x)`.

    `a = gdtria(p, b, x)` returns the inverse with respect to the parameter `a`
    of `p = gdtr(a, b, x)`, the cumulative distribution function of the gamma
    distribution.

    Parameters
    ----------
    p : array_like
        Probability values.
    b : array_like
        `b` parameter values of `gdtr(a, b, x)`.  `b` is the "shape" parameter
        of the gamma distribution.
    x : array_like
        Nonnegative real values, from the domain of the gamma distribution.
    out : ndarray, optional
        If a fourth argument is given, it must be a numpy.ndarray whose size
        matches the broadcast result of `a`, `b` and `x`.  `out` is then the
        array returned by the function.

    Returns
    -------
    a : ndarray
        Values of the `a` parameter such that `p = gdtr(a, b, x)`.  `1/a`
        is the "scale" parameter of the gamma distribution.

    See Also
    --------
    gdtr : CDF of the gamma distribution.
    gdtrib : Inverse with respect to `b` of `gdtr(a, b, x)`.
    gdtrix : Inverse with respect to `x` of `gdtr(a, b, x)`.

    Examples
    --------
    First evaluate `gdtr`.

    >>> p = gdtr(1.2, 3.4, 5.6)
    >>> print(p)
    0.94378087442

    Verify the inverse.

    >>> gdtria(p, 3.4, 5.6)
    1.2
    """)

add_newdoc("scipy.special", "gdtrib",
    """
    gdtrib(a, p, x, out=None)

    Inverse with respect to `b` of `gdtr(a, b, x)`.

    `b = gdtrib(a, p, x)` returns the inverse with respect to the parameter `b`
    of `p = gdtr(a, b, x)`, the cumulative distribution function of the gamma
    distribution.

    Parameters
    ----------
    a : array_like
        `a` parameter values of `gdtr(a, b, x)`. `1/a` is the "scale"
        parameter of the gamma distribution.
    p : array_like
        Probability values.
    x : array_like
        Nonnegative real values, from the domain of the gamma distribution.
    out : ndarray, optional
        If a fourth argument is given, it must be a numpy.ndarray whose size
        matches the broadcast result of `a`, `b` and `x`.  `out` is then the
        array returned by the function.

    Returns
    -------
    b : ndarray
        Values of the `b` parameter such that `p = gdtr(a, b, x)`.  `b` is
        the "shape" parameter of the gamma distribution.

    See Also
    --------
    gdtr : CDF of the gamma distribution.
    gdtria : Inverse with respect to `a` of `gdtr(a, b, x)`.
    gdtrix : Inverse with respect to `x` of `gdtr(a, b, x)`.

    Examples
    --------
    First evaluate `gdtr`.

    >>> p = gdtr(1.2, 3.4, 5.6)
    >>> print(p)
    0.94378087442

    Verify the inverse.

    >>> gdtrib(1.2, p, 5.6)
    3.3999999999723882
    """)

add_newdoc("scipy.special", "gdtrix",
    """
    gdtrix(a, b, p, out=None)

    Inverse with respect to `x` of `gdtr(a, b, x)`.

    `x = gdtrix(a, b, p)` returns the inverse with respect to the parameter `x`
    of `p = gdtr(a, b, x)`, the cumulative distribution function of the gamma
    distribution. This is also known as the p'th quantile of the distribution.

    Parameters
    ----------
    a : array_like
        `a` parameter values of `gdtr(a, b, x)`.  `1/a` is the "scale"
        parameter of the gamma distribution.
    b : array_like
        `b` parameter values of `gdtr(a, b, x)`.  `b` is the "shape" parameter
        of the gamma distribution.
    p : array_like
        Probability values.
    out : ndarray, optional
        If a fourth argument is given, it must be a numpy.ndarray whose size
        matches the broadcast result of `a`, `b` and `x`.  `out` is then the
        array returned by the function.

    Returns
    -------
    x : ndarray
        Values of the `x` parameter such that `p = gdtr(a, b, x)`.

    See Also
    --------
    gdtr : CDF of the gamma distribution.
    gdtria : Inverse with respect to `a` of `gdtr(a, b, x)`.
    gdtrib : Inverse with respect to `b` of `gdtr(a, b, x)`.

    Examples
    --------
    First evaluate `gdtr`.

    >>> p = gdtr(1.2, 3.4, 5.6)
    >>> print(p)
    0.94378087442

    Verify the inverse.

    >>> gdtrix(1.2, 3.4, p)
    5.5999999999999996
    """)

add_newdoc("scipy.special", "hankel1",
    """
    y=hankel1(v,z) returns the Hankel function of the first kind for real order v and complex argument z.
    """)

add_newdoc("scipy.special", "hankel1e",
    """
    y=hankel1e(v,z) returns the exponentially scaled Hankel function of the first
    kind for real order v and complex argument z:
    hankel1e(v,z) = hankel1(v,z) * exp(-1j * z)
    """)

add_newdoc("scipy.special", "hankel2",
    """
    y=hankel2(v,z) returns the Hankel function of the second kind for real order v and complex argument z.
    """)

add_newdoc("scipy.special", "hankel2e",
    """
    y=hankel2e(v,z) returns the exponentially scaled Hankel function of the second
    kind for real order v and complex argument z:
    hankel1e(v,z) = hankel1(v,z) * exp(1j * z)
    """)

add_newdoc("scipy.special", "hyp1f1",
    """
    y=hyp1f1(a,b,x) returns the confluent hypergeometeric function
    ( 1F1(a,b;x) ) evaluated at the values a, b, and x.
    """)

add_newdoc("scipy.special", "hyp1f2",
    """
    (y,err)=hyp1f2(a,b,c,x) returns (y,err) with the hypergeometric function 1F2 in y and an error estimate in err.
    """)

add_newdoc("scipy.special", "hyp2f0",
    """
    (y,err)=hyp2f0(a,b,x,type) returns (y,err) with the hypergeometric function 2F0 in y and an error estimate in err.  The input type determines a convergence factor and
    can be either 1 or 2.
    """)

add_newdoc("scipy.special", "hyp2f1",
    """
    y=hyp2f1(a,b,c,z) returns the Gauss hypergeometric function
    ( 2F1(a,b;c;z) ).
    """)

add_newdoc("scipy.special", "hyp3f0",
    """
    (y,err)=hyp3f0(a,b,c,x) returns (y,err) with the hypergeometric function 3F0 in y and an error estimate in err.
    """)

add_newdoc("scipy.special", "hyperu",
    """
    y=hyperu(a,b,x) returns the confluent hypergeometric function of the
    second kind U(a,b,x).
    """)

add_newdoc("scipy.special", "i0",
    """
    y=i0(x) returns the modified Bessel function of order 0 at x.
    """)

add_newdoc("scipy.special", "i0e",
    """
    y=i0e(x) returns the exponentially scaled modified Bessel function
    of order 0 at x.  i0e(x) = exp(-abs(x)) * i0(x).
    """)

add_newdoc("scipy.special", "i1",
    """
    y=i1(x) returns the modified Bessel function of order 1 at x.
    """)

add_newdoc("scipy.special", "i1e",
    """
    y=i1e(x) returns the exponentially scaled modified Bessel function
    of order 0 at x.  i1e(x) = exp(-abs(x)) * i1(x).
    """)

add_newdoc("scipy.special", "it2i0k0",
    """
    (ii0,ik0)=it2i0k0(x) returns the integrals int((i0(t)-1)/t,t=0..x) and
    int(k0(t)/t,t=x..infinitity).
    """)

add_newdoc("scipy.special", "it2j0y0",
    """
    (ij0,iy0)=it2j0y0(x) returns the integrals int((1-j0(t))/t,t=0..x) and
    int(y0(t)/t,t=x..infinitity).
    """)

add_newdoc("scipy.special", "it2struve0",
    """
    y=it2struve0(x) returns the integral of the Struve function of order 0
    divided by t from x to infinity:  integral(H0(t)/t, t=x..inf).
    """)

add_newdoc("scipy.special", "itairy",
    """
    (Apt,Bpt,Ant,Bnt)=itairy(x) calculates the integral of Airy functions from 0 to x
    for positive (Apt, Bpt) and negative (Ant, Bnt) arguments.
    """)

add_newdoc("scipy.special", "iti0k0",
    """
    (ii0,ik0)=iti0k0(x) returns simple integrals from 0 to x of the zeroth order
    modified Bessel functions i0 and k0.
    """)

add_newdoc("scipy.special", "itj0y0",
    """
    (ij0,iy0)=itj0y0(x) returns simple integrals from 0 to x of the zeroth order
    Bessel functions j0 and y0.
    """)

add_newdoc("scipy.special", "itmodstruve0",
    """
    y=itmodstruve0(x) returns the integral of the modified Struve function
    of order 0 from 0 to x:  integral(L0(t), t=0..x).
    """)

add_newdoc("scipy.special", "itstruve0",
    """
    y=itstruve0(x) returns the integral of the Struve function of order 0
    from 0 to x:  integral(H0(t), t=0..x).
    """)

add_newdoc("scipy.special", "iv",
    """
    y=iv(v,z) returns the modified Bessel function of real order v of
    z.  If z is of real type and negative, v must be integer valued.
    """)

add_newdoc("scipy.special", "ive",
    """
    y=ive(v,z) returns the exponentially scaled modified Bessel function of
    real order v and complex z: ive(v,z) = iv(v,z) * exp(-abs(z.real))
    """)

add_newdoc("scipy.special", "j0",
    """
    y=j0(x) returns the Bessel function of order 0 at x.
    """)

add_newdoc("scipy.special", "j1",
    """
    y=j1(x) returns the Bessel function of order 1 at x.
    """)

add_newdoc("scipy.special", "jn",
    """
    y=jn(n,x) returns the Bessel function of integer order n at  x.
    """)

add_newdoc("scipy.special", "jv",
    """
    y=jv(v,z) returns the Bessel function of real order v at complex z.
    """)

add_newdoc("scipy.special", "jve",
    """
    y=jve(v,z) returns the exponentially scaled Bessel function of real order
    v at complex z: jve(v,z) = jv(v,z) * exp(-abs(z.imag))
    """)

add_newdoc("scipy.special", "k0",
    """
    y=k0(x) returns the modified Bessel function of the second kind (sometimes called the third kind) of
    order 0 at x.
    """)

add_newdoc("scipy.special", "k0e",
    """
    y=k0e(x) returns the exponentially scaled modified Bessel function
    of the second kind (sometimes called the third kind) of order 0 at x.  k0e(x) = exp(x) * k0(x).
    """)

add_newdoc("scipy.special", "k1",
    """
    y=i1(x) returns the modified Bessel function of the second kind (sometimes called the third kind) of
    order 1 at x.
    """)

add_newdoc("scipy.special", "k1e",
    """
    y=k1e(x) returns the exponentially scaled modified Bessel function
    of the second kind (sometimes called the third kind) of order 1 at x.  k1e(x) = exp(x) * k1(x)
    """)

add_newdoc("scipy.special", "kei",
    """
    y=kei(x) returns the Kelvin function ker x
    """)

add_newdoc("scipy.special", "keip",
    """
    y=keip(x) returns the derivative of the Kelvin function kei x
    """)

add_newdoc("scipy.special", "kelvin",
    """
    (Be, Ke, Bep, Kep)=kelvin(x) returns the tuple (Be, Ke, Bep, Kep) which contains
    complex numbers representing the real and imaginary Kelvin functions
    and their derivatives evaluated at x.  For example,
    kelvin(x)[0].real = ber x and kelvin(x)[0].imag = bei x with similar
    relationships for ker and kei.
    """)

add_newdoc("scipy.special", "ker",
    """
    y=ker(x) returns the Kelvin function ker x
    """)

add_newdoc("scipy.special", "kerp",
    """
    y=kerp(x) returns the derivative of the Kelvin function ker x
    """)

add_newdoc("scipy.special", "kn",
    """
    y=kn(n,x) returns the modified Bessel function of the second kind (sometimes called the third kind) for
    integer order n at x.
    """)

add_newdoc("scipy.special", "kolmogi",
    """
    y=kolmogi(p) returns y such that kolmogorov(y) = p
    """)

add_newdoc("scipy.special", "kolmogorov",
    """
    p=kolmogorov(y) returns the complementary cumulative distribution
    function of Kolmogorov's limiting distribution (Kn* for large n)
    of a two-sided test for equality between an empirical and a theoretical
    distribution. It is equal to the (limit as n->infinity of the) probability
    that sqrt(n) * max absolute deviation > y.
    """)

add_newdoc("scipy.special", "kv",
    """
    y=kv(v,z) returns the modified Bessel function of the second kind (sometimes called the third kind) for
    real order v at complex z.
    """)

add_newdoc("scipy.special", "kve",
    """
    y=kve(v,z) returns the exponentially scaled, modified Bessel function
    of the second kind (sometimes called the third kind) for real order v at complex z: kve(v,z) = kv(v,z) * exp(z)
    """)

add_newdoc("scipy.special", "log1p",
    """
    y=log1p(x) calculates log(1+x) for use when x is near zero.
    """)

add_newdoc('scipy.special', 'logit',
    """
    Logit ufunc for ndarrays.

    The logit function is defined as logit(p) = log(p/(1-p)).
    Note that logit(0) = -inf, logit(1) = inf, and logit(p)
    for p<0 or p>1 yields nan.

    Parameters
    ----------
    x : ndarray
        The ndarray to apply logit to element-wise.

    Returns
    -------
    out : ndarray
        An ndarray of the same shape as x. Its entries
        are logit of the corresponding entry of x.

    Notes
    -----
    As a ufunc logit takes a number of optional
    keywork arguments. For more information
    see `ufuncs <http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_
    """)

add_newdoc("scipy.special", "lpmv",
    """
    y=lpmv(m,v,x) returns the associated legendre function of integer order
    m and real degree v (s.t. v>-m-1 or v<m): ``|x| <= 1``.
    """)

add_newdoc("scipy.special", "mathieu_a",
    """
    lmbda=mathieu_a(m,q) returns the characteristic value for the even solution,
    ce_m(z,q), of Mathieu's equation
    """)

add_newdoc("scipy.special", "mathieu_b",
    """
    lmbda=mathieu_b(m,q) returns the characteristic value for the odd solution,
    se_m(z,q), of Mathieu's equation
    """)

add_newdoc("scipy.special", "mathieu_cem",
    """
    (y,yp)=mathieu_cem(m,q,x) returns the even Mathieu function, ce_m(x,q),
    of order m and parameter q evaluated at x (given in degrees).
    Also returns the derivative with respect to x of ce_m(x,q)
    """)

add_newdoc("scipy.special", "mathieu_modcem1",
    """
    (y,yp)=mathieu_modcem1(m,q,x) evaluates the even modified Mathieu function
    of the first kind, Mc1m(x,q), and its derivative at x for order m and
    parameter q.
    """)

add_newdoc("scipy.special", "mathieu_modcem2",
    """
    (y,yp)=mathieu_modcem2(m,q,x) evaluates the even modified Mathieu function
    of the second kind, Mc2m(x,q), and its derivative at x (given in degrees)
    for order m and parameter q.
    """)

add_newdoc("scipy.special", "mathieu_modsem1",
    """
    (y,yp)=mathieu_modsem1(m,q,x) evaluates the odd modified Mathieu function
    of the first kind, Ms1m(x,q), and its derivative at x (given in degrees)
    for order m and parameter q.
    """)

add_newdoc("scipy.special", "mathieu_modsem2",
    """
    (y,yp)=mathieu_modsem2(m,q,x) evaluates the odd modified Mathieu function
    of the second kind, Ms2m(x,q), and its derivative at x (given in degrees)
    for order m and parameter q.
    """)

add_newdoc("scipy.special", "mathieu_sem",
    """
    (y,yp)=mathieu_sem(m,q,x) returns the odd Mathieu function, se_m(x,q),
    of order m and parameter q evaluated at x (given in degrees).
    Also returns the derivative with respect to x of se_m(x,q).
    """)

add_newdoc("scipy.special", "modfresnelm",
    """
    (fm,km)=modfresnelp(x) returns the modified Fresnel integrals ``F_-(x)`` and ``K_-(x)``
    as ``fp=integral(exp(-1j*t*t),t=x..inf)`` and ``kp=1/sqrt(pi)*exp(1j*(x*x+pi/4))*fp``
    """)

add_newdoc("scipy.special", "modfresnelp",
    """
    (fp,kp)=modfresnelp(x) returns the modified Fresnel integrals F_+(x) and K_+(x)
    as fp=integral(exp(1j*t*t),t=x..inf) and kp=1/sqrt(pi)*exp(-1j*(x*x+pi/4))*fp
    """)

add_newdoc("scipy.special", "modstruve",
    """
    y=modstruve(v,x) returns the modified Struve function Lv(x) of order
    v at x, x must be positive unless v is an integer and it is recommended
    that ``|v| <= 20``.
    """)

add_newdoc("scipy.special", "nbdtr",
    """
    y=nbdtr(k,n,p) returns the sum of the terms 0 through k of the
    negative binomial distribution: sum((n+j-1)Cj p**n (1-p)**j,j=0..k).
    In a sequence of Bernoulli trials this is the probability that k or
    fewer failures precede the nth success.
    """)

add_newdoc("scipy.special", "nbdtrc",
    """
    y=nbdtrc(k,n,p) returns the sum of the terms k+1 to infinity of the
    negative binomial distribution.
    """)

add_newdoc("scipy.special", "nbdtri",
    """
    p=nbdtri(k,n,y) finds the argument p such that nbdtr(k,n,p)=y.
    """)

add_newdoc("scipy.special", "nbdtrik",
    """
    k=nbdtrik(y,n,p) finds the argument k such that nbdtr(k,n,p)=y.
    """)

add_newdoc("scipy.special", "nbdtrin",
    """
    n=nbdtrin(k,y,p) finds the argument n such that nbdtr(k,n,p)=y.
    """)

add_newdoc("scipy.special", "ncfdtr",
    """
    """)

add_newdoc("scipy.special", "ncfdtri",
    """
    """)

add_newdoc("scipy.special", "ncfdtrifn",
    """
    """)

add_newdoc("scipy.special", "ncfdtridfd",
    """
    """)

add_newdoc("scipy.special", "ncfdtridfn",
    """
    """)

add_newdoc("scipy.special", "ncfdtrinc",
    """
    """)

add_newdoc("scipy.special", "nctdtr",
    """
    """)

add_newdoc("scipy.special", "nctdtridf",
    """
    """)

add_newdoc("scipy.special", "nctdtrinc",
    """
    """)

add_newdoc("scipy.special", "nctdtrit",
    """
    """)

add_newdoc("scipy.special", "ndtr",
    """
    y=ndtr(x) returns the area under the standard Gaussian probability
    density function, integrated from minus infinity to x:
    1/sqrt(2*pi) * integral(exp(-t**2 / 2),t=-inf..x)
    """)

add_newdoc("scipy.special", "nrdtrimn",
    """
    """)

add_newdoc("scipy.special", "nrdtrisd",
    """
    """)

add_newdoc("scipy.special", "log_ndtr",
    """
    y=log_ndtr(x) returns the log of the area under the standard Gaussian probability
    density function, integrated from minus infinity to x:
    1/sqrt(2*pi) * integral(exp(-t**2 / 2),t=-inf..x)
    """)

add_newdoc("scipy.special", "ndtri",
    """
    x=ndtri(y) returns the argument x for which the area udnder the
    Gaussian probability density function (integrated from minus infinity
    to x) is equal to y.
    """)

add_newdoc("scipy.special", "obl_ang1",
    """
    (s,sp)=obl_ang1(m,n,c,x) computes the oblate sheroidal angular function
    of the first kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``.
    """)

add_newdoc("scipy.special", "obl_ang1_cv",
    """
    (s,sp)=obl_ang1_cv(m,n,c,cv,x) computes the oblate sheroidal angular function
    of the first kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed
    characteristic value.
    """)

add_newdoc("scipy.special", "obl_cv",
    """
    cv=obl_cv(m,n,c) computes the characteristic value of oblate spheroidal
    wave functions of order m,n (n>=m) and spheroidal parameter c.
    """)

add_newdoc("scipy.special", "obl_rad1",
    """
    (s,sp)=obl_rad1(m,n,c,x) computes the oblate sheroidal radial function
    of the first kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``.
    """)

add_newdoc("scipy.special", "obl_rad1_cv",
    """
    (s,sp)=obl_rad1_cv(m,n,c,cv,x) computes the oblate sheroidal radial function
    of the first kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed
    characteristic value.
    """)

add_newdoc("scipy.special", "obl_rad2",
    """
    (s,sp)=obl_rad2(m,n,c,x) computes the oblate sheroidal radial function
    of the second kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``.
    """)

add_newdoc("scipy.special", "obl_rad2_cv",
    """
    (s,sp)=obl_rad2_cv(m,n,c,cv,x) computes the oblate sheroidal radial function
    of the second kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed
    characteristic value.
    """)

add_newdoc("scipy.special", "pbdv",
    """
    (d,dp)=pbdv(v,x) returns (d,dp) with the parabolic cylinder function Dv(x) in
    d and the derivative, Dv'(x) in dp.
    """)

add_newdoc("scipy.special", "pbvv",
    """
    (v,vp)=pbvv(v,x) returns (v,vp) with the parabolic cylinder function Vv(x) in
    v and the derivative, Vv'(x) in vp.
    """)

add_newdoc("scipy.special", "pbwa",
    """
    (w,wp)=pbwa(a,x) returns (w,wp) with the parabolic cylinder function W(a,x) in
    w and the derivative, W'(a,x) in wp.  May not be accurate for large (>5)
    arguments in a and/or x.
    """)

add_newdoc("scipy.special", "pdtr",
    """
    y=pdtr(k,m) returns the sum of the first k terms of the Poisson
    distribution: sum(exp(-m) * m**j / j!, j=0..k) = gammaincc( k+1, m).
    Arguments must both be positive and k an integer.
    """)

add_newdoc("scipy.special", "pdtrc",
    """
    y=pdtrc(k,m) returns the sum of the terms from k+1 to infinity of the
    Poisson distribution: sum(exp(-m) * m**j / j!, j=k+1..inf) = gammainc( k+1, m).
    Arguments must both be positive and k an integer.
    """)

add_newdoc("scipy.special", "pdtri",
    """
    m=pdtri(k,y) returns the Poisson variable m such that the sum
    from 0 to k of the Poisson density is equal to the given probability
    y:  calculated by gammaincinv( k+1, y).  k must be a nonnegative integer and
    y between 0 and 1.
    """)

add_newdoc("scipy.special", "pdtrik",
    """
    k=pdtrik(p,m) returns the quantile k such that pdtr(k,m)=p
    """)

add_newdoc("scipy.special", "pro_ang1",
    """
    (s,sp)=pro_ang1(m,n,c,x) computes the prolate sheroidal angular function
    of the first kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``.
    """)

add_newdoc("scipy.special", "pro_ang1_cv",
    """
    (s,sp)=pro_ang1_cv(m,n,c,cv,x) computes the prolate sheroidal angular function
    of the first kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed
    characteristic value.
    """)

add_newdoc("scipy.special", "pro_cv",
    """
    cv=pro_cv(m,n,c) computes the characteristic value of prolate spheroidal
    wave functions of order m,n (n>=m) and spheroidal parameter c.
    """)

add_newdoc("scipy.special", "pro_rad1",
    """
    (s,sp)=pro_rad1(m,n,c,x) computes the prolate sheroidal radial function
    of the first kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``.
    """)

add_newdoc("scipy.special", "pro_rad1_cv",
    """
    (s,sp)=pro_rad1_cv(m,n,c,cv,x) computes the prolate sheroidal radial function
    of the first kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed
    characteristic value.
    """)

add_newdoc("scipy.special", "pro_rad2",
    """
    (s,sp)=pro_rad2(m,n,c,x) computes the prolate sheroidal radial function
    of the second kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and |x|<1.0.
    """)

add_newdoc("scipy.special", "pro_rad2_cv",
    """
    (s,sp)=pro_rad2_cv(m,n,c,cv,x) computes the prolate sheroidal radial function
    of the second kind and its derivative (with respect to x) for mode parameters
    m>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed
    characteristic value.
    """)

add_newdoc("scipy.special", "psi",
    """
    y=psi(z) is the derivative of the logarithm of the gamma function
    evaluated at z (also called the digamma function).
    """)

add_newdoc("scipy.special", "radian",
    """
    y=radian(d,m,s) returns the angle given in (d)egrees, (m)inutes, and
    (s)econds in radians.
    """)

add_newdoc("scipy.special", "rgamma",
    """
    y=rgamma(z) returns one divided by the gamma function of x.
    """)

add_newdoc("scipy.special", "round",
    """
    y=Returns the nearest integer to x as a double precision
    floating point result.  If x ends in 0.5 exactly, the
    nearest even integer is chosen.
    """)

add_newdoc("scipy.special", "shichi",
    """
    (shi,chi)=shichi(x) returns the hyperbolic sine and cosine integrals:
    integral(sinh(t)/t,t=0..x) and eul + ln x +
    integral((cosh(t)-1)/t,t=0..x) where eul is Euler's Constant.
    """)

add_newdoc("scipy.special", "sici",
    """
    (si,ci)=sici(x) returns in si the integral of the sinc function from 0 to x:
    integral(sin(t)/t,t=0..x).  It returns in ci the cosine integral: eul + ln x +
    integral((cos(t) - 1)/t,t=0..x).
    """)

add_newdoc("scipy.special", "sindg",
    """
    y=sindg(x) calculates the sine of the angle x given in degrees.
    """)

add_newdoc("scipy.special", "smirnov",
    """
    y=smirnov(n,e) returns the exact Kolmogorov-Smirnov complementary
    cumulative distribution function (Dn+ or Dn-) for a one-sided test of
    equality between an empirical and a theoretical distribution. It is equal
    to the probability that the maximum difference between a theoretical
    distribution and an empirical one based on n samples is greater than e.
    """)

add_newdoc("scipy.special", "smirnovi",
    """
    e=smirnovi(n,y) returns e such that smirnov(n,e) = y.
    """)

add_newdoc("scipy.special", "spence",
    """
    y=spence(x) returns the dilogarithm integral: -integral(log t /
    (t-1),t=1..x)
    """)

add_newdoc("scipy.special", "stdtr",
    """
    p=stdtr(df,t) returns the integral from minus infinity to t of the Student t
    distribution with df > 0 degrees of freedom:
    gamma((df+1)/2)/(sqrt(df*pi)*gamma(df/2)) * integral((1+x**2/df)**(-df/2-1/2),
    x=-inf..t)
    """)

add_newdoc("scipy.special", "stdtridf",
    """
    t=stdtridf(p,t) returns the argument df such that stdtr(df,t) is equal to p.
    """)

add_newdoc("scipy.special", "stdtrit",
    """
    t=stdtrit(df,p) returns the argument t such that stdtr(df,t) is equal to p.
    """)

add_newdoc("scipy.special", "struve",
    """
    y=struve(v,x) returns the Struve function Hv(x) of order v at x, x
    must be positive unless v is an integer.
    """)

add_newdoc("scipy.special", "tandg",
    """
    y=tandg(x) calculates the tangent of the angle x given in degrees.
    """)

add_newdoc("scipy.special", "tklmbda",
    """
    """)

add_newdoc("scipy.special", "wofz",
    """
    y=wofz(z) returns the value of the fadeeva function for complex argument
    z: exp(-z**2)*erfc(-i*z)

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva
    """)

add_newdoc("scipy.special", "xlogy",
    """
    xlogy(x, y)

    Compute ``x*log(y)`` so that the result is 0 if `x = 0`.

    Parameters
    ----------
    x : array_like
        Multiplier
    y : array_like
        Argument

    Returns
    -------
    z : array_like
        Computed x*log(y)

    """)

add_newdoc("scipy.special", "xlog1py",
    """
    xlog1py(x, y)

    Compute ``x*log1p(y)`` so that the result is 0 if `x = 0`.

    Parameters
    ----------
    x : array_like
        Multiplier
    y : array_like
        Argument

    Returns
    -------
    z : array_like
        Computed x*log1p(y)

    """)

add_newdoc("scipy.special", "y0",
    """
    y=y0(x) returns the Bessel function of the second kind of order 0 at x.
    """)

add_newdoc("scipy.special", "y1",
    """
    y=y1(x) returns the Bessel function of the second kind of order 1 at x.
    """)

add_newdoc("scipy.special", "yn",
    """
    y=yn(n,x) returns the Bessel function of the second kind of integer
    order n at x.
    """)

add_newdoc("scipy.special", "yv",
    """
    y=yv(v,z) returns the Bessel function of the second kind of real
    order v at complex z.
    """)

add_newdoc("scipy.special", "yve",
    """
    y=yve(v,z) returns the exponentially scaled Bessel function of the second
    kind of real order v at complex z: yve(v,z) = yv(v,z) * exp(-abs(z.imag))
    """)

add_newdoc("scipy.special", "zeta",
    """
    y=zeta(x,q) returns the Riemann zeta function of two arguments:
    sum((k+q)**(-x),k=0..inf)
    """)

add_newdoc("scipy.special", "zetac",
    """
    y=zetac(x) returns 1.0 - the Riemann zeta function: sum(k**(-x), k=2..inf)
    """)

add_newdoc("scipy.special", "_struve_asymp_large_z",
    """
    Function for testing struve & modstruve
    """)

add_newdoc("scipy.special", "_struve_power_series",
    """
    Function for testing struve & modstruve
    """)

add_newdoc("scipy.special", "_struve_bessel_series",
    """
    Function for testing struve & modstruve
    """)

