"""
========================================
Special functions (:mod:`scipy.special`)
========================================

.. module:: scipy.special

Nearly all of the functions below are universal functions and follow
broadcasting and automatic array-looping rules. Exceptions are noted.

Error handling
==============

Errors are handled by returning nans, or other appropriate values.
Some of the special function routines will emit warnings when an error
occurs.  By default this is disabled.  To enable such messages use
``errprint(1)``, and to disable such messages use ``errprint(0)``.

Example:

    >>> print scipy.special.bdtr(-1,10,0.3)
    >>> scipy.special.errprint(1)
    >>> print scipy.special.bdtr(-1,10,0.3)

.. autosummary::
   :toctree: generated/

   errprint
   SpecialFunctionWarning -- Warning that can be issued with ``errprint(True)``

Available functions
===================

Airy functions
--------------

.. autosummary::
   :toctree: generated/

   airy     -- Airy functions and their derivatives.
   airye    -- Exponentially scaled Airy functions
   ai_zeros -- [+]Zeros of Airy functions Ai(x) and Ai'(x)
   bi_zeros -- [+]Zeros of Airy functions Bi(x) and Bi'(x)


Elliptic Functions and Integrals
--------------------------------

.. autosummary::
   :toctree: generated/

   ellipj    -- Jacobian elliptic functions
   ellipk    -- Complete elliptic integral of the first kind.
   ellipkm1  -- ellipkm1(x) == ellipk(1 - x)
   ellipkinc -- Incomplete elliptic integral of the first kind.
   ellipe    -- Complete elliptic integral of the second kind.
   ellipeinc -- Incomplete elliptic integral of the second kind.

Bessel Functions
----------------

.. autosummary::
   :toctree: generated/

   jv       -- Bessel function of real-valued order and complex argument.
   jve      -- Exponentially scaled Bessel function.
   yn       -- Bessel function of second kind (integer order).
   yv       -- Bessel function of the second kind (real-valued order).
   yve      -- Exponentially scaled Bessel function of the second kind.
   kn       -- Modified Bessel function of the second kind (integer order).
   kv       -- Modified Bessel function of the second kind (real order).
   kve      -- Exponentially scaled modified Bessel function of the second kind.
   iv       -- Modified Bessel function.
   ive      -- Exponentially scaled modified Bessel function.
   hankel1  -- Hankel function of the first kind.
   hankel1e -- Exponentially scaled Hankel function of the first kind.
   hankel2  -- Hankel function of the second kind.
   hankel2e -- Exponentially scaled Hankel function of the second kind.

The following is not an universal function:

.. autosummary::
   :toctree: generated/

   lmbda       -- [+]Sequence of lambda functions with arbitrary order v.

Zeros of Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^^

These are not universal functions:

.. autosummary::
   :toctree: generated/

   jnjnp_zeros -- [+]Zeros of integer-order Bessel functions and derivatives sorted in order.
   jnyn_zeros  -- [+]Zeros of integer-order Bessel functions and derivatives as separate arrays.
   jn_zeros    -- [+]Zeros of Jn(x)
   jnp_zeros   -- [+]Zeros of Jn'(x)
   yn_zeros    -- [+]Zeros of Yn(x)
   ynp_zeros   -- [+]Zeros of Yn'(x)
   y0_zeros    -- [+]Complex zeros: Y0(z0)=0 and values of Y0'(z0)
   y1_zeros    -- [+]Complex zeros: Y1(z1)=0 and values of Y1'(z1)
   y1p_zeros   -- [+]Complex zeros of Y1'(z1')=0 and values of Y1(z1')

Faster versions of common Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   j0       -- Bessel function of order 0.
   j1       -- Bessel function of order 1.
   y0       -- Bessel function of second kind of order 0.
   y1       -- Bessel function of second kind of order 1.
   i0       -- Modified Bessel function of order 0.
   i0e      -- Exponentially scaled modified Bessel function of order 0.
   i1       -- Modified Bessel function of order 1.
   i1e      -- Exponentially scaled modified Bessel function of order 1.
   k0       -- Modified Bessel function of the second kind of order 0.
   k0e      -- Exponentially scaled modified Bessel function of the second kind of order 0.
   k1       -- Modified Bessel function of the second kind of order 1.
   k1e      -- Exponentially scaled modified Bessel function of the second kind of order 1.

Integrals of Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   itj0y0     -- Basic integrals of j0 and y0 from 0 to x.
   it2j0y0    -- Integrals of (1-j0(t))/t from 0 to x and y0(t)/t from x to inf.
   iti0k0     -- Basic integrals of i0 and k0 from 0 to x.
   it2i0k0    -- Integrals of (i0(t)-1)/t from 0 to x and k0(t)/t from x to inf.
   besselpoly -- Integral of a Bessel function: Jv(2* a* x) * x[+]lambda from x=0 to 1.

Derivatives of Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   jvp     -- Nth derivative of Jv(v,z)
   yvp     -- Nth derivative of Yv(v,z)
   kvp     -- Nth derivative of Kv(v,z)
   ivp     -- Nth derivative of Iv(v,z)
   h1vp    -- Nth derivative of H1v(v,z)
   h2vp    -- Nth derivative of H2v(v,z)

Spherical Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^

These are not universal functions:

.. autosummary::
   :toctree: generated/

   sph_jn   -- [+]Sequence of spherical Bessel functions, jn(z)
   sph_yn   -- [+]Sequence of spherical Bessel functions, yn(z)
   sph_jnyn -- [+]Sequence of spherical Bessel functions, jn(z) and yn(z)
   sph_in   -- [+]Sequence of spherical Bessel functions, in(z)
   sph_kn   -- [+]Sequence of spherical Bessel functions, kn(z)
   sph_inkn -- [+]Sequence of spherical Bessel functions, in(z) and kn(z)

Riccati-Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^

These are not universal functions:

.. autosummary::
   :toctree: generated/

   riccati_jn -- [+]Sequence of Ricatti-Bessel functions of first kind.
   riccati_yn -- [+]Sequence of Ricatti-Bessel functions of second kind.

Struve Functions
----------------

.. autosummary::
   :toctree: generated/

   struve       -- Struve function --- Hv(x)
   modstruve    -- Modified Struve function --- Lv(x)
   itstruve0    -- Integral of H0(t) from 0 to x
   it2struve0   -- Integral of H0(t)/t from x to Inf.
   itmodstruve0 -- Integral of L0(t) from 0 to x.


Raw Statistical Functions
-------------------------

.. seealso:: :mod:`scipy.stats`: Friendly versions of these functions.

.. autosummary::
   :toctree: generated/

   bdtr       -- Sum of terms 0 through k of the binomial pdf.
   bdtrc      -- Sum of terms k+1 through n of the binomial pdf.
   bdtri      -- Inverse of bdtr
   btdtr      -- Integral from 0 to x of beta pdf.
   btdtri     -- Quantiles of beta distribution
   fdtr       -- Integral from 0 to x of F pdf.
   fdtrc      -- Integral from x to infinity under F pdf.
   fdtri      -- Inverse of fdtrc
   gdtr       -- Integral from 0 to x of gamma pdf.
   gdtrc      -- Integral from x to infinity under gamma pdf.
   gdtria     -- Inverse with respect to `a` of gdtr.
   gdtrib     -- Inverse with respect to `b` of gdtr.
   gdtrix     -- Inverse with respect to `x` of gdtr.
   nbdtr      -- Sum of terms 0 through k of the negative binomial pdf.
   nbdtrc     -- Sum of terms k+1 to infinity under negative binomial pdf.
   nbdtri     -- Inverse of nbdtr
   ncfdtr     -- CDF of non-central t distribution.
   ncfdtridfd -- Find degrees of freedom (denominator) of noncentral F distribution.
   ncfdtridfn -- Find degrees of freedom (numerator) of noncentral F distribution.
   ncfdtri    -- Inverse CDF of noncentral F distribution.
   ncfdtrinc  -- Find noncentrality parameter of noncentral F distribution.
   nctdtr     -- CDF of noncentral t distribution.
   nctdtridf  -- Find degrees of freedom of noncentral t distribution.
   nctdtrit   -- Inverse CDF of noncentral t distribution.
   nctdtrinc  -- Find noncentrality parameter of noncentral t distribution.
   nrdtrimn   -- Find mean of normal distribution from cdf and std.
   nrdtrisd   -- Find std of normal distribution from cdf and mean.
   pdtr       -- Sum of terms 0 through k of the Poisson pdf.
   pdtrc      -- Sum of terms k+1 to infinity of the Poisson pdf.
   pdtri      -- Inverse of pdtr
   stdtr      -- Integral from -infinity to t of the Student-t pdf.
   stdtridf   --
   stdtrit    --
   chdtr      -- Integral from 0 to x of the Chi-square pdf.
   chdtrc     -- Integral from x to infnity of Chi-square pdf.
   chdtri     -- Inverse of chdtrc.
   ndtr       -- Integral from -infinity to x of standard normal pdf
   ndtri      -- Inverse of ndtr (quantiles)
   smirnov    -- Kolmogorov-Smirnov complementary CDF for one-sided test statistic (Dn+ or Dn-)
   smirnovi   -- Inverse of smirnov.
   kolmogorov -- The complementary CDF of the (scaled) two-sided test statistic (Kn*) valid for large n.
   kolmogi    -- Inverse of kolmogorov
   tklmbda    -- Tukey-Lambda CDF
   logit      --
   expit      --
   boxcox     -- Compute the Box-Cox transformation.
   boxcox1p   -- Compute the Box-Cox transformation.


Information Theory Functions
----------------------------

.. autosummary::
   :toctree: generated/

   entr         -- entr(x) = -x*log(x)
   rel_entr     -- rel_entr(x, y) = x*log(x/y)
   kl_div       -- kl_div(x, y) = x*log(x/y) - x + y
   huber        -- Huber loss function.
   pseudo_huber -- Pseudo-Huber loss function.


Gamma and Related Functions
---------------------------

.. autosummary::
   :toctree: generated/

   gamma        -- Gamma function.
   gammaln      -- Log of the absolute value of the gamma function.
   gammasgn     -- Sign of the gamma function.
   gammainc     -- Incomplete gamma integral.
   gammaincinv  -- Inverse of gammainc.
   gammaincc    -- Complemented incomplete gamma integral.
   gammainccinv -- Inverse of gammaincc.
   beta         -- Beta function.
   betaln       -- Log of the absolute value of the beta function.
   betainc      -- Incomplete beta integral.
   betaincinv   -- Inverse of betainc.
   psi          -- Logarithmic derivative of the gamma function.
   rgamma       -- One divided by the gamma function.
   polygamma    -- Nth derivative of psi function.
   multigammaln -- Log of the multivariate gamma.
   digamma      -- Digamma function (derivative of the logarithm of gamma).


Error Function and Fresnel Integrals
------------------------------------

.. autosummary::
   :toctree: generated/

   erf           -- Error function.
   erfc          -- Complemented error function (1- erf(x))
   erfcx         -- Scaled complemented error function exp(x**2)*erfc(x)
   erfi          -- Imaginary error function, -i erf(i x)
   erfinv        -- Inverse of error function
   erfcinv       -- Inverse of erfc
   wofz          -- Fadeeva function.
   dawsn         -- Dawson's integral.
   fresnel       -- Fresnel sine and cosine integrals.
   fresnel_zeros -- Complex zeros of both Fresnel integrals
   modfresnelp   -- Modified Fresnel integrals F_+(x) and K_+(x)
   modfresnelm   -- Modified Fresnel integrals F_-(x) and K_-(x)

These are not universal functions:

.. autosummary::
   :toctree: generated/

   erf_zeros      -- [+]Complex zeros of erf(z)
   fresnelc_zeros -- [+]Complex zeros of Fresnel cosine integrals
   fresnels_zeros -- [+]Complex zeros of Fresnel sine integrals

Legendre Functions
------------------

.. autosummary::
   :toctree: generated/

   lpmv     -- Associated Legendre Function of arbitrary non-negative degree v.
   sph_harm -- Spherical Harmonics (complex-valued) Y^m_n(theta,phi)

These are not universal functions:

.. autosummary::
   :toctree: generated/

   clpmn    -- [+]Associated Legendre Function of the first kind for complex arguments.
   lpn      -- [+]Legendre Functions (polynomials) of the first kind
   lqn      -- [+]Legendre Functions of the second kind.
   lpmn     -- [+]Associated Legendre Function of the first kind for real arguments.
   lqmn     -- [+]Associated Legendre Function of the second kind.

Orthogonal polynomials
----------------------

The following functions evaluate values of orthogonal polynomials:

.. autosummary::
   :toctree: generated/

   assoc_laguerre
   eval_legendre
   eval_chebyt
   eval_chebyu
   eval_chebyc
   eval_chebys
   eval_jacobi
   eval_laguerre
   eval_genlaguerre
   eval_hermite
   eval_hermitenorm
   eval_gegenbauer
   eval_sh_legendre
   eval_sh_chebyt
   eval_sh_chebyu
   eval_sh_jacobi

The functions below, in turn, return the polynomial coefficients in
:ref:`orthopoly1d` objects, which function similarly as :ref:`numpy.poly1d`.
The :ref:`orthopoly1d` class also has an attribute ``weights`` which returns
the roots, weights, and total weights for the appropriate form of Gaussian
quadrature.  These are returned in an ``n x 3`` array with roots in the first
column, weights in the second column, and total weights in the final column.
Note that ``orthopoly1d`` objects are converted to ``poly1d`` when doing
arithmetic, and lose information of the original orthogonal polynomial.

.. autosummary::
   :toctree: generated/

   legendre    -- [+]Legendre polynomial P_n(x) (lpn -- for function).
   chebyt      -- [+]Chebyshev polynomial T_n(x)
   chebyu      -- [+]Chebyshev polynomial U_n(x)
   chebyc      -- [+]Chebyshev polynomial C_n(x)
   chebys      -- [+]Chebyshev polynomial S_n(x)
   jacobi      -- [+]Jacobi polynomial P^(alpha,beta)_n(x)
   laguerre    -- [+]Laguerre polynomial, L_n(x)
   genlaguerre -- [+]Generalized (Associated) Laguerre polynomial, L^alpha_n(x)
   hermite     -- [+]Hermite polynomial H_n(x)
   hermitenorm -- [+]Normalized Hermite polynomial, He_n(x)
   gegenbauer  -- [+]Gegenbauer (Ultraspherical) polynomials, C^(alpha)_n(x)
   sh_legendre -- [+]shifted Legendre polynomial, P*_n(x)
   sh_chebyt   -- [+]shifted Chebyshev polynomial, T*_n(x)
   sh_chebyu   -- [+]shifted Chebyshev polynomial, U*_n(x)
   sh_jacobi   -- [+]shifted Jacobi polynomial, J*_n(x) = G^(p,q)_n(x)

.. warning::

   Computing values of high-order polynomials (around ``order > 20``) using
   polynomial coefficients is numerically unstable. To evaluate polynomial
   values, the ``eval_*`` functions should be used instead.



Hypergeometric Functions
------------------------

.. autosummary::
   :toctree: generated/

   hyp2f1   -- Gauss hypergeometric function (2F1)
   hyp1f1   -- Confluent hypergeometric function (1F1)
   hyperu   -- Confluent hypergeometric function (U)
   hyp0f1   -- Confluent hypergeometric limit function (0F1)
   hyp2f0   -- Hypergeometric function (2F0)
   hyp1f2   -- Hypergeometric function (1F2)
   hyp3f0   -- Hypergeometric function (3F0)


Parabolic Cylinder Functions
----------------------------

.. autosummary::
   :toctree: generated/

   pbdv     -- Parabolic cylinder function Dv(x) and derivative.
   pbvv     -- Parabolic cylinder function Vv(x) and derivative.
   pbwa     -- Parabolic cylinder function W(a,x) and derivative.

These are not universal functions:

.. autosummary::
   :toctree: generated/

   pbdv_seq -- [+]Sequence of parabolic cylinder functions Dv(x)
   pbvv_seq -- [+]Sequence of parabolic cylinder functions Vv(x)
   pbdn_seq -- [+]Sequence of parabolic cylinder functions Dn(z), complex z

Mathieu and Related Functions
-----------------------------

.. autosummary::
   :toctree: generated/

   mathieu_a       -- Characteristic values for even solution (ce_m)
   mathieu_b       -- Characteristic values for odd solution (se_m)

These are not universal functions:

.. autosummary::
   :toctree: generated/

   mathieu_even_coef -- [+]sequence of expansion coefficients for even solution
   mathieu_odd_coef  -- [+]sequence of expansion coefficients for odd solution

The following return both function and first derivative:

.. autosummary::
   :toctree: generated/

   mathieu_cem     -- Even Mathieu function
   mathieu_sem     -- Odd Mathieu function
   mathieu_modcem1 -- Even modified Mathieu function of the first kind
   mathieu_modcem2 -- Even modified Mathieu function of the second kind
   mathieu_modsem1 -- Odd modified Mathieu function of the first kind
   mathieu_modsem2 -- Odd modified Mathieu function of the second kind

Spheroidal Wave Functions
-------------------------

.. autosummary::
   :toctree: generated/

   pro_ang1   -- Prolate spheroidal angular function of the first kind
   pro_rad1   -- Prolate spheroidal radial function of the first kind
   pro_rad2   -- Prolate spheroidal radial function of the second kind
   obl_ang1   -- Oblate spheroidal angular function of the first kind
   obl_rad1   -- Oblate spheroidal radial function of the first kind
   obl_rad2   -- Oblate spheroidal radial function of the second kind
   pro_cv     -- Compute characteristic value for prolate functions
   obl_cv     -- Compute characteristic value for oblate functions
   pro_cv_seq -- Compute sequence of prolate characteristic values
   obl_cv_seq -- Compute sequence of oblate characteristic values

The following functions require pre-computed characteristic value:

.. autosummary::
   :toctree: generated/

   pro_ang1_cv -- Prolate spheroidal angular function of the first kind
   pro_rad1_cv -- Prolate spheroidal radial function of the first kind
   pro_rad2_cv -- Prolate spheroidal radial function of the second kind
   obl_ang1_cv -- Oblate spheroidal angular function of the first kind
   obl_rad1_cv -- Oblate spheroidal radial function of the first kind
   obl_rad2_cv -- Oblate spheroidal radial function of the second kind

Kelvin Functions
----------------

.. autosummary::
   :toctree: generated/

   kelvin       -- All Kelvin functions (order 0) and derivatives.
   kelvin_zeros -- [+]Zeros of All Kelvin functions (order 0) and derivatives
   ber          -- Kelvin function ber x
   bei          -- Kelvin function bei x
   berp         -- Derivative of Kelvin function ber x
   beip         -- Derivative of Kelvin function bei x
   ker          -- Kelvin function ker x
   kei          -- Kelvin function kei x
   kerp         -- Derivative of Kelvin function ker x
   keip         -- Derivative of Kelvin function kei x

These are not universal functions:

.. autosummary::
   :toctree: generated/

   ber_zeros    -- [+]Zeros of Kelvin function bei x
   bei_zeros    -- [+]Zeros of Kelvin function ber x
   berp_zeros   -- [+]Zeros of derivative of Kelvin function ber x
   beip_zeros   -- [+]Zeros of derivative of Kelvin function bei x
   ker_zeros    -- [+]Zeros of Kelvin function kei x
   kei_zeros    -- [+]Zeros of Kelvin function ker x
   kerp_zeros   -- [+]Zeros of derivative of Kelvin function ker x
   keip_zeros   -- [+]Zeros of derivative of Kelvin function kei x

Combinatorics
-------------

.. autosummary::
    :toctree: generated/

    comb    -- [+]Combinations of N things taken k at a time, "N choose k"
    perm    -- [+]Permutations of N things taken k at a time, "k-permutations of N"

Other Special Functions
-----------------------

.. autosummary::
   :toctree: generated/

   agm          -- Arithmetic-Geometric Mean
   bernoulli    -- Bernoulli numbers
   binom        -- Binomial coefficient.
   diric        -- Dirichlet function (periodic sinc)
   euler        -- Euler numbers
   expn         -- Exponential integral.
   exp1         -- Exponential integral of order 1 (for complex argument)
   expi         -- Another exponential integral -- Ei(x)
   factorial    -- The factorial function, n! = special.gamma(n+1)
   factorial2   -- Double factorial, (n!)!
   factorialk   -- [+](...((n!)!)!...)! where there are k '!'
   shichi       -- Hyperbolic sine and cosine integrals.
   sici         -- Integral of the sinc and "cosinc" functions.
   spence       -- Dilogarithm integral.
   lambertw     -- Lambert W function
   zeta         -- Riemann zeta function of two arguments.
   zetac        -- Standard Riemann zeta function minus 1.

Convenience Functions
---------------------

.. autosummary::
   :toctree: generated/

   cbrt     -- Cube root.
   exp10    -- 10 raised to the x power.
   exp2     -- 2 raised to the x power.
   radian   -- radian angle given degrees, minutes, and seconds.
   cosdg    -- cosine of the angle given in degrees.
   sindg    -- sine of the angle given in degrees.
   tandg    -- tangent of the angle given in degrees.
   cotdg    -- cotangent of the angle given in degrees.
   log1p    -- log(1+x)
   expm1    -- exp(x)-1
   cosm1    -- cos(x)-1
   round    -- round the argument to the nearest integer. If argument ends in 0.5 exactly, pick the nearest even integer.
   xlogy    -- x*log(y)
   xlog1py  -- x*log1p(y)

.. [+] in the description indicates a function which is not a universal
.. function and does not follow broadcasting and automatic
.. array-looping rules.

"""

from __future__ import division, print_function, absolute_import

from ._ufuncs import *

from .basic import *
from . import specfun
from . import orthogonal
from .orthogonal import *
from .spfun_stats import multigammaln
from ._ellip_harm import ellip_harm, ellip_harm_2, ellip_normal
from .lambertw import lambertw


__all__ = [s for s in dir() if not s.startswith('_')]

from numpy.dual import register_func
register_func('i0',i0)
del register_func

from numpy.testing import Tester
test = Tester().test
