"""
 Special Functions

  general_function -- a class that wraps a Python function taking scalar
                         arguments into a generalized function which
                         can handle arrays of arguments using the broadcast
                         rules of Numeric Python.

  Airy Functions

    airy     -- Airy functions and their derivatives.
    airye    -- Exponentially scaled Airy functions

  Elliptic Functions and Integrals

    ellipj    -- Jacobian elliptic functions
    ellipk    -- Complete elliptic integral of the first kind.
    ellipkinc -- Incomplete elliptic integral of the first kind.
    ellipe    -- Complete elliptic integral of the second kind.
    ellipeinc -- Incomplete elliptic integral of the second kind.

  Bessel Functions

    jn       -- Bessel function of integer order and real argument.
    jv       -- Bessel function of real-valued order and complex argument.
    jve      -- Exponentially scaled Bessel function.
    yn       -- Bessel function of second kind (integer order).
    yv       -- Bessel function of the second kind (real-valued order).
    yve      -- Exponentially scaled Bessel function of the second kind.
    kn       -- Modified Bessel function of the third kind (integer order).
    kv       -- Modified Bessel function of the third kind (real order).
    kve      -- Exponentially scaled modified Bessel function of the
                  third kind.
    iv       -- Modified Bessel function.
    ive      -- Exponentially scaled modified Bessel function.
    hankel1  -- Hankel function of the first kind.
    hankel1e -- Exponentially scaled Hankel function of the first kind.
    hankel2  -- Hankel function of the second kind.
    hankel2e -- Exponentially scaled Hankel function of the second kind.

  Faster versions of common Bessel Functions.
    
    j0       -- Bessel function of order 0.
    j1       -- Bessel function of order 1.
    y0       -- Bessel function of second kind of order 0.
    y1       -- Bessel function of second kind of order 1.
    i0       -- Modified Bessel function of order 0.
    i0e      -- Exponentially scaled modified Bessel function of order 0.
    i1       -- Modified Bessel function of order 1.
    i1e      -- Exponentially scaled modified Bessel function of order 1.
    k0       -- Modified Bessel function of the third kind of order 0.
    k0e      -- Exponentially scaled modified Bessel function of the
                  third kind of order 0.
    k1       -- Modified Bessel function of the third kind of order 1.
    k1e      -- Exponentially scaled modified Bessel function of the
                  third kind of order 1.

  Integrals of Bessel Functions.

    itj0y0   -- Basic integrals of j0 and y0 from 0 to x.
    it2j0y0  -- Integrals of (1-j0(t))/t from 0 to x and
                  y0(t)/t from x to inf.
    iti0k0   -- Basic integrals of i0 and k0 from 0 to x.
    it2i0k0  -- Integrals of (i0(t)-1)/t from 0 to x and
                  k0(t)/t from x to inf.

  Derivatives of Bessel Functions.

    djv     -- Nth derivative of Jv(v,z)
    dyv     -- Nth derivative of Yv(v,z)
    dkv     -- Nth derivative of Kv(v,z)
    div     -- Nth derivative of Iv(v,z)
    dh1v    -- Nth derivative of H1v(v,z)
    dh2v    -- Nth derivative of H2v(v,z)
                  
  Raw Statistical Functions (Friendly versions in scipy.stats)

    bdtr       -- Sum of terms 0 through k of of the binomial pdf.
    bdtrc      -- Sum of terms k+1 through n of the binomial pdf.
    bdtri      -- Inverse of bdtr
    btdtr      -- Integral from 0 to x of beta pdf.
    btdtri     -- Quantiles of beta distribution
    fdtr       -- Integral from 0 to x of F pdf.
    fdtrc      -- Integral from x to infinity under F pdf.
    fdtri      -- Inverse of fdtrc
    gdtr       -- Integral from 0 to x of gamma pdf.
    gdtrc      -- Integral from x to infinity under gamma pdf.
    gdtri      -- Quantiles of gamma distribution
    nbdtr      -- Sum of terms 0 through k of the negative binomial pdf.
    nbdtrc     -- Sum of terms k+1 to infinity under negative binomial pdf.
    nbdtri     -- Inverse of nbdtr
    pdtr       -- Sum of terms 0 through k of the Poisson pdf.
    pdtrc      -- Sum of terms k+1 to infinity of the Poisson pdf.
    pdtri      -- Inverse of pdtr
    stdtr      -- Integral from -infinity to t of the Student-t pdf.
    stdtri     -- Inverse of stdtr (quantiles)
    chdtr      -- Integral from 0 to x of the Chi-square pdf.
    chdtrc     -- Integral from x to infnity of Chi-square pdf.
    chdtri     -- Inverse of chdtrc.
    ndtr       -- Integral from -infinity to x of standard normal pdf
    ndtri      -- Inverse of ndtr (quantiles)
    smirnov    -- Kolmogorov-Smirnov complementary CDF for one-sided
                    test statistic (Dn+ or Dn-)
    smirnovi   -- Inverse of smirnov.
    kolmogorov -- The complementary CDF of the (scaled) two-sided test
                          statistic (Kn*) valid for large n.
    kolmogi    -- Inverse of kolmogorov
    tklmbda    -- Tukey-Lambda CDF
      
  Gamma and Related Functions

    gamma        -- Gamma function.
    gammaln      -- Log of the absolute value of the gamma function.
    gammainc     -- Incomplete gamma integral.
    gammaincinv  -- Inverse of gammainc.
    gammaincc    -- Complemented incomplete gamma integral.
    gammainccinv -- Inverse of gammaincc.
    beta         -- Beta function.
    betaln       -- Log of the absolute value of the beta function.
    betainc      -- Incomplete beta integral.
    betaincinv   -- Inverse of betainc.
    betaincinva  -- Inverse (in first argument, a) of betainc
    betaincinvb  -- Inverse (in first argument, b) of betainc        
    psi(digamma) -- Logarithmic derivative of the gamma function.
    rgamma       -- One divided by the gamma function.
    polygamma    -- Nth derivative of psi function.

  Error Function and Fresnel Integrals
  
    erf        -- Error function.
    erfc       -- Complemented error function (1- erf(x))
    ervinv     -- Inverse of error function
    ervcinv    -- Inverse of erfc
    fresnel    -- Fresnel sine and cosine integrals.

  HyperGeometric Functions

    hyp2f1   -- Gauss hypergeometric function (2F1)
    hyp1f1   -- Confluent hypergeometric function (1F1)
    hyperu   -- Confluent hypergeometric function (U)
    hyp0f1   -- Confluent hypergeometric limit function (0F1)
    hyp2f0   -- Hypergeometric function (2F0) 
    hyp1f2   -- Hypergeometric function (1F2)
    hyp3f0   -- Hypergeometric function (3F0)

  Mathieu and Related Functions (and derivatives)

    mathieu_a       -- Characteristic values for even solution (ce_m)
    mathieu_b       -- Characteristic values for odd solution (se_m)
    mathieu_A       -- **sequence of expansion coefficients for even solution
    mathieu_B       -- **sequence of expansion coefficients for odd solution
       ** All the following return both function and first derivative **
    mathieu_cem     -- Even Mathieu function
    mathieu_sem     -- Odd Mathieu function
    mathieu_modcem1 -- Even modified Mathieu function of the first kind
    mathieu_modcem2 -- Even modified Mathieu function of the second kind
    mathieu_modsem1 -- Odd modified Mathieu function of the first kind
    mathieu_modsem2 -- Odd modified Mathieu function of the second kind

  Other Special Functions

    expn         -- Exponential integral.
    exp1         -- Exponential integral of order 1 (for complex argument)
    expi         -- Another exponential integral -- Ei(x)
    wofz         -- Fadeeva function.
    dawsn        -- Dawson's integral.
    shichi       -- Hyperbolic sine and cosine integrals.
    sici         -- Integral of the sinc and "cosinc" functions.
    spence       -- Dilogarithm integral.
    struve       -- Struve function --- Hv(x)
    modstruve    -- Modified struve function --- Lv(x)
    itstruve0    -- Integral of H0(t) from 0 to x
    it2struve0   -- Integral of H0(t)/t from x to Inf.
    itmodstruve0 -- Integral of L0(t) from 0 to x.
    kelvin       -- Kelvin functions (order 0) and derivatives.
    zeta         -- Riemann zeta function of two arguments.
    zetac        -- 1.0 - standard Riemann zeta function.
    besselpoly   -- Integral of a bessel function times x**lambda.

  Convenience Functions

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
    round    -- round the argument to the nearest integer. If argument
                 ends in 0.5 exactly, pick the nearest even integer.

  ** by the description indicates a function whis is not a universal
  function and does not follow broadcasting and automatic
  array-looping rules.
 """
from special import *
import specfun



