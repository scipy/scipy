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

  Statistical Functions

    binomcdf       -- Sum of terms 0 through k of of the binomial pdf.
    binomcdfc      -- Sum of terms k+1 through n of the binomial pdf.
    binomcdfinv    -- Inverse of binomcdf.
    betacdf        -- Integral from 0 to x of beta pdf.
    betaq          -- Quantiles of beta distribution
    fcdf           -- Integral from 0 to x of F pdf.
    fcdfc          -- Integral from x to infinity under F pdf.
    fp             -- Inverse of fcdfc
    gammacdf       -- Integral from 0 to x of gamma pdf.
    gammacdfc      -- Integral from x to infinity under gamma pdf.
    gammaq         -- Quantiles of gamma distribution
    negbinomcdf    -- Sum of terms 0 through k of the negative binomial pdf.
    negbinomcdfc   -- Sum of terms k+1 to infinity under negative binomial pdf.
    negbinomcdfinv -- Inverse of negbinomcdf
    poissoncdf     -- Sum of terms 0 through k of the Poisson pdf.
    poissoncdfc    -- Sum of terms k+1 to infinity of the Poisson pdf.
    poissoncdfinv  -- Inverse of poissoncdf
    studentcdf     -- Integral from -infinity to t of the Student-t pdf.
    studentq       -- Inverse of studentcdf (quantiles)
    chi2cdf        -- Integral from 0 to x of the Chi-square pdf.
    chi2cdfc       -- Integral from x to infnity of Chi-square pdf.
    chi2p          -- Inverse of chi2cdfc.
    normalcdf      -- Integral from -infinity to x of Gaussian pdf
    normalq        -- Inverse of normalcdf (quantiles)
    erf            -- Error function.
    erfc           -- Complemented error function (1- erf(x)) 
    smirnovcdfc    -- Exact Smirnov statistic for one-sided test.
    smirnovp       -- Inverse of smirnov.
    kolmogorovcdfc -- Kolmogorov's limiting distribution of a two-sided test.
    kolmogorovp    -- Inverse of kolmogorov.
      
  Gamma and Related Functions

    gamma        -- Gamma function.
    gammaln      -- Log of the absolute value of the gamma function.
    gammainc     -- Incomplete gamma integral.
    gammaincc    -- Complemented incomplete gamma integral.
    gammainccinv -- Inverse of gammaincc.
    beta         -- Beta function.
    betaln       -- Log of the absolute value of the beta function.
    betainc      -- Incomplete beta integral.
    betaincinv   -- Inverse of incbet.
    psi          -- Logarithmic derivative of the gamma function.
    rgamma       -- One divided by the gamma function.

  HyperGeometric Functions

    hyp2f1   -- Gauss hypergeometric function (2F1)
    hyp1f1   -- Confluent hypergeometric function (1F1)
    hyp2f0   -- Hypergeometric function (2F0) 
    hyp1f2   -- Hypergeometric function (1F2)
    hyp3f0   -- Hypergeometric function (3F0)

  Other Special Functions

    expn     -- Exponential integral.
    wofz     -- Fadeeva function.
    fresnl   -- Fresnel sine and cosine integrals.
    dawsn    -- Dawson's integral.
    shichi   -- Hyperbolic sine and cosine integrals.
    sici     -- Integral of the sinc and "cosinc" functions.
    spence   -- Dilogarithm integral.
    struve   -- Struve function.
    zeta     -- Riemann zeta function of two arguments.
    zetac    -- Riemann zeta function.
    besselpoly -- Integral of a bessel function times x**lambda.

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

 """
__all__ = []
import scipy
scipy.names2all(__all__, ['special'], globals())
del scipy


