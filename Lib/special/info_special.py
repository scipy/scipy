"""
Special Functions
=================

  Airy Functions

    airy     -- Airy functions and their derivatives.
    airye    -- Exponentially scaled Airy functions
    ai_zeros -- **Zeros of Airy functions Ai(x) and Ai'(x)
    bi_zeros -- **Zeros of Airy functions Bi(x) and Bi'(x)

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

    lmbda       -- **Sequence of lambda functions with arbitrary order v.

  Zeros of Bessel Functions

    jnjnp_zeros -- **Zeros of integer-order Bessel functions and derivatives
                       sorted in order.
    jnyn_zeros  -- **Zeros of integer-order Bessel functions and derivatives
                       as separate arrays.
    jn_zeros    -- **Zeros of Jn(x)
    jnp_zeros   -- **Zeros of Jn'(x)
    yn_zeros    -- **Zeros of Yn(x)
    ynp_zeros   -- **Zeros of Yn'(x)
    y0_zeros    -- **Complex zeros: Y0(z0)=0 and values of Y0'(z0)
    y1_zeros    -- **Complex zeros: Y1(z1)=0 and values of Y1'(z1)
    y1p_zeros   -- **Complex zeros of Y1'(z1')=0 and values of Y1(z1')

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

    itj0y0     -- Basic integrals of j0 and y0 from 0 to x.
    it2j0y0    -- Integrals of (1-j0(t))/t from 0 to x and
                    y0(t)/t from x to inf.
    iti0k0     -- Basic integrals of i0 and k0 from 0 to x.
    it2i0k0    -- Integrals of (i0(t)-1)/t from 0 to x and
                    k0(t)/t from x to inf.
    besselpoly -- Integral of a bessel function: Jv(2*a*x) * x**lambda
                    from x=0 to 1.

  Derivatives of Bessel Functions.

    jvp     -- Nth derivative of Jv(v,z)
    yvp     -- Nth derivative of Yv(v,z)
    kvp     -- Nth derivative of Kv(v,z)
    ivp     -- Nth derivative of Iv(v,z)
    h1vp    -- Nth derivative of H1v(v,z)
    h2vp    -- Nth derivative of H2v(v,z)    

  Spherical Bessel Functions

    sph_jn   -- **Sequence of spherical Bessel functions, jn(z)
    sph_yn   -- **Sequence of spherical Bessel functions, yn(z)
    sph_jnyn -- **Sequence of spherical Bessel functions, jn(z) and yn(z)
    sph_in   -- **Sequence of spherical Bessel functions, in(z)
    sph_kn   -- **Sequence of spherical Bessel functions, kn(z)
    sph_inkn -- **Sequence of spherical Bessel functions, in(z) and kn(z)

  Ricatti-Bessel Functions

    ricatti_jn -- **Sequence of Ricatti-Bessel functions of first kind.
    ricatti_yn -- **Sequence of Ricatti-Bessel functions of second kind.

  Struve Functions

    struve       -- Struve function --- Hv(x)
    modstruve    -- Modified struve function --- Lv(x)
    itstruve0    -- Integral of H0(t) from 0 to x
    it2struve0   -- Integral of H0(t)/t from x to Inf.
    itmodstruve0 -- Integral of L0(t) from 0 to x.
    
        
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
  
    erf           -- Error function.
    erfc          -- Complemented error function (1- erf(x))
    erfinv        -- Inverse of error function
    erfcinv       -- Inverse of erfc
    erf_zeros     -- **Complex zeros of erf(z)
    fresnel       -- Fresnel sine and cosine integrals.
    fresnel_zeros -- Complex zeros of both Fresnel integrals
    fresnelc_zeros -- **Complex zeros of fresnel cosine integrals
    fresnels_zeros -- **Complex zeros of fresnel sine integrals
    modfresnelp   -- Modified Fresnel integrals F_+(x) and K_+(x)
    modfresnelm   -- Modified Fresnel integrals F_-(x) and K_-(x)

  Legendre Functions

    lpn      -- **Legendre Functions (polynomials) of the first kind
    lqn      -- **Legendre Functions of the second kind.
    lpmn     -- **Associated Legendre Function of the first kind.
    lqmn     -- **Associated Legendre Function of the second kind.
    lpmv     -- Associated Legendre Function of arbitrary non-negative
                   degree v.
    sph_harm -- Spherical Harmonics (complex-valued) Y^m_n(theta,phi)

  Orthogonal polynomials  --- 15 types
   ** These functions all return a polynomial class which can then be
      evaluated:  vals = chebyt(n)(x)
      This class also has an attribute 'weights' which
      return the roots, weights, and total weights for the appropriate 
      form of Gaussian quadrature.  These are returned in an n x 3 array with roots
      in the first column, weights in the second column, and total weights in the final
      column

    legendre    -- **Legendre polynomial P_n(x) (lpn -- for function).
    chebyt      -- **Chebyshev polynomial T_n(x)
    chebyu      -- **Chebyshev polynomial U_n(x)
    chebyc      -- **Chebyshev polynomial C_n(x)
    chebys      -- **Chebyshev polynomial S_n(x)
    jacobi      -- **Jacobi polynomial P^(alpha,beta)_n(x)
    laguerre    -- **Laguerre polynomial, L_n(x)
    genlaguerre -- **Generalized (Associated) Laguerre polynomial, L^alpha_n(x)
    hermite     -- **Hermite polynomial H_n(x)
    hermitenorm -- **Normalized Hermite polynomial, He_n(x)
    gegenbauer  -- **Gegenbauer (Ultraspherical) polynomials, C^(alpha)_n(x)
    sh_legendre -- **shifted Legendre polynomial, P*_n(x)
    sh_chebyt   -- **shifted Chebyshev polynomial, T*_n(x)
    sh_chebyu   -- **shifted Chebyshev polynomial, U*_n(x)
    sh_jacobi   -- **shifted Jacobi polynomial, J*_n(x) = G^(p,q)_n(x)
      
  HyperGeometric Functions

    hyp2f1   -- Gauss hypergeometric function (2F1)
    hyp1f1   -- Confluent hypergeometric function (1F1)
    hyperu   -- Confluent hypergeometric function (U)
    hyp0f1   -- Confluent hypergeometric limit function (0F1)
    hyp2f0   -- Hypergeometric function (2F0) 
    hyp1f2   -- Hypergeometric function (1F2)
    hyp3f0   -- Hypergeometric function (3F0)

  Parabolic Cylinder Functions

    pbdv     -- Parabolic cylinder function Dv(x) and derivative.
    pbvv     -- Parabolic cylinder function Vv(x) and derivative. 
    pbwa     -- Parabolic cylinder function W(a,x) and derivative.
    pbdv_seq -- **Sequence of parabolic cylinder functions Dv(x)
    pbvv_seq -- **Sequence of parabolic cylinder functions Vv(x)
    pbdn_seq -- **Sequence of parabolic cylinder functions Dn(z), complex z

  mathieu and Related Functions (and derivatives)

    mathieu_a       -- Characteristic values for even solution (ce_m)
    mathieu_b       -- Characteristic values for odd solution (se_m)
    mathieu_even_coef -- **sequence of expansion coefficients for even solution
    mathieu_odd_coef  -- **sequence of expansion coefficients for odd solution
       ** All the following return both function and first derivative **
    mathieu_cem     -- Even mathieu function
    mathieu_sem     -- Odd mathieu function
    mathieu_modcem1 -- Even modified mathieu function of the first kind
    mathieu_modcem2 -- Even modified mathieu function of the second kind
    mathieu_modsem1 -- Odd modified mathieu function of the first kind
    mathieu_modsem2 -- Odd modified mathieu function of the second kind

  Spheroidal Wave Functions

    pro_ang1   -- Prolate spheroidal angular function of the first kind
    pro_rad1   -- Prolate spheroidal radial function of the first kind
    pro_rad2   -- Prolate spheroidal radial function of the second kind
    obl_ang1   -- Oblate spheroidal angluar function of the first kind
    obl_rad1   -- Oblate spheroidal radial function of the first kind
    obl_rad2   -- Oblate spheroidal radial function of the second kind
    pro_cv     -- Compute characteristic value for prolate functions
    obl_cv     -- Compute characteristic value for oblate functions
    pro_cv_seq -- Compute sequence of prolate characteristic values
    obl_cv_seq -- Compute sequence of oblate characteristic values    
     ** The following functions require pre-computed characteristic values **
    pro_ang1_cv -- Prolate spheroidal angular function of the first kind
    pro_rad1_cv -- Prolate spheroidal radial function of the first kind
    pro_rad2_cv -- Prolate spheroidal radial function of the second kind
    obl_ang1_cv -- Oblate spheroidal angluar function of the first kind
    obl_rad1_cv -- Oblate spheroidal radial function of the first kind
    obl_rad2_cv -- Oblate spheroidal radial function of the second kind

  Kelvin Functions
  
    kelvin       -- All Kelvin functions (order 0) and derivatives.
    kelvin_zeros -- **Zeros of All Kelvin functions (order 0) and derivatives
    ber          -- Kelvin function ber x
    bei          -- Kelvin function bei x
    berp         -- Derivative of Kelvin function ber x
    beip         -- Derivative of Kelvin function bei x
    ker          -- Kelvin function ker x
    kei          -- Kelvin function kei x
    kerp         -- Derivative of Kelvin function ker x
    keip         -- Derivative of Kelvin function kei x
    ber_zeros    -- **Zeros of Kelvin function bei x
    bei_zeros    -- **Zeros of Kelvin function ber x
    berp_zeros   -- **Zeros of derivative of Kelvin function ber x
    beip_zeros   -- **Zeros of derivative of Kelvin function bei x
    ker_zeros    -- **Zeros of Kelvin function kei x
    kei_zeros    -- **Zeros of Kelvin function ker x
    kerp_zeros   -- **Zeros of derivative of Kelvin function ker x
    keip_zeros   -- **Zeros of derivative of Kelvin function kei x
  
  Other Special Functions

    expn         -- Exponential integral.
    exp1         -- Exponential integral of order 1 (for complex argument)
    expi         -- Another exponential integral -- Ei(x)
    wofz         -- Fadeeva function.
    dawsn        -- Dawson's integral.
    shichi       -- Hyperbolic sine and cosine integrals.
    sici         -- Integral of the sinc and "cosinc" functions.
    spence       -- Dilogarithm integral.
    zeta         -- Riemann zeta function of two arguments.
    zetac        -- 1.0 - standard Riemann zeta function.

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

  ** in the description indicates a function which is not a universal
  function and does not follow broadcasting and automatic
  array-looping rules.

   Error handling:

      Errors are handled by returning nans, or other appropriate values.
      Some of the special function routines will print an error message
      when an error occurs.  By default this printing
      is disabled.  To enable such messages use errprint(1)
      To disable such messages use errprint(0).

      Example:
      >>> print scipy.special.bdtr(-1,10,0.3)
      >>> scipy.special.errprint(1)
      >>> print scipy.special.bdtr(-1,10,0.3)
 
"""

postpone_import = 1
global_symbols = ['isinf','isfinite','isnan']
