.. _continuous-random-variables:

====================================
Continuous Statistical Distributions
====================================

Overview
========

All distributions will have location (L) and Scale (S) parameters
along with any shape parameters needed, the names for the shape
parameters will vary. Standard form for the distributions will be
given where :math:`L=0.0` and :math:`S=1.0.` The nonstandard forms can be obtained for the various functions using
(note :math:`U` is a standard uniform random variate).


======================================  =============================================================================  =========================================================================================================================================
Function Name                           Standard Function                                                              Transformation
======================================  =============================================================================  =========================================================================================================================================
Cumulative Distribution Function (CDF)  :math:`F\left(x\right)`                                                        :math:`F\left(x;L,S\right)=F\left(\frac{\left(x-L\right)}{S}\right)`
Probability Density Function (PDF)      :math:`f\left(x\right)=F^{\prime}\left(x\right)`                               :math:`f\left(x;L,S\right)=\frac{1}{S}f\left(\frac{\left(x-L\right)}{S}\right)`
Percent Point Function (PPF)            :math:`G\left(q\right)=F^{-1}\left(q\right)`                                   :math:`G\left(q;L,S\right)=L+SG\left(q\right)`
Probability Sparsity Function (PSF)     :math:`g\left(q\right)=G^{\prime}\left(q\right)`                               :math:`g\left(q;L,S\right)=Sg\left(q\right)`
Hazard Function (HF)                    :math:`h_{a}\left(x\right)=\frac{f\left(x\right)}{1-F\left(x\right)}`          :math:`h_{a}\left(x;L,S\right)=\frac{1}{S}h_{a}\left(\frac{\left(x-L\right)}{S}\right)`
Cumulative Hazard Function (CHF)        :math:`H_{a}\left(x\right)=` :math:`\log\frac{1}{1-F\left(x\right)}`           :math:`H_{a}\left(x;L,S\right)=H_{a}\left(\frac{\left(x-L\right)}{S}\right)`
Survival Function (SF)                  :math:`S\left(x\right)=1-F\left(x\right)`                                      :math:`S\left(x;L,S\right)=S\left(\frac{\left(x-L\right)}{S}\right)`
Inverse Survival Function (ISF)         :math:`Z\left(\alpha\right)=S^{-1}\left(\alpha\right)=G\left(1-\alpha\right)`  :math:`Z\left(\alpha;L,S\right)=L+SZ\left(\alpha\right)`
Moment Generating Function (MGF)        :math:`M_{Y}\left(t\right)=E\left[e^{Yt}\right]`                               :math:`M_{X}\left(t\right)=e^{Lt}M_{Y}\left(St\right)`
Random Variates                         :math:`Y=G\left(U\right)`                                                      :math:`X=L+SY`
(Differential) Entropy                  :math:`h\left[Y\right]=-\int f\left(y\right)\log f\left(y\right)dy`            :math:`h\left[X\right]=h\left[Y\right]+\log S`
(Non-central) Moments                   :math:`\mu_{n}^{\prime}=E\left[Y^{n}\right]`                                   :math:`E\left[X^{n}\right]=L^{n}\sum_{k=0}^{N}\left(\begin{array}{c} n\\ k\end{array}\right)\left(\frac{S}{L}\right)^{k}\mu_{k}^{\prime}`
Central Moments                         :math:`\mu_{n}=E\left[\left(Y-\mu\right)^{n}\right]`                           :math:`E\left[\left(X-\mu_{X}\right)^{n}\right]=S^{n}\mu_{n}`
mean (mode, median), var                :math:`\mu,\,\mu_{2}`                                                          :math:`L+S\mu,\, S^{2}\mu_{2}`
skewness                                :math:`\gamma_{1}=\frac{\mu_{3}}{\left(\mu_{2}\right)^{3/2}}`                  :math:`\gamma_{1}`
kurtosis                                :math:`\gamma_{2}=\frac{\mu_{4}}{\left(\mu_{2}\right)^{2}}-3`                  :math:`\gamma_{2}`
======================================  =============================================================================  =========================================================================================================================================


Moments
-------

Non-central moments are defined using the PDF

.. math::

   \mu_{n}^{\prime}=\int_{-\infty}^{\infty}x^{n}f\left(x\right)dx.

Note, that these can always be computed using the PPF. Substitute :math:`x=G\left(q\right)` in the above equation and get

.. math::

   \mu_{n}^{\prime}=\int_{0}^{1}G^{n}\left(q\right)dq

which may be easier to compute numerically. Note that :math:`q=F\left(x\right)` so that :math:`dq=f\left(x\right)dx.` Central moments are computed similarly :math:`\mu=\mu_{1}^{\prime}`

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{n} & = & \int_{-\infty}^{\infty}\left(x-\mu\right)^{n}f\left(x\right)dx\\  & = & \int_{0}^{1}\left(G\left(q\right)-\mu\right)^{n}dq\\  & = & \sum_{k=0}^{n}\left(\begin{array}{c} n\\ k\end{array}\right)\left(-\mu\right)^{k}\mu_{n-k}^{\prime}\end{eqnarray*}

In particular

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{3} & = & \mu_{3}^{\prime}-3\mu\mu_{2}^{\prime}+2\mu^{3}\\  & = & \mu_{3}^{\prime}-3\mu\mu_{2}-\mu^{3}\\ \mu_{4} & = & \mu_{4}^{\prime}-4\mu\mu_{3}^{\prime}+6\mu^{2}\mu_{2}^{\prime}-3\mu^{4}\\  & = & \mu_{4}^{\prime}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}\end{eqnarray*}

Skewness is defined as

.. math::

     \gamma_{1}=\sqrt{\beta_{1}}=\frac{\mu_{3}}{\mu_{2}^{3/2}}

while (Fisher) kurtosis is

.. math::

     \gamma_{2}=\frac{\mu_{4}}{\mu_{2}^{2}}-3,

so that a normal distribution has a kurtosis of zero.


Median and mode
---------------

The median, :math:`m_{n}` is defined as the point at which half of the density is on one side
and half on the other. In other words, :math:`F\left(m_{n}\right)=\frac{1}{2}` so that

.. math::

     m_{n}=G\left(\frac{1}{2}\right).

In addition, the mode, :math:`m_{d}` , is defined as the value for which the probability density function
reaches it's peak

.. math::

     m_{d}=\arg\max_{x}f\left(x\right).


Fitting data
------------

To fit data to a distribution, maximizing the likelihood function is
common. Alternatively, some distributions have well-known minimum
variance unbiased estimators. These will be chosen by default, but the
likelihood function will always be available for minimizing.

If :math:`f\left(x;\boldsymbol{\theta}\right)` is the PDF of a random-variable where :math:`\boldsymbol{\theta}` is a vector of parameters ( *e.g.* :math:`L` and :math:`S` ), then for a collection of :math:`N` independent samples from this distribution, the joint distribution the
random vector :math:`\mathbf{x}` is

.. math::

     f\left(\mathbf{x};\boldsymbol{\theta}\right)=\prod_{i=1}^{N}f\left(x_{i};\boldsymbol{\theta}\right).

The maximum likelihood estimate of the parameters :math:`\boldsymbol{\theta}` are the parameters which maximize this function with :math:`\mathbf{x}` fixed and given by the data:

.. math::
   :nowrap:

    \begin{eqnarray*} \boldsymbol{\theta}_{es} & = & \arg\max_{\boldsymbol{\theta}}f\left(\mathbf{x};\boldsymbol{\theta}\right)\\  & = & \arg\min_{\boldsymbol{\theta}}l_{\mathbf{x}}\left(\boldsymbol{\theta}\right).\end{eqnarray*}

Where

.. math::
   :nowrap:

    \begin{eqnarray*} l_{\mathbf{x}}\left(\boldsymbol{\theta}\right) & = & -\sum_{i=1}^{N}\log f\left(x_{i};\boldsymbol{\theta}\right)\\  & = & -N\overline{\log f\left(x_{i};\boldsymbol{\theta}\right)}\end{eqnarray*}

Note that if :math:`\boldsymbol{\theta}` includes only shape parameters, the location and scale-parameters can
be fit by replacing :math:`x_{i}` with :math:`\left(x_{i}-L\right)/S` in the log-likelihood function adding :math:`N\log S` and minimizing, thus

.. math::
   :nowrap:

    \begin{eqnarray*} l_{\mathbf{x}}\left(L,S;\boldsymbol{\theta}\right) & = & N\log S-\sum_{i=1}^{N}\log f\left(\frac{x_{i}-L}{S};\boldsymbol{\theta}\right)\\  & = & N\log S+l_{\frac{\mathbf{x}-S}{L}}\left(\boldsymbol{\theta}\right)\end{eqnarray*}

If desired, sample estimates for :math:`L` and :math:`S` (not necessarily maximum likelihood estimates) can be obtained from
samples estimates of the mean and variance using

.. math::
   :nowrap:

    \begin{eqnarray*} \hat{S} & = & \sqrt{\frac{\hat{\mu}_{2}}{\mu_{2}}}\\ \hat{L} & = & \hat{\mu}-\hat{S}\mu\end{eqnarray*}

where :math:`\mu` and :math:`\mu_{2}` are assumed known as the mean and variance of the **untransformed** distribution (when :math:`L=0` and :math:`S=1` ) and

.. math::
   :nowrap:

    \begin{eqnarray*} \hat{\mu} & = & \frac{1}{N}\sum_{i=1}^{N}x_{i}=\bar{\mathbf{x}}\\ \hat{\mu}_{2} & = & \frac{1}{N-1}\sum_{i=1}^{N}\left(x_{i}-\hat{\mu}\right)^{2}=\frac{N}{N-1}\overline{\left(\mathbf{x}-\bar{\mathbf{x}}\right)^{2}}\end{eqnarray*}


Standard notation for mean
--------------------------

We will use

.. math::

    \overline{y\left(\mathbf{x}\right)}=\frac{1}{N}\sum_{i=1}^{N}y\left(x_{i}\right)

where :math:`N` should be clear from context as the number of samples :math:`x_{i}`

References
----------

-  Documentation for ranlib, rv2, cdflib

-  Eric Weisstein's world of mathematics http://mathworld.wolfram.com/,
   http://mathworld.wolfram.com/topics/StatisticalDistributions.html

-  Documentation to Regress+ by Michael McLaughlin item Engineering and
   Statistics Handbook (NIST),
   https://www.itl.nist.gov/div898/handbook/

-  Documentation for DATAPLOT from NIST,
   https://www.itl.nist.gov/div898/software/dataplot/distribu.htm

-  Norman Johnson, Samuel Kotz, and N. Balakrishnan Continuous
   Univariate Distributions, second edition, Volumes I and II, Wiley &
   Sons, 1994.


In the tutorials several special functions appear repeatedly and are listed here.

===============================================================  ======================================================================================  =============================================================================================================================
Symbol                                                           Description                                                                             Definition
===============================================================  ======================================================================================  =============================================================================================================================
:math:`\gamma\left(s, x\right)`                                  lower incomplete Gamma function                                                         :math:`\int_0^x t^{s-1} e^{-t} dt`
:math:`\Gamma\left(s, x\right)`                                  upper incomplete Gamma function                                                         :math:`\int_x^\infty t^{s-1} e^{-t} dt`
:math:`B\left(x;a,b\right)`                                      incomplete Beta function                                                                :math:`\int_{0}^{x} t^{a-1}\left(1-t\right)^{b-1} dt`
:math:`I\left(x;a,b\right)`                                      regularized incomplete Beta function                                                    :math:`\frac{\Gamma\left(a+b\right)}{\Gamma\left(a\right)\Gamma\left(b\right)} \int_{0}^{x} t^{a-1}\left(1-t\right)^{b-1} dt`
:math:`\phi\left(x\right)`                                       PDF for normal distribution                                                             :math:`\frac{1}{\sqrt{2\pi}}e^{-x^{2}/2}`
:math:`\Phi\left(x\right)`                                       CDF for normal distribution                                                             :math:`\int_{-\infty}^{x}\phi\left(t\right) dt = \frac{1}{2}+\frac{1}{2}\mathrm{erf}\left(\frac{x}{\sqrt{2}}\right)`
:math:`\psi\left(z\right)`                                       digamma function                                                                        :math:`\frac{d}{dz} \log\left(\Gamma\left(z\right)\right)`
:math:`\psi_{n}\left(z\right)`                                   polygamma function                                                                      :math:`\frac{d^{n+1}}{dz^{n+1}}\log\left(\Gamma\left(z\right)\right)`
:math:`I_{\nu}\left(y\right)`                                    modified Bessel function of the first kind
:math:`\mathrm{Ei}(\mathrm{z})`                                  exponential integral                                                                    :math:`-\int_{-x}^\infty \frac{e^{-t}}{t} dt`
:math:`\zeta\left(n\right)`                                      Riemann zeta function                                                                   :math:`\sum_{k=1}^{\infty} \frac{1}{k^{n}}`
:math:`\zeta\left(n,z\right)`                                    Hurwitz zeta function                                                                   :math:`\sum_{k=0}^{\infty} \frac{1}{\left(k+z\right)^{n}}`
:math:`\,{}_{p}F_{q}(a_{1},\ldots,a_{p};b_{1},\ldots,b_{q};z)`   Hypergeometric function                                                                 :math:`\sum_{n=0}^{\infty} {\frac{(a_{1})_{n}\cdots(a_{p})_{n}}{(b_{1})_{n}\cdots(b_{q})_{n}}} \,{\frac{z^{n}}{n!}}`
===============================================================  ======================================================================================  =============================================================================================================================


Continuous Distributions in `scipy.stats`
-----------------------------------------
.. toctree::
   :maxdepth: 1

   continuous_alpha
   continuous_anglit
   continuous_arcsine
   continuous_beta
   continuous_betaprime
   continuous_bradford
   continuous_burr
   continuous_burr12
   continuous_cauchy
   continuous_skewcauchy
   continuous_chi
   continuous_chi2
   continuous_cosine
   continuous_dgamma
   continuous_dweibull
   continuous_erlang
   continuous_expon
   continuous_exponweib
   continuous_exponpow
   continuous_fatiguelife
   continuous_fisk
   continuous_foldcauchy
   continuous_foldnorm
   continuous_f
   continuous_gamma
   continuous_genlogistic
   continuous_genpareto
   continuous_genexpon
   continuous_genextreme
   continuous_gengamma
   continuous_genhalflogistic
   continuous_genhyperbolic
   continuous_geninvgauss
   continuous_gennorm
   continuous_gibrat
   continuous_gompertz
   continuous_gumbel_r
   continuous_gumbel_l
   continuous_halfcauchy
   continuous_halfnorm
   continuous_halflogistic
   continuous_hypsecant
   continuous_gausshyper
   continuous_invgamma
   continuous_invgauss
   continuous_invweibull
   continuous_johnsonsb
   continuous_johnsonsu
   continuous_ksone
   continuous_kstwo
   continuous_kstwobign
   continuous_laplace
   continuous_laplace_asymmetric
   continuous_levy_l
   continuous_levy
   continuous_logistic
   continuous_loglaplace
   continuous_loggamma
   continuous_lognorm
   continuous_loguniform
   continuous_maxwell
   continuous_mielke
   continuous_nakagami
   continuous_ncx2
   continuous_ncf
   continuous_nct
   continuous_norm
   continuous_norminvgauss
   continuous_pareto
   continuous_lomax
   continuous_powerlognorm
   continuous_powernorm
   continuous_powerlaw
   continuous_rdist
   continuous_rayleigh
   continuous_rice
   continuous_recipinvgauss
   continuous_semicircular
   continuous_studentized_range
   continuous_t
   continuous_trapezoid
   continuous_triang
   continuous_truncexpon
   continuous_truncnorm
   continuous_truncweibull_min
   continuous_tukeylambda
   continuous_uniform
   continuous_vonmises
   continuous_wald
   continuous_weibull_max
   continuous_weibull_min
   continuous_wrapcauchy
