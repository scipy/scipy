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


======================================  ==============================================================================================================================  =========================================================================================================================================
Function Name                           Standard Function                                                                                                               Transformation
======================================  ==============================================================================================================================  =========================================================================================================================================
Cumulative Distribution Function (CDF)  :math:`F\left(x\right)`                                                                                                         :math:`F\left(x;L,S\right)=F\left(\frac{\left(x-L\right)}{S}\right)`
Probability Density Function (PDF)      :math:`f\left(x\right)=F^{\prime}\left(x\right)`                                                                                :math:`f\left(x;L,S\right)=\frac{1}{S}f\left(\frac{\left(x-L\right)}{S}\right)`
Percent Point Function (PPF)            :math:`G\left(q\right)=F^{-1}\left(q\right)`                                                                                    :math:`G\left(q;L,S\right)=L+SG\left(q\right)`
Probability Sparsity Function (PSF)     :math:`g\left(q\right)=G^{\prime}\left(q\right)`                                                                                :math:`g\left(q;L,S\right)=Sg\left(q\right)`
Hazard Function (HF)                    :math:`h_{a}\left(x\right)=\frac{f\left(x\right)}{1-F\left(x\right)}`                                                           :math:`h_{a}\left(x;L,S\right)=\frac{1}{S}h_{a}\left(\frac{\left(x-L\right)}{S}\right)`
Cumulative Hazard Functon (CHF)         :math:`H_{a}\left(x\right)=` :math:`\log\frac{1}{1-F\left(x\right)}`                                                            :math:`H_{a}\left(x;L,S\right)=H_{a}\left(\frac{\left(x-L\right)}{S}\right)`
Survival Function (SF)                  :math:`S\left(x\right)=1-F\left(x\right)`                                                                                       :math:`S\left(x;L,S\right)=S\left(\frac{\left(x-L\right)}{S}\right)`
Inverse Survival Function (ISF)         :math:`Z\left(\alpha\right)=S^{-1}\left(\alpha\right)=G\left(1-\alpha\right)`                                                   :math:`Z\left(\alpha;L,S\right)=L+SZ\left(\alpha\right)`
Moment Generating Function (MGF)        :math:`M_{Y}\left(t\right)=E\left[e^{Yt}\right]`                                                                                :math:`M_{X}\left(t\right)=e^{Lt}M_{Y}\left(St\right)`
Random Variates                         :math:`Y=G\left(U\right)`                                                                                                       :math:`X=L+SY`
(Differential) Entropy                  :math:`h\left[Y\right]=-\int f\left(y\right)\log f\left(y\right)dy`                                                             :math:`h\left[X\right]=h\left[Y\right]+\log S`
(Non-central) Moments                   :math:`\mu_{n}^{\prime}=E\left[Y^{n}\right]`                                                                                    :math:`E\left[X^{n}\right]=L^{n}\sum_{k=0}^{N}\left(\begin{array}{c} n\\ k\end{array}\right)\left(\frac{S}{L}\right)^{k}\mu_{k}^{\prime}`
Central Moments                         :math:`\mu_{n}=E\left[\left(Y-\mu\right)^{n}\right]`                                                                            :math:`E\left[\left(X-\mu_{X}\right)^{n}\right]=S^{n}\mu_{n}`
mean (mode, median), var                :math:`\mu,\,\mu_{2}`                                                                                                           :math:`L+S\mu,\, S^{2}\mu_{2}`
skewness, kurtosis                      :math:`\gamma_{1}=\frac{\mu_{3}}{\left(\mu_{2}\right)^{3/2}},\,` :math:`\gamma_{2}=\frac{\mu_{4}}{\left(\mu_{2}\right)^{2}}-3`  :math:`\gamma_{1},\,\gamma_{2}`
======================================  ==============================================================================================================================  =========================================================================================================================================


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

-  Eric Weisstein~s world of mathematics http://mathworld.wolfram.com/,
   http://mathworld.wolfram.com/topics/StatisticalDistributions.html

-  Documentation to Regress+ by Michael McLaughlin item Engineering and
   Statistics Handbook (NIST),
   http://www.itl.nist.gov/div898/handbook/index.htm

-  Documentation for DATAPLOT from NIST,
   http://www.itl.nist.gov/div898/software/dataplot/distribu.htm

-  Norman Johnson, Samuel Kotz, and N. Balakrishnan Continuous
   Univariate Distributions, second edition, Volumes I and II, Wiley &
   Sons, 1994.


Alpha
=====

One shape parameters :math:`\alpha>0` (parameter :math:`\beta` in DATAPLOT
is a scale-parameter). Standard form is :math:`x>0:`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha\right) & = & \frac{1}{x^{2}\Phi\left(\alpha\right)\sqrt{2\pi}}\exp\left(-\frac{1}{2}\left(\alpha-\frac{1}{x}\right)^{2}\right)\\ F\left(x;\alpha\right) & = & \frac{\Phi\left(\alpha-\frac{1}{x}\right)}{\Phi\left(\alpha\right)}\\ G\left(q;\alpha\right) & = & \left[\alpha-\Phi^{-1}\left(q\Phi\left(\alpha\right)\right)\right]^{-1}\end{eqnarray*}

.. math::

     M\left(t\right)=\frac{1}{\Phi\left(a\right)\sqrt{2\pi}}\int_{0}^{\infty}\frac{e^{xt}}{x^{2}}\exp\left(-\frac{1}{2}\left(\alpha-\frac{1}{x}\right)^{2}\right)dx

No moments?

.. math::

     l_{\mathbf{x}}\left(\alpha\right)=N\log\left[\Phi\left(\alpha\right)\sqrt{2\pi}\right]+2N\overline{\log\mathbf{x}}+\frac{N}{2}\alpha^{2}-\alpha\overline{\mathbf{x}^{-1}}+\frac{1}{2}\overline{\mathbf{x}^{-2}}

Implementation: `scipy.stats.alpha`


Anglit
======

Defined over :math:`x\in\left[-\frac{\pi}{4},\frac{\pi}{4}\right]`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \sin\left(2x+\frac{\pi}{2}\right)=\cos\left(2x\right)\\ F\left(x\right) & = & \sin^{2}\left(x+\frac{\pi}{4}\right)\\ G\left(q\right) & = & \arcsin\left(\sqrt{q}\right)-\frac{\pi}{4}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 0\\ \mu_{2} & = & \frac{\pi^{2}}{16}-\frac{1}{2}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & -2\frac{\pi^{4}-96}{\left(\pi^{2}-8\right)^{2}}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & 1-\log2\\  & \approx & 0.30685281944005469058\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & \int_{-\frac{\pi}{4}}^{\frac{\pi}{4}}\cos\left(2x\right)e^{xt}dx\\  & = & \frac{4\cosh\left(\frac{\pi t}{4}\right)}{t^{2}+4}\end{eqnarray*}

.. math::

     l_{\mathbf{x}}\left(\cdot\right)=-N\overline{\log\left[\cos\left(2\mathbf{x}\right)\right]}

Implementation: `scipy.stats.anglit`


Arcsine
=======

Defined over :math:`x\in\left(0,1\right)` . To get the JKB definition put :math:`x=\frac{u+1}{2}.` i.e. :math:`L=-1` and :math:`S=2.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{\pi\sqrt{x\left(1-x\right)}}\\ F\left(x\right) & = & \frac{2}{\pi}\arcsin\left(\sqrt{x}\right)\\ G\left(q\right) & = & \sin^{2}\left(\frac{\pi}{2}q\right)\end{eqnarray*}

.. math::

     M\left(t\right)=E^{t/2}I_{0}\left(\frac{t}{2}\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{n}^{\prime} & = & \frac{1}{\pi}\int_{0}^{1}dx\, x^{n-1/2}\left(1-x\right)^{-1/2}\\  & = & \frac{1}{\pi}B\left(\frac{1}{2},n+\frac{1}{2}\right)=\frac{\left(2n-1\right)!!}{2^{n}n!}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{1}{2}\\ \mu_{2} & = & \frac{1}{8}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & -\frac{3}{2}\end{eqnarray*}

.. math::

     h\left[X\right]\approx-0.24156447527049044468

.. math::

     l_{\mathbf{x}}\left(\cdot\right)=N\log\pi+\frac{N}{2}\overline{\log\mathbf{x}}+\frac{N}{2}\overline{\log\left(1-\mathbf{x}\right)}

Implementation: `scipy.stats.arcsine`


Beta
====

Two shape parameters

.. math::

     a,b>0

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,b\right) & = & \frac{\Gamma\left(a+b\right)}{\Gamma\left(a\right)\Gamma\left(b\right)}x^{a-1}\left(1-x\right)^{b-1}I_{\left(0,1\right)}\left(x\right)\\ F\left(x;a,b\right) & = & \int_{0}^{x}f\left(y;a,b\right)dy=I\left(x,a,b\right)\\ G\left(\alpha;a,b\right) & = & I^{-1}\left(\alpha;a,b\right)\\ M\left(t\right) & = & \frac{\Gamma\left(a\right)\Gamma\left(b\right)}{\Gamma\left(a+b\right)}\,_{1}F_{1}\left(a;a+b;t\right)\\ \mu & = & \frac{a}{a+b}\\ \mu_{2} & = & \frac{ab\left(a+b+1\right)}{\left(a+b\right)^{2}}\\ \gamma_{1} & = & 2\frac{b-a}{a+b+2}\sqrt{\frac{a+b+1}{ab}}\\ \gamma_{2} & = & \frac{6\left(a^{3}+a^{2}\left(1-2b\right)+b^{2}\left(b+1\right)-2ab\left(b+2\right)\right)}{ab\left(a+b+2\right)\left(a+b+3\right)}\\ m_{d} & = & \frac{\left(a-1\right)}{\left(a+b-2\right)}\, a+b\neq2\end{eqnarray*}

:math:`f\left(x;a,1\right)` is also called the Power-function distribution.

.. math::

     l_{\mathbf{x}}\left(a,b\right)=-N\log\Gamma\left(a+b\right)+N\log\Gamma\left(a\right)+N\log\Gamma\left(b\right)-N\left(a-1\right)\overline{\log\mathbf{x}}-N\left(b-1\right)\overline{\log\left(1-\mathbf{x}\right)}

All of the :math:`x_{i}\in\left[0,1\right]`

Implementation: `scipy.stats.beta`


Beta Prime
==========

Defined over :math:`0<x<\infty.` :math:`\alpha,\beta>0.` (Note the CDF evaluation uses Eq. 3.194.1 on pg. 313 of Gradshteyn &
Ryzhik (sixth edition).

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha,\beta\right) & = & \frac{\Gamma\left(\alpha+\beta\right)}{\Gamma\left(\alpha\right)\Gamma\left(\beta\right)}x^{\alpha-1}\left(1+x\right)^{-\alpha-\beta}\\ F\left(x;\alpha,\beta\right) & = & \frac{\Gamma\left(\alpha+\beta\right)}{\alpha\Gamma\left(\alpha\right)\Gamma\left(\beta\right)}x^{\alpha}\,_{2}F_{1}\left(\alpha+\beta,\alpha;1+\alpha;-x\right)\\ G\left(q;\alpha,\beta\right) & = & F^{-1}\left(x;\alpha,\beta\right)\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\left\{ \begin{array}{ccc} \frac{\Gamma\left(n+\alpha\right)\Gamma\left(\beta-n\right)}{\Gamma\left(\alpha\right)\Gamma\left(\beta\right)}=\frac{\left(\alpha\right)_{n}}{\left(\beta-n\right)_{n}} &  & \beta>n\\ \infty &  & \mathrm{otherwise}\end{array}\right.

Therefore,

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\alpha}{\beta-1}\quad\beta>1\\ \mu_{2} & = & \frac{\alpha\left(\alpha+1\right)}{\left(\beta-2\right)\left(\beta-1\right)}-\frac{\alpha^{2}}{\left(\beta-1\right)^{2}}\quad\beta>2\\ \gamma_{1} & = & \frac{\frac{\alpha\left(\alpha+1\right)\left(\alpha+2\right)}{\left(\beta-3\right)\left(\beta-2\right)\left(\beta-1\right)}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\quad\beta>3\\ \gamma_{2} & = & \frac{\mu_{4}}{\mu_{2}^{2}}-3\\ \mu_{4} & = & \frac{\alpha\left(\alpha+1\right)\left(\alpha+2\right)\left(\alpha+3\right)}{\left(\beta-4\right)\left(\beta-3\right)\left(\beta-2\right)\left(\beta-1\right)}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}\quad\beta>4\end{eqnarray*}

Implementation: `scipy.stats.betaprime`


Bradford
========

.. math::
   :nowrap:

    \begin{eqnarray*} c & > & 0\\ k & = & \log\left(1+c\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{c}{k\left(1+cx\right)}I_{\left(0,1\right)}\left(x\right)\\ F\left(x;c\right) & = & \frac{\log\left(1+cx\right)}{k}\\ G\left(\alpha\; c\right) & = & \frac{\left(1+c\right)^{\alpha}-1}{c}\\ M\left(t\right) & = & \frac{1}{k}e^{-t/c}\left[\mathrm{Ei}\left(t+\frac{t}{c}\right)-\mathrm{Ei}\left(\frac{t}{c}\right)\right]\\ \mu & = & \frac{c-k}{ck}\\ \mu_{2} & = & \frac{\left(c+2\right)k-2c}{2ck^{2}}\\ \gamma_{1} & = & \frac{\sqrt{2}\left(12c^{2}-9kc\left(c+2\right)+2k^{2}\left(c\left(c+3\right)+3\right)\right)}{\sqrt{c\left(c\left(k-2\right)+2k\right)}\left(3c\left(k-2\right)+6k\right)}\\ \gamma_{2} & = & \frac{c^{3}\left(k-3\right)\left(k\left(3k-16\right)+24\right)+12kc^{2}\left(k-4\right)\left(k-3\right)+6ck^{2}\left(3k-14\right)+12k^{3}}{3c\left(c\left(k-2\right)+2k\right)^{2}}\\ m_{d} & = & 0\\ m_{n} & = & \sqrt{1+c}-1\end{eqnarray*}

where :math:`\mathrm{Ei}\left(\mathrm{z}\right)` is the exponential integral function. Also

.. math::

     h\left[X\right]=\frac{1}{2}\log\left(1+c\right)-\log\left(\frac{c}{\log\left(1+c\right)}\right)

Implementation: `scipy.stats.bradford`


Burr
====

.. math::
   :nowrap:

    \begin{eqnarray*} c & > & 0\\ d & > & 0\\ k & = & \Gamma\left(d\right)\Gamma\left(1-\frac{2}{c}\right)\Gamma\left(\frac{2}{c}+d\right)-\Gamma^{2}\left(1-\frac{1}{c}\right)\Gamma^{2}\left(\frac{1}{c}+d\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c,d\right) & = & \frac{cd}{x^{c+1}\left(1+x^{-c}\right)^{d+1}}I_{\left(0,\infty\right)}\left(x\right)\\ F\left(x;c,d\right) & = & \left(1+x^{-c}\right)^{-d}\\ G\left(\alpha;c,d\right) & = & \left(\alpha^{-1/d}-1\right)^{-1/c}\\ \mu & = & \frac{\Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+d\right)}{\Gamma\left(d\right)}\\ \mu_{2} & = & \frac{k}{\Gamma^{2}\left(d\right)}\\ \gamma_{1} & = & \frac{1}{\sqrt{k^{3}}}\left[2\Gamma^{3}\left(1-\frac{1}{c}\right)\Gamma^{3}\left(\frac{1}{c}+d\right)+\Gamma^{2}\left(d\right)\Gamma\left(1-\frac{3}{c}\right)\Gamma\left(\frac{3}{c}+d\right)\right.\\  &  & \left.-3\Gamma\left(d\right)\Gamma\left(1-\frac{2}{c}\right)\Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+d\right)\Gamma\left(\frac{2}{c}+d\right)\right]\\ \gamma_{2} & = & -3+\frac{1}{k^{2}}\left[6\Gamma\left(d\right)\Gamma\left(1-\frac{2}{c}\right)\Gamma^{2}\left(1-\frac{1}{c}\right)\Gamma^{2}\left(\frac{1}{c}+d\right)\Gamma\left(\frac{2}{c}+d\right)\right.\\  &  & -3\Gamma^{4}\left(1-\frac{1}{c}\right)\Gamma^{4}\left(\frac{1}{c}+d\right)+\Gamma^{3}\left(d\right)\Gamma\left(1-\frac{4}{c}\right)\Gamma\left(\frac{4}{c}+d\right)\\  &  & \left.-4\Gamma^{2}\left(d\right)\Gamma\left(1-\frac{3}{c}\right)\Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+d\right)\Gamma\left(\frac{3}{c}+d\right)\right]\\ m_{d} & = & \left(\frac{cd-1}{c+1}\right)^{1/c}\,\mathrm{if }cd>1\,\mathrm{otherwise }0\\ m_{n} & = & \left(2^{1/d}-1\right)^{-1/c}\end{eqnarray*}

Implementation: `scipy.stats.burr`


Cauchy
======

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{\pi\left(1+x^{2}\right)}\\ F\left(x\right) & = & \frac{1}{2}+\frac{1}{\pi}\tan^{-1}x\\ G\left(\alpha\right) & = & \tan\left(\pi\alpha-\frac{\pi}{2}\right)\\ m_{d} & = & 0\\ m_{n} & = & 0\end{eqnarray*}

No finite moments. This is the t distribution with one degree of
freedom.

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(4\pi\right)\\  & \approx & 2.5310242469692907930.\end{eqnarray*}

Implementation: `scipy.stats.cauchy`


Chi
===

Generated by taking the (positive) square-root of chi-squared
variates.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\nu\right) & = & \frac{x^{\nu-1}e^{-x^{2}/2}}{2^{\nu/2-1}\Gamma\left(\frac{\nu}{2}\right)}I_{\left(0,\infty\right)}\left(x\right)\\ F\left(x;\nu\right) & = & \Gamma\left(\frac{\nu}{2},\frac{x^{2}}{2}\right)\\ G\left(\alpha;\nu\right) & = & \sqrt{2\Gamma^{-1}\left(\frac{\nu}{2},\alpha\right)}\end{eqnarray*}

.. math::

     M\left(t\right)=\Gamma\left(\frac{v}{2}\right)\,_{1}F_{1}\left(\frac{v}{2};\frac{1}{2};\frac{t^{2}}{2}\right)+\frac{t}{\sqrt{2}}\Gamma\left(\frac{1+\nu}{2}\right)\,_{1}F_{1}\left(\frac{1+\nu}{2};\frac{3}{2};\frac{t^{2}}{2}\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\sqrt{2}\Gamma\left(\frac{\nu+1}{2}\right)}{\Gamma\left(\frac{\nu}{2}\right)}\\ \mu_{2} & = & \nu-\mu^{2}\\ \gamma_{1} & = & \frac{2\mu^{3}+\mu\left(1-2\nu\right)}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{2\nu\left(1-\nu\right)-6\mu^{4}+4\mu^{2}\left(2\nu-1\right)}{\mu_{2}^{2}}\\ m_{d} & = & \sqrt{\nu-1}\quad\nu\geq1\\ m_{n} & = & \sqrt{2\Gamma^{-1}\left(\frac{\nu}{2},\frac{1}{2}\right)}\end{eqnarray*}

Implementation: `scipy.stats.chi`


Chi-squared
===========

This is the gamma distribution with :math:`L=0.0` and :math:`S=2.0` and :math:`\alpha=\nu/2` where :math:`\nu` is called the degrees of freedom. If :math:`Z_{1}\ldots Z_{\nu}` are all standard normal distributions, then :math:`W=\sum_{k}Z_{k}^{2}` has (standard) chi-square distribution with :math:`\nu` degrees of freedom.

The standard form (most often used in standard form only) is :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha\right) & = & \frac{1}{2\Gamma\left(\frac{\nu}{2}\right)}\left(\frac{x}{2}\right)^{\nu/2-1}e^{-x/2}\\ F\left(x;\alpha\right) & = & \Gamma\left(\frac{\nu}{2},\frac{x}{2}\right)\\ G\left(q;\alpha\right) & = & 2\Gamma^{-1}\left(\frac{\nu}{2},q\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\frac{\Gamma\left(\frac{\nu}{2}\right)}{\left(\frac{1}{2}-t\right)^{\nu/2}}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \nu\\ \mu_{2} & = & 2\nu\\ \gamma_{1} & = & \frac{2\sqrt{2}}{\sqrt{\nu}}\\ \gamma_{2} & = & \frac{12}{\nu}\\ m_{d} & = & \frac{\nu}{2}-1\end{eqnarray*}

Implementation: `scipy.stats.chi2`


Cosine
======

Approximation to the normal distribution.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{2\pi}\left[1+\cos x\right]I_{\left[-\pi,\pi\right]}\left(x\right)\\ F\left(x\right) & = & \frac{1}{2\pi}\left[\pi+x+\sin x\right]I_{\left[-\pi,\pi\right]}\left(x\right)+I_{\left(\pi,\infty\right)}\left(x\right)\\ G\left(\alpha\right) & = & F^{-1}\left(\alpha\right)\\ M\left(t\right) & = & \frac{\sinh\left(\pi t\right)}{\pi t\left(1+t^{2}\right)}\\ \mu=m_{d}=m_{n} & = & 0\\ \mu_{2} & = & \frac{\pi^{2}}{3}-2\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & \frac{-6\left(\pi^{4}-90\right)}{5\left(\pi^{2}-6\right)^{2}}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(4\pi\right)-1\\  & \approx & 1.5310242469692907930.\end{eqnarray*}

Implementation: `scipy.stats.cosine`


Double Gamma
============

The double gamma is the signed version of the Gamma distribution. For :math:`\alpha>0:`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha\right) & = & \frac{1}{2\Gamma\left(\alpha\right)}\left|x\right|^{\alpha-1}e^{-\left|x\right|}\\ F\left(x;\alpha\right) & = & \left\{ \begin{array}{ccc} \frac{1}{2}-\frac{1}{2}\Gamma\left(\alpha,\left|x\right|\right) &  & x\leq0\\ \frac{1}{2}+\frac{1}{2}\Gamma\left(\alpha,\left|x\right|\right) &  & x>0\end{array}\right.\\ G\left(q;\alpha\right) & = & \left\{ \begin{array}{ccc} -\Gamma^{-1}\left(\alpha,\left|2q-1\right|\right) &  & q\leq\frac{1}{2}\\ \Gamma^{-1}\left(\alpha,\left|2q-1\right|\right) &  & q>\frac{1}{2}\end{array}\right.\end{eqnarray*}

.. math::

     M\left(t\right)=\frac{1}{2\left(1-t\right)^{a}}+\frac{1}{2\left(1+t\right)^{a}}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu=m_{n} & = & 0\\ \mu_{2} & = & \alpha\left(\alpha+1\right)\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & \frac{\left(\alpha+2\right)\left(\alpha+3\right)}{\alpha\left(\alpha+1\right)}-3\\ m_{d} & = & \mathrm{NA}\end{eqnarray*}

Implementation: `scipy.stats.dgamma`


Double Weibull
==============

This is a signed form of the Weibull distribution.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{c}{2}\left|x\right|^{c-1}\exp\left(-\left|x\right|^{c}\right)\\ F\left(x;c\right) & = & \left\{ \begin{array}{ccc} \frac{1}{2}\exp\left(-\left|x\right|^{c}\right) &  & x\leq0\\ 1-\frac{1}{2}\exp\left(-\left|x\right|^{c}\right) &  & x>0\end{array}\right.\\ G\left(q;c\right) & = & \left\{ \begin{array}{ccc} -\log^{1/c}\left(\frac{1}{2q}\right) &  & q\leq\frac{1}{2}\\ \log^{1/c}\left(\frac{1}{2q-1}\right) &  & q>\frac{1}{2}\end{array}\right.\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\mu_{n}=\begin{cases} \Gamma\left(1+\frac{n}{c}\right) & n\mathrm{ even}\\ 0 & n\mathrm{ odd}\end{cases}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=\mu & = & 0\\ \mu_{2} & = & \Gamma\left(\frac{c+2}{c}\right)\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & \frac{\Gamma\left(1+\frac{4}{c}\right)}{\Gamma^{2}\left(1+\frac{2}{c}\right)}\\ m_{d} & = & \mathrm{NA bimodal}\end{eqnarray*}

Implementation: `scipy.stats.dweibull`


Erlang
======

This is just the Gamma distribution with shape parameter :math:`\alpha=n` an integer.

Implementation: `scipy.stats.erlang`


Exponential
===========

This is a special case of the Gamma (and Erlang) distributions with
shape parameter :math:`\left(\alpha=1\right)` and the same location and scale parameters. The standard form is
therefore ( :math:`x\geq0` )

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & e^{-x}\\ F\left(x\right) & = & \Gamma\left(1,x\right)=1-e^{-x}\\ G\left(q\right) & = & -\log\left(1-q\right)\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=n!

.. math::

     M\left(t\right)=\frac{1}{1-t}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 1\\ \mu_{2} & = & 1\\ \gamma_{1} & = & 2\\ \gamma_{2} & = & 6\\ m_{d} & = & 0\end{eqnarray*}

.. math::

     h\left[X\right]=1.

Implementation: `scipy.stats.expon`


Exponentiated Weibull
=====================

Two positive shape parameters :math:`a` and :math:`c` and :math:`x\in\left(0,\infty\right)`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,c\right) & = & ac\left[1-\exp\left(-x^{c}\right)\right]^{a-1}\exp\left(-x^{c}\right)x^{c-1}\\ F\left(x;a,c\right) & = & \left[1-\exp\left(-x^{c}\right)\right]^{a}\\ G\left(q;a,c\right) & = & \left[-\log\left(1-q^{1/a}\right)\right]^{1/c}\end{eqnarray*}

Implementation: `scipy.stats.exponweib`


Exponential Power
=================

One positive shape parameter :math:`b` . Defined for :math:`x\geq0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;b\right) & = & ebx^{b-1}\exp\left[x^{b}-e^{x^{b}}\right]\\ F\left(x;b\right) & = & 1-\exp\left[1-e^{x^{b}}\right]\\ G\left(q;b\right) & = & \log^{1/b}\left[1-\log\left(1-q\right)\right]\end{eqnarray*}

Implementation: `scipy.stats.exponpow`


Fatigue Life (Birnbaum-Saunders)
================================

This distribution's pdf is the average of the inverse-Gaussian :math:`\left(\mu=1\right)` and reciprocal inverse-Gaussian pdf :math:`\left(\mu=1\right)` . We follow the notation of JKB here with :math:`\beta=S.` for :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{x+1}{2c\sqrt{2\pi x^{3}}}\exp\left(-\frac{\left(x-1\right)^{2}}{2xc^{2}}\right)\\ F\left(x;c\right) & = & \Phi\left(\frac{1}{c}\left(\sqrt{x}-\frac{1}{\sqrt{x}}\right)\right)\\ G\left(q;c\right) & = & \frac{1}{4}\left[c\Phi^{-1}\left(q\right)+\sqrt{c^{2}\left(\Phi^{-1}\left(q\right)\right)^{2}+4}\right]^{2}\end{eqnarray*}

.. math::

     M\left(t\right)=c\sqrt{2\pi}\exp\left[\frac{1}{c^{2}}\left(1-\sqrt{1-2c^{2}t}\right)\right]\left(1+\frac{1}{\sqrt{1-2c^{2}t}}\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{c^{2}}{2}+1\\ \mu_{2} & = & c^{2}\left(\frac{5}{4}c^{2}+1\right)\\ \gamma_{1} & = & \frac{4c\sqrt{11c^{2}+6}}{\left(5c^{2}+4\right)^{3/2}}\\ \gamma_{2} & = & \frac{6c^{2}\left(93c^{2}+41\right)}{\left(5c^{2}+4\right)^{2}}\end{eqnarray*}

Implementation: `scipy.stats.fatiguelife`


Fisk (Log Logistic)
===================

Special case of the Burr distribution with :math:`d=1`

.. math::
   :nowrap:

    \begin{eqnarray*} c & > & 0\\ k & = & \Gamma\left(1-\frac{2}{c}\right)\Gamma\left(\frac{2}{c}+1\right)-\Gamma^{2}\left(1-\frac{1}{c}\right)\Gamma^{2}\left(\frac{1}{c}+1\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c,d\right) & = & \frac{cx^{c-1}}{\left(1+x^{c}\right)^{2}}I_{\left(0,\infty\right)}\left(x\right)\\ F\left(x;c,d\right) & = & \left(1+x^{-c}\right)^{-1}\\ G\left(\alpha;c,d\right) & = & \left(\alpha^{-1}-1\right)^{-1/c}\\ \mu & = & \Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+1\right)\\ \mu_{2} & = & k\\ \gamma_{1} & = & \frac{1}{\sqrt{k^{3}}}\left[2\Gamma^{3}\left(1-\frac{1}{c}\right)\Gamma^{3}\left(\frac{1}{c}+1\right)+\Gamma\left(1-\frac{3}{c}\right)\Gamma\left(\frac{3}{c}+1\right)\right.\\  &  & \left.-3\Gamma\left(1-\frac{2}{c}\right)\Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+1\right)\Gamma\left(\frac{2}{c}+1\right)\right]\\ \gamma_{2} & = & -3+\frac{1}{k^{2}}\left[6\Gamma\left(1-\frac{2}{c}\right)\Gamma^{2}\left(1-\frac{1}{c}\right)\Gamma^{2}\left(\frac{1}{c}+1\right)\Gamma\left(\frac{2}{c}+1\right)\right.\\  &  & -3\Gamma^{4}\left(1-\frac{1}{c}\right)\Gamma^{4}\left(\frac{1}{c}+1\right)+\Gamma\left(1-\frac{4}{c}\right)\Gamma\left(\frac{4}{c}+1\right)\\  &  & \left.-4\Gamma\left(1-\frac{3}{c}\right)\Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+1\right)\Gamma\left(\frac{3}{c}+1\right)\right]\\ m_{d} & = & \left(\frac{c-1}{c+1}\right)^{1/c}\,\mathrm{if }c>1\,\mathrm{otherwise }0\\ m_{n} & = & 1\end{eqnarray*}

.. math::

     h\left[X\right]=2-\log c.

Implementation: `scipy.stats.fisk`


Folded Cauchy
=============

This formula can be expressed in terms of the standard formulas for
the Cauchy distribution (call the cdf :math:`C\left(x\right)` and the pdf :math:`d\left(x\right)` ). if :math:`Y` is cauchy then :math:`\left|Y\right|` is folded cauchy. Note that :math:`x\geq0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{1}{\pi\left(1+\left(x-c\right)^{2}\right)}+\frac{1}{\pi\left(1+\left(x+c\right)^{2}\right)}\\ F\left(x;c\right) & = & \frac{1}{\pi}\tan^{-1}\left(x-c\right)+\frac{1}{\pi}\tan^{-1}\left(x+c\right)\\ G\left(q;c\right) & = & F^{-1}\left(x;c\right)\end{eqnarray*}

No moments

Implementation: `scipy.stats.foldcauchy`


Folded Normal
=============

If :math:`Z` is Normal with mean :math:`L` and :math:`\sigma=S` , then :math:`\left|Z\right|` is a folded normal with shape parameter :math:`c=\left|L\right|/S` , location parameter :math:`0` and scale parameter :math:`S` . This is a special case of the non-central chi distribution with one-
degree of freedom and non-centrality parameter :math:`c^{2}.` Note that :math:`c\geq0` . The standard form of the folded normal is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \sqrt{\frac{2}{\pi}}\cosh\left(cx\right)\exp\left(-\frac{x^{2}+c^{2}}{2}\right)\\ F\left(x;c\right) & = & \Phi\left(x-c\right)-\Phi\left(-x-c\right)=\Phi\left(x-c\right)+\Phi\left(x+c\right)-1\\ G\left(\alpha;c\right) & = & F^{-1}\left(x;c\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\exp\left[\frac{t}{2}\left(t-2c\right)\right]\left(1+e^{2ct}\right)

.. math::
   :nowrap:

    \begin{eqnarray*} k & = & \mathrm{erf}\left(\frac{c}{\sqrt{2}}\right)\\ p & = & \exp\left(-\frac{c^{2}}{2}\right)\\ \mu & = & \sqrt{\frac{2}{\pi}}p+ck\\ \mu_{2} & = & c^{2}+1-\mu^{2}\\ \gamma_{1} & = & \frac{\sqrt{\frac{2}{\pi}}p^{3}\left(4-\frac{\pi}{p^{2}}\left(2c^{2}+1\right)\right)+2ck\left(6p^{2}+3cpk\sqrt{2\pi}+\pi c\left(k^{2}-1\right)\right)}{\pi\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{c^{4}+6c^{2}+3+6\left(c^{2}+1\right)\mu^{2}-3\mu^{4}-4p\mu\left(\sqrt{\frac{2}{\pi}}\left(c^{2}+2\right)+\frac{ck}{p}\left(c^{2}+3\right)\right)}{\mu_{2}^{2}}\end{eqnarray*}

Implementation: `scipy.stats.foldnorm`


Fratio (or F)
=============

Defined for :math:`x>0` . The distribution of :math:`\left(X_{1}/X_{2}\right)\left(\nu_{2}/\nu_{1}\right)` if :math:`X_{1}` is chi-squared with :math:`v_{1}` degrees of freedom and :math:`X_{2}` is chi-squared with :math:`v_{2}` degrees of freedom.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\nu_{1},\nu_{2}\right) & = & \frac{\nu_{2}^{\nu_{2}/2}\nu_{1}^{\nu_{1}/2}x^{\nu_{1}/2-1}}{\left(\nu_{2}+\nu_{1}x\right)^{\left(\nu_{1}+\nu_{2}\right)/2}B\left(\frac{\nu_{1}}{2},\frac{\nu_{2}}{2}\right)}\\ F\left(x;v_{1},v_{2}\right) & = & I\left(\frac{\nu_{1}}{2},\frac{\nu_{2}}{2},\frac{\nu_{2}x}{\nu_{2}+\nu_{1}x}\right)\\ G\left(q;\nu_{1},\nu_{2}\right) & = & \left[\frac{\nu_{2}}{I^{-1}\left(\nu_{1}/2,\nu_{2}/2,q\right)}-\frac{\nu_{1}}{\nu_{2}}\right]^{-1}.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\nu_{2}}{\nu_{2}-2}\quad\nu_{2}>2\\ \mu_{2} & = & \frac{2\nu_{2}^{2}\left(\nu_{1}+\nu_{2}-2\right)}{\nu_{1}\left(\nu_{2}-2\right)^{2}\left(\nu_{2}-4\right)}\quad v_{2}>4\\ \gamma_{1} & = & \frac{2\left(2\nu_{1}+\nu_{2}-2\right)}{\nu_{2}-6}\sqrt{\frac{2\left(\nu_{2}-4\right)}{\nu_{1}\left(\nu_{1}+\nu_{2}-2\right)}}\quad\nu_{2}>6\\ \gamma_{2} & = & \frac{3\left[8+\left(\nu_{2}-6\right)\gamma_{1}^{2}\right]}{2\nu-16}\quad\nu_{2}>8\end{eqnarray*}

Implementation: `scipy.stats.f`


Fréchet (ExtremeLB, Extreme Value II, Weibull minimum)
=======================================================

A type of extreme-value distribution with a lower bound. Defined for :math:`x>0` and :math:`c>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & cx^{c-1}\exp\left(-x^{c}\right)\\ F\left(x;c\right) & = & 1-\exp\left(-x^{c}\right)\\ G\left(q;c\right) & = & \left[-\log\left(1-q\right)\right]^{1/c}\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\Gamma\left(1+\frac{n}{c}\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \Gamma\left(1+\frac{1}{c}\right)\\ \mu_{2} & = & \Gamma\left(1+\frac{2}{c}\right)-\Gamma^{2}\left(1-\frac{1}{c}\right)\\ \gamma_{1} & = & \frac{\Gamma\left(1+\frac{3}{c}\right)-3\Gamma\left(1+\frac{2}{c}\right)\Gamma\left(1+\frac{1}{c}\right)+2\Gamma^{3}\left(1+\frac{1}{c}\right)}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\Gamma\left(1+\frac{4}{c}\right)-4\Gamma\left(1+\frac{1}{c}\right)\Gamma\left(1+\frac{3}{c}\right)+6\Gamma^{2}\left(1+\frac{1}{c}\right)\Gamma\left(1+\frac{2}{c}\right)-\Gamma^{4}\left(1+\frac{1}{c}\right)}{\mu_{2}^{2}}-3\\ m_{d} & = & \left(\frac{c}{1+c}\right)^{1/c}\\ m_{n} & = & G\left(\frac{1}{2};c\right)\end{eqnarray*}

.. math::

     h\left[X\right]=-\frac{\gamma}{c}-\log\left(c\right)+\gamma+1

where :math:`\gamma` is Euler's constant and equal to

.. math::

     \gamma\approx0.57721566490153286061.

Implementation: `scipy.stats.frechet_r`


Fréchet (left-skewed, Extreme Value Type III, Weibull maximum)
===============================================================

Defined for :math:`x<0` and :math:`c>0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & c\left(-x\right)^{c-1}\exp\left(-\left(-x\right)^{c}\right)\\ F\left(x;c\right) & = & \exp\left(-\left(-x\right)^{c}\right)\\ G\left(q;c\right) & = & -\left(-\log q\right)^{1/c}\end{eqnarray*}

The mean is the negative of the right-skewed Frechet distribution
given above, and the other statistical parameters can be computed from

.. math::

     \mu_{n}^{\prime}=\left(-1\right)^{n}\Gamma\left(1+\frac{n}{c}\right).

.. math::

     h\left[X\right]=-\frac{\gamma}{c}-\log\left(c\right)+\gamma+1

where :math:`\gamma` is Euler's constant and equal to

.. math::

     \gamma\approx0.57721566490153286061.

Implementation: `scipy.stats.frechet_l`


Gamma
=====

The standard form for the gamma distribution is :math:`\left(\alpha>0\right)` valid for :math:`x\geq0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha\right) & = & \frac{1}{\Gamma\left(\alpha\right)}x^{\alpha-1}e^{-x}\\ F\left(x;\alpha\right) & = & \Gamma\left(\alpha,x\right)\\ G\left(q;\alpha\right) & = & \Gamma^{-1}\left(\alpha,q\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\frac{1}{\left(1-t\right)^{\alpha}}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \alpha\\ \mu_{2} & = & \alpha\\ \gamma_{1} & = & \frac{2}{\sqrt{\alpha}}\\ \gamma_{2} & = & \frac{6}{\alpha}\\ m_{d} & = & \alpha-1\end{eqnarray*}

.. math::

     h\left[X\right]=\Psi\left(a\right)\left[1-a\right]+a+\log\Gamma\left(a\right)

where

.. math::

     \Psi\left(a\right)=\frac{\Gamma^{\prime}\left(a\right)}{\Gamma\left(a\right)}.

Implementation: `scipy.stats.gamma`


Generalized Logistic
====================

Has been used in the analysis of extreme values. Has one shape
parameter :math:`c>0.` And :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{c\exp\left(-x\right)}{\left[1+\exp\left(-x\right)\right]^{c+1}}\\ F\left(x;c\right) & = & \frac{1}{\left[1+\exp\left(-x\right)\right]^{c}}\\ G\left(q;c\right) & = & -\log\left(q^{-1/c}-1\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\frac{c}{1-t}\,_{2}F_{1}\left(1+c,\,1-t\,;\,2-t\,;-1\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \gamma+\psi_{0}\left(c\right)\\ \mu_{2} & = & \frac{\pi^{2}}{6}+\psi_{1}\left(c\right)\\ \gamma_{1} & = & \frac{\psi_{2}\left(c\right)+2\zeta\left(3\right)}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\left(\frac{\pi^{4}}{15}+\psi_{3}\left(c\right)\right)}{\mu_{2}^{2}}\\ m_{d} & = & \log c\\ m_{n} & = & -\log\left(2^{1/c}-1\right)\end{eqnarray*}

Note that the polygamma function is

.. math::
   :nowrap:

    \begin{eqnarray*} \psi_{n}\left(z\right) & = & \frac{d^{n+1}}{dz^{n+1}}\log\Gamma\left(z\right)\\  & = & \left(-1\right)^{n+1}n!\sum_{k=0}^{\infty}\frac{1}{\left(z+k\right)^{n+1}}\\  & = & \left(-1\right)^{n+1}n!\zeta\left(n+1,z\right)\end{eqnarray*}

where :math:`\zeta\left(k,x\right)` is a generalization of the Riemann zeta function called the Hurwitz
zeta function Note that :math:`\zeta\left(n\right)\equiv\zeta\left(n,1\right)`

Implementation: `scipy.stats.genlogistic`


Generalized Pareto
==================

Shape parameter :math:`c\neq0` and defined for :math:`x\geq0` for all :math:`c` and :math:`x<\frac{1}{\left|c\right|}` if :math:`c` is negative.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \left(1+cx\right)^{-1-\frac{1}{c}}\\ F\left(x;c\right) & = & 1-\frac{1}{\left(1+cx\right)^{1/c}}\\ G\left(q;c\right) & = & \frac{1}{c}\left[\left(\frac{1}{1-q}\right)^{c}-1\right]\end{eqnarray*}

.. math::

     M\left(t\right)=\left\{ \begin{array}{cc} \left(-\frac{t}{c}\right)^{\frac{1}{c}}e^{-\frac{t}{c}}\left[\Gamma\left(1-\frac{1}{c}\right)+\Gamma\left(-\frac{1}{c},-\frac{t}{c}\right)-\pi\csc\left(\frac{\pi}{c}\right)/\Gamma\left(\frac{1}{c}\right)\right] & c>0\\ \left(\frac{\left|c\right|}{t}\right)^{1/\left|c\right|}\Gamma\left[\frac{1}{\left|c\right|},\frac{t}{\left|c\right|}\right] & c<0\end{array}\right.

.. math::

     \mu_{n}^{\prime}=\frac{\left(-1\right)^{n}}{c^{n}}\sum_{k=0}^{n}\left(\begin{array}{c} n\\ k\end{array}\right)\frac{\left(-1\right)^{k}}{1-ck}\quad cn<1

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{1}^{\prime} & = & \frac{1}{1-c}\quad c<1\\ \mu_{2}^{\prime} & = & \frac{2}{\left(1-2c\right)\left(1-c\right)}\quad c<\frac{1}{2}\\ \mu_{3}^{\prime} & = & \frac{6}{\left(1-c\right)\left(1-2c\right)\left(1-3c\right)}\quad c<\frac{1}{3}\\ \mu_{4}^{\prime} & = & \frac{24}{\left(1-c\right)\left(1-2c\right)\left(1-3c\right)\left(1-4c\right)}\quad c<\frac{1}{4}\end{eqnarray*}

Thus,

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu_{1}^{\prime}\\ \mu_{2} & = & \mu_{2}^{\prime}-\mu^{2}\\ \gamma_{1} & = & \frac{\mu_{3}^{\prime}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\mu_{4}^{\prime}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\end{eqnarray*}

.. math::

     h\left[X\right]=1+c\quad c>0

Implementation: `scipy.stats.genpareto`


Generalized Exponential
=======================

Three positive shape parameters for :math:`x\geq0.` Note that :math:`a,b,` and :math:`c` are all :math:`>0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,b,c\right) & = & \left(a+b\left(1-e^{-cx}\right)\right)\exp\left[ax-bx+\frac{b}{c}\left(1-e^{-cx}\right)\right]\\ F\left(x;a,b,c\right) & = & 1-\exp\left[ax-bx+\frac{b}{c}\left(1-e^{-cx}\right)\right]\\ G\left(q;a,b,c\right) & = & F^{-1}\end{eqnarray*}

Implementation: `scipy.stats.genexpon`


Generalized Extreme Value
=========================

Extreme value distributions with shape parameter :math:`c` .

For :math:`c>0` defined on :math:`-\infty<x\leq1/c.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \exp\left[-\left(1-cx\right)^{1/c}\right]\left(1-cx\right)^{1/c-1}\\ F\left(x;c\right) & = & \exp\left[-\left(1-cx\right)^{1/c}\right]\\ G\left(q;c\right) & = & \frac{1}{c}\left[1-\left(-\log q\right)^{c}\right]\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\frac{1}{c^{n}}\sum_{k=0}^{n}\left(\begin{array}{c} n\\ k\end{array}\right)\left(-1\right)^{k}\Gamma\left(ck+1\right)\quad cn>-1

So,

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{1}^{\prime} & = & \frac{1}{c}\left(1-\Gamma\left(1+c\right)\right)\quad c>-1\\ \mu_{2}^{\prime} & = & \frac{1}{c^{2}}\left(1-2\Gamma\left(1+c\right)+\Gamma\left(1+2c\right)\right)\quad c>-\frac{1}{2}\\ \mu_{3}^{\prime} & = & \frac{1}{c^{3}}\left(1-3\Gamma\left(1+c\right)+3\Gamma\left(1+2c\right)-\Gamma\left(1+3c\right)\right)\quad c>-\frac{1}{3}\\ \mu_{4}^{\prime} & = & \frac{1}{c^{4}}\left(1-4\Gamma\left(1+c\right)+6\Gamma\left(1+2c\right)-4\Gamma\left(1+3c\right)+\Gamma\left(1+4c\right)\right)\quad c>-\frac{1}{4}\end{eqnarray*}

For :math:`c<0` defined on :math:`\frac{1}{c}\leq x<\infty.` For :math:`c=0` defined over all space

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;0\right) & = & \exp\left[-e^{-x}\right]e^{-x}\\ F\left(x;0\right) & = & \exp\left[-e^{-x}\right]\\ G\left(q;0\right) & = & -\log\left(-\log q\right)\end{eqnarray*}

This is just the (left-skewed) Gumbel distribution for c=0.

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \gamma=-\psi_{0}\left(1\right)\\ \mu_{2} & = & \frac{\pi^{2}}{6}\\ \gamma_{1} & = & \frac{12\sqrt{6}}{\pi^{3}}\zeta\left(3\right)\\ \gamma_{2} & = & \frac{12}{5}\end{eqnarray*}

Implementation: `scipy.stats.genextreme`


Generalized Gamma
=================

A general probability form that reduces to many common distributions: :math:`x>0` :math:`a>0` and :math:`c\neq0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,c\right) & = & \frac{\left|c\right|x^{ca-1}}{\Gamma\left(a\right)}\exp\left(-x^{c}\right)\\ F\left(x;a,c\right) & = & \begin{array}{cc} \frac{\Gamma\left(a,x^{c}\right)}{\Gamma\left(a\right)} & c>0\\ 1-\frac{\Gamma\left(a,x^{c}\right)}{\Gamma\left(a\right)} & c<0\end{array}\\ G\left(q;a,c\right) & = & \left\{ \Gamma^{-1}\left[a,\Gamma\left(a\right)q\right]\right\} ^{1/c}\quad c>0\\  &  & \left\{ \Gamma^{-1}\left[a,\Gamma\left(a\right)\left(1-q\right)\right]\right\} ^{1/c}\quad c<0\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\frac{\Gamma\left(a+\frac{n}{c}\right)}{\Gamma\left(a\right)}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\Gamma\left(a+\frac{1}{c}\right)}{\Gamma\left(a\right)}\\ \mu_{2} & = & \frac{\Gamma\left(a+\frac{2}{c}\right)}{\Gamma\left(a\right)}-\mu^{2}\\ \gamma_{1} & = & \frac{\Gamma\left(a+\frac{3}{c}\right)/\Gamma\left(a\right)-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\Gamma\left(a+\frac{4}{c}\right)/\Gamma\left(a\right)-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\\ m_{d} & = & \left(\frac{ac-1}{c}\right)^{1/c}.\end{eqnarray*}

Special cases are Weibull :math:`\left(a=1\right)` , half-normal :math:`\left(a=1/2,c=2\right)` and ordinary gamma distributions :math:`c=1.` If :math:`c=-1` then it is the inverted gamma distribution.

.. math::

     h\left[X\right]=a-a\Psi\left(a\right)+\frac{1}{c}\Psi\left(a\right)+\log\Gamma\left(a\right)-\log\left|c\right|.

Implementation: `scipy.stats.gengamma`


Generalized Half-Logistic
=========================

For :math:`x\in\left[0,1/c\right]` and :math:`c>0` we have

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{2\left(1-cx\right)^{\frac{1}{c}-1}}{\left(1+\left(1-cx\right)^{1/c}\right)^{2}}\\ F\left(x;c\right) & = & \frac{1-\left(1-cx\right)^{1/c}}{1+\left(1-cx\right)^{1/c}}\\ G\left(q;c\right) & = & \frac{1}{c}\left[1-\left(\frac{1-q}{1+q}\right)^{c}\right]\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & 2-\left(2c+1\right)\log2.\end{eqnarray*}

Implementation: `scipy.stats.genhalflogistic`


Gilbrat
=======

Special case of the log-normal with :math:`\sigma=1` and :math:`S=1.0` (typically also :math:`L=0.0` )

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\sigma\right) & = & \frac{1}{x\sqrt{2\pi}}\exp\left[-\frac{1}{2}\left(\log x\right)^{2}\right]\\ F\left(x;\sigma\right) & = & \Phi\left(\log x\right)=\frac{1}{2}\left[1+\mathrm{erf}\left(\frac{\log x}{\sqrt{2}}\right)\right]\\ G\left(q;\sigma\right) & = & \exp\left\{ \Phi^{-1}\left(q\right)\right\} \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \sqrt{e}\\ \mu_{2} & = & e\left[e-1\right]\\ \gamma_{1} & = & \sqrt{e-1}\left(2+e\right)\\ \gamma_{2} & = & e^{4}+2e^{3}+3e^{2}-6\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(\sqrt{2\pi e}\right)\\  & \approx & 1.4189385332046727418\end{eqnarray*}

Implementation: `scipy.stats.gilbrat`


Gompertz (Truncated Gumbel)
===========================

For :math:`x\geq0` and :math:`c>0` . In JKB the two shape parameters :math:`b,a` are reduced to the single shape-parameter :math:`c=b/a` . As :math:`a` is just a scale parameter when :math:`a\neq0` . If :math:`a=0,` the distribution reduces to the exponential distribution scaled by :math:`1/b.` Thus, the standard form is given as

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & ce^{x}\exp\left[-c\left(e^{x}-1\right)\right]\\ F\left(x;c\right) & = & 1-\exp\left[-c\left(e^{x}-1\right)\right]\\ G\left(q;c\right) & = & \log\left[1-\frac{1}{c}\log\left(1-q\right)\right]\end{eqnarray*}

.. math::

     h\left[X\right]=1-\log\left(c\right)-e^{c}\mathrm{Ei}\left(1,c\right),

where

.. math::

     \mathrm{Ei}\left(n,x\right)=\int_{1}^{\infty}t^{-n}\exp\left(-xt\right)dt

Implementation: `scipy.stats.gompertz`


Gumbel (LogWeibull, Fisher-Tippetts, Type I Extreme Value)
==========================================================

One of a clase of extreme value distributions (right-skewed).

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \exp\left(-\left(x+e^{-x}\right)\right)\\ F\left(x\right) & = & \exp\left(-e^{-x}\right)\\ G\left(q\right) & = & -\log\left(-\log\left(q\right)\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\Gamma\left(1-t\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \gamma=-\psi_{0}\left(1\right)\\ \mu_{2} & = & \frac{\pi^{2}}{6}\\ \gamma_{1} & = & \frac{12\sqrt{6}}{\pi^{3}}\zeta\left(3\right)\\ \gamma_{2} & = & \frac{12}{5}\\ m_{d} & = & 0\\ m_{n} & = & -\log\left(\log2\right)\end{eqnarray*}

.. math::

     h\left[X\right]\approx1.0608407169541684911

Implementation: `scipy.stats.gumbel_r`


Gumbel Left-skewed (for minimum order statistic)
================================================

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \exp\left(x-e^{x}\right)\\ F\left(x\right) & = & 1-\exp\left(-e^{x}\right)\\ G\left(q\right) & = & \log\left(-\log\left(1-q\right)\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\Gamma\left(1+t\right)

Note, that :math:`\mu` is negative the mean for the right-skewed distribution. Similar for
median and mode. All other moments are the same.

.. math::

     h\left[X\right]\approx1.0608407169541684911.

Implementation: `scipy.stats.gumbel_l`


HalfCauchy
==========

If :math:`Z` is Hyperbolic Secant distributed then :math:`e^{Z}` is Half-Cauchy distributed. Also, if :math:`W` is (standard) Cauchy distributed, then :math:`\left|W\right|` is Half-Cauchy distributed. Special case of the Folded Cauchy
distribution with :math:`c=0.` The standard form is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{2}{\pi\left(1+x^{2}\right)}I_{[0,\infty)}\left(x\right)\\ F\left(x\right) & = & \frac{2}{\pi}\arctan\left(x\right)I_{\left[0,\infty\right]}\left(x\right)\\ G\left(q\right) & = & \tan\left(\frac{\pi}{2}q\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\cos t+\frac{2}{\pi}\left[\mathrm{Si}\left(t\right)\cos t-\mathrm{Ci}\left(\mathrm{-}t\right)\sin t\right]

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d} & = & 0\\ m_{n} & = & \tan\left(\frac{\pi}{4}\right)\end{eqnarray*}

No moments, as the integrals diverge.

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(2\pi\right)\\  & \approx & 1.8378770664093454836.\end{eqnarray*}

Implementation: `scipy.stats.halfcauchy`


HalfNormal
==========

This is a special case of the chi distribution with :math:`L=a` and :math:`S=b` and :math:`\nu=1.` This is also a special case of the folded normal with shape parameter :math:`c=0` and :math:`S=S.` If :math:`Z` is (standard) normally distributed then, :math:`\left|Z\right|` is half-normal. The standard form is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \sqrt{\frac{2}{\pi}}e^{-x^{2}/2}I_{\left(0,\infty\right)}\left(x\right)\\ F\left(x\right) & = & 2\Phi\left(x\right)-1\\ G\left(q\right) & = & \Phi^{-1}\left(\frac{1+q}{2}\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\sqrt{2\pi}e^{t^{2}/2}\Phi\left(t\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \sqrt{\frac{2}{\pi}}\\ \mu_{2} & = & 1-\frac{2}{\pi}\\ \gamma_{1} & = & \frac{\sqrt{2}\left(4-\pi\right)}{\left(\pi-2\right)^{3/2}}\\ \gamma_{2} & = & \frac{8\left(\pi-3\right)}{\left(\pi-2\right)^{2}}\\ m_{d} & = & 0\\ m_{n} & = & \Phi^{-1}\left(\frac{3}{4}\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(\sqrt{\frac{\pi e}{2}}\right)\\  & \approx & 0.72579135264472743239.\end{eqnarray*}

Implementation: `scipy.stats.halfnorm`


Half-Logistic
=============

In the limit as :math:`c\rightarrow\infty` for the generalized half-logistic we have the half-logistic defined
over :math:`x\geq0.` Also, the distribution of :math:`\left|X\right|` where :math:`X` has logistic distribtution.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{2e^{-x}}{\left(1+e^{-x}\right)^{2}}=\frac{1}{2}\mathrm{sech}^{2}\left(\frac{x}{2}\right)\\ F\left(x\right) & = & \frac{1-e^{-x}}{1+e^{-x}}=\tanh\left(\frac{x}{2}\right)\\ G\left(q\right) & = & \log\left(\frac{1+q}{1-q}\right)=2\mathrm{arctanh}\left(q\right)\end{eqnarray*}

.. math::

     M\left(t\right)=1-t\psi_{0}\left(\frac{1}{2}-\frac{t}{2}\right)+t\psi_{0}\left(1-\frac{t}{2}\right)

.. math::

     \mu_{n}^{\prime}=2\left(1-2^{1-n}\right)n!\zeta\left(n\right)\quad n\neq1

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{1}^{\prime} & = & 2\log\left(2\right)\\ \mu_{2}^{\prime} & = & 2\zeta\left(2\right)=\frac{\pi^{2}}{3}\\ \mu_{3}^{\prime} & = & 9\zeta\left(3\right)\\ \mu_{4}^{\prime} & = & 42\zeta\left(4\right)=\frac{7\pi^{4}}{15}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & 2-\log\left(2\right)\\  & \approx & 1.3068528194400546906.\end{eqnarray*}

Implementation: `scipy.stats.halflogistic`


Hyperbolic Secant
=================

Related to the logistic distribution and used in lifetime analysis.
Standard form is (defined over all :math:`x` )

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{\pi}\mathrm{sech}\left(x\right)\\ F\left(x\right) & = & \frac{2}{\pi}\arctan\left(e^{x}\right)\\ G\left(q\right) & = & \log\left(\tan\left(\frac{\pi}{2}q\right)\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\sec\left(\frac{\pi}{2}t\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{n}^{\prime} & = & \frac{1+\left(-1\right)^{n}}{2\pi2^{2n}}n!\left[\zeta\left(n+1,\frac{1}{4}\right)-\zeta\left(n+1,\frac{3}{4}\right)\right]\\  & = & \left\{ \begin{array}{cc} 0 & n\mathrm{ odd}\\ C_{n/2}\frac{\pi^{n}}{2^{n}} & n\mathrm{ even}\end{array}\right.\end{eqnarray*}

where :math:`C_{m}` is an integer given by

.. math::
   :nowrap:

    \begin{eqnarray*} C_{m} & = & \frac{\left(2m\right)!\left[\zeta\left(2m+1,\frac{1}{4}\right)-\zeta\left(2m+1,\frac{3}{4}\right)\right]}{\pi^{2m+1}2^{2m}}\\  & = & 4\left(-1\right)^{m-1}\frac{16^{m}}{2m+1}B_{2m+1}\left(\frac{1}{4}\right)\end{eqnarray*}

where :math:`B_{2m+1}\left(\frac{1}{4}\right)` is the Bernoulli polynomial of order :math:`2m+1` evaluated at :math:`1/4.` Thus

.. math::

     \mu_{n}^{\prime}=\left\{ \begin{array}{cc} 0 & n\mathrm{ odd}\\ 4\left(-1\right)^{n/2-1}\frac{\left(2\pi\right)^{n}}{n+1}B_{n+1}\left(\frac{1}{4}\right) & n\mathrm{ even}\end{array}\right.

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=m_{n}=\mu & = & 0\\ \mu_{2} & = & \frac{\pi^{2}}{4}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & 2\end{eqnarray*}

.. math::

     h\left[X\right]=\log\left(2\pi\right).

Implementation: `scipy.stats.hypsecant`


Gauss Hypergeometric
====================

:math:`x\in\left[0,1\right]` , :math:`\alpha>0,\,\beta>0`

.. math::

     C^{-1}=B\left(\alpha,\beta\right)\,_{2}F_{1}\left(\gamma,\alpha;\alpha+\beta;-z\right)

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha,\beta,\gamma,z\right) & = & Cx^{\alpha-1}\frac{\left(1-x\right)^{\beta-1}}{\left(1+zx\right)^{\gamma}}\\ \mu_{n}^{\prime} & = & \frac{B\left(n+\alpha,\beta\right)}{B\left(\alpha,\beta\right)}\frac{\,_{2}F_{1}\left(\gamma,\alpha+n;\alpha+\beta+n;-z\right)}{\,_{2}F_{1}\left(\gamma,\alpha;\alpha+\beta;-z\right)}\end{eqnarray*}

Implementation: `scipy.stats.gausshyper`


Inverted Gamma
==============

Special case of the generalized Gamma distribution with :math:`c=-1` and :math:`a>0` , :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a\right) & = & \frac{x^{-a-1}}{\Gamma\left(a\right)}\exp\left(-\frac{1}{x}\right)\\ F\left(x;a\right) & = & \frac{\Gamma\left(a,\frac{1}{x}\right)}{\Gamma\left(a\right)}\\ G\left(q;a\right) & = & \left\{ \Gamma^{-1}\left[a,\Gamma\left(a\right)q\right]\right\} ^{-1}\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\frac{\Gamma\left(a-n\right)}{\Gamma\left(a\right)}\quad a>n

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{1}{a-1}\quad a>1\\ \mu_{2} & = & \frac{1}{\left(a-2\right)\left(a-1\right)}-\mu^{2}\quad a>2\\ \gamma_{1} & = & \frac{\frac{1}{\left(a-3\right)\left(a-2\right)\left(a-1\right)}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\frac{1}{\left(a-4\right)\left(a-3\right)\left(a-2\right)\left(a-1\right)}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\end{eqnarray*}

.. math::

     m_{d}=\frac{1}{a+1}

.. math::

     h\left[X\right]=a-\left(a+1\right)\Psi\left(a\right)+\log\Gamma\left(a\right).

Implementation: `scipy.stats.invgamma`


Inverse Normal (Inverse Gaussian)
=================================

The standard form involves the shape parameter :math:`\mu` (in most definitions, :math:`L=0.0` is used). (In terms of the regress documentation :math:`\mu=A/B` ) and :math:`B=S` and :math:`L` is not a parameter in that distribution. A standard form is :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\mu\right) & = & \frac{1}{\sqrt{2\pi x^{3}}}\exp\left(-\frac{\left(x-\mu\right)^{2}}{2x\mu^{2}}\right).\\ F\left(x;\mu\right) & = & \Phi\left(\frac{1}{\sqrt{x}}\frac{x-\mu}{\mu}\right)+\exp\left(\frac{2}{\mu}\right)\Phi\left(-\frac{1}{\sqrt{x}}\frac{x+\mu}{\mu}\right)\\ G\left(q;\mu\right) & = & F^{-1}\left(q;\mu\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu\\ \mu_{2} & = & \mu^{3}\\ \gamma_{1} & = & 3\sqrt{\mu}\\ \gamma_{2} & = & 15\mu\\ m_{d} & = & \frac{\mu}{2}\left(\sqrt{9\mu^{2}+4}-3\mu\right)\end{eqnarray*}

This is related to the canonical form or JKB "two-parameter "inverse Gaussian when written in it's full form with scale parameter :math:`S` and location parameter :math:`L` by taking :math:`L=0` and :math:`S\equiv\lambda,` then :math:`\mu S` is equal to :math:`\mu_{2}` where :math:`\mu_{2}` is the parameter used by JKB. We prefer this form because of it's
consistent use of the scale parameter. Notice that in JKB the skew :math:`\left(\sqrt{\beta_{1}}\right)` and the kurtosis ( :math:`\beta_{2}-3` ) are both functions only of :math:`\mu_{2}/\lambda=\mu S/S=\mu` as shown here, while the variance and mean of the standard form here
are transformed appropriately.

Implementation: `scipy.stats.invgauss`


Inverted Weibull
================

Shape parameter :math:`c>0` and :math:`x>0` . Then

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & cx^{-c-1}\exp\left(-x^{-c}\right)\\ F\left(x;c\right) & = & \exp\left(-x^{-c}\right)\\ G\left(q;c\right) & = & \left(-\log q\right)^{-1/c}\end{eqnarray*}

.. math::

     h\left[X\right]=1+\gamma+\frac{\gamma}{c}-\log\left(c\right)

where :math:`\gamma` is Euler's constant.

Implementation: `scipy.stats.invweibull`


Johnson SB
==========

Defined for :math:`x\in\left(0,1\right)` with two shape parameters :math:`a` and :math:`b>0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,b\right) & = & \frac{b}{x\left(1-x\right)}\phi\left(a+b\log\frac{x}{1-x}\right)\\ F\left(x;a,b\right) & = & \Phi\left(a+b\log\frac{x}{1-x}\right)\\ G\left(q;a,b\right) & = & \frac{1}{1+\exp\left[-\frac{1}{b}\left(\Phi^{-1}\left(q\right)-a\right)\right]}\end{eqnarray*}

Implementation: `scipy.stats.johnsonsb`


Johnson SU
==========

Defined for all :math:`x` with two shape parameters :math:`a` and :math:`b>0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,b\right) & = & \frac{b}{\sqrt{x^{2}+1}}\phi\left(a+b\log\left(x+\sqrt{x^{2}+1}\right)\right)\\ F\left(x;a,b\right) & = & \Phi\left(a+b\log\left(x+\sqrt{x^{2}+1}\right)\right)\\ G\left(q;a,b\right) & = & \sinh\left[\frac{\Phi^{-1}\left(q\right)-a}{b}\right]\end{eqnarray*}

Implementation: `scipy.stats.johnsonsu`


KSone
=====

Implementation: `scipy.stats.ksone`


KStwo
=====

Implementation: `scipy.stats.kstwobign`


Laplace (Double Exponential, Bilateral Exponential)
===================================================

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{2}e^{-\left|x\right|}\\ F\left(x\right) & = & \left\{ \begin{array}{ccc} \frac{1}{2}e^{x} &  & x\leq0\\ 1-\frac{1}{2}e^{-x} &  & x>0\end{array}\right.\\ G\left(q\right) & = & \left\{ \begin{array}{ccc} \log\left(2q\right) &  & q\leq\frac{1}{2}\\ -\log\left(2-2q\right) &  & q>\frac{1}{2}\end{array}\right.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=m_{n}=\mu & = & 0\\ \mu_{2} & = & 2\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & 3\end{eqnarray*}

The ML estimator of the location parameter is

.. math::

     \hat{L}=\mathrm{median}\left(X_{i}\right)

where :math:`X_{i}` is a sequence of :math:`N` mutually independent Laplace RV's and the median is some number
between the :math:`\frac{1}{2}N\mathrm{th}` and the :math:`(N/2+1)\mathrm{th}` order statistic ( *e.g.* take the average of these two) when :math:`N` is even. Also,

.. math::

     \hat{S}=\frac{1}{N}\sum_{j=1}^{N}\left|X_{j}-\hat{L}\right|.

Replace :math:`\hat{L}` with :math:`L` if it is known. If :math:`L` is known then this estimator is distributed as :math:`\left(2N\right)^{-1}S\cdot\chi_{2N}^{2}` .

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(2e\right)\\  & \approx & 1.6931471805599453094.\end{eqnarray*}

Implementation: `scipy.stats.laplace`


Left-skewed Lévy
=================

Special case of Lévy-stable distribution with :math:`\alpha=\frac{1}{2}` and :math:`\beta=-1` the support is :math:`x<0` . In standard form

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{\left|x\right|\sqrt{2\pi\left|x\right|}}\exp\left(-\frac{1}{2\left|x\right|}\right)\\ F\left(x\right) & = & 2\Phi\left(\frac{1}{\sqrt{\left|x\right|}}\right)-1\\ G\left(q\right) & = & -\left[\Phi^{-1}\left(\frac{q+1}{2}\right)\right]^{-2}.\end{eqnarray*}

No moments.

Implementation: `scipy.stats.levy_l`


Lévy
=====

A special case of Lévy-stable distributions with :math:`\alpha=\frac{1}{2}` and :math:`\beta=1` . In standard form it is defined for :math:`x>0` as

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{x\sqrt{2\pi x}}\exp\left(-\frac{1}{2x}\right)\\ F\left(x\right) & = & 2\left[1-\Phi\left(\frac{1}{\sqrt{x}}\right)\right]\\ G\left(q\right) & = & \left[\Phi^{-1}\left(1-\frac{q}{2}\right)\right]^{-2}.\end{eqnarray*}

It has no finite moments.

Implementation: `scipy.stats.levy`


Logistic (Sech-squared)
=======================

A special case of the Generalized Logistic distribution with :math:`c=1.` Defined for :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{\exp\left(-x\right)}{\left[1+\exp\left(-x\right)\right]^{2}}\\ F\left(x\right) & = & \frac{1}{1+\exp\left(-x\right)}\\ G\left(q\right) & = & -\log\left(1/q-1\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \gamma+\psi_{0}\left(1\right)=0\\ \mu_{2} & = & \frac{\pi^{2}}{6}+\psi_{1}\left(1\right)=\frac{\pi^{2}}{3}\\ \gamma_{1} & = & \frac{\psi_{2}\left(c\right)+2\zeta\left(3\right)}{\mu_{2}^{3/2}}=0\\ \gamma_{2} & = & \frac{\left(\frac{\pi^{4}}{15}+\psi_{3}\left(c\right)\right)}{\mu_{2}^{2}}=\frac{6}{5}\\ m_{d} & = & \log1=0\\ m_{n} & = & -\log\left(2-1\right)=0\end{eqnarray*}

.. math::

     h\left[X\right]=1.

Implementation: `scipy.stats.logistic`


Log Double Exponential (Log-Laplace)
====================================

Defined over :math:`x>0` with :math:`c>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \left\{ \begin{array}{ccc} \frac{c}{2}x^{c-1} &  & 0<x<1\\ \frac{c}{2}x^{-c-1} &  & x\geq1\end{array}\right.\\ F\left(x;c\right) & = & \left\{ \begin{array}{ccc} \frac{1}{2}x^{c} &  & 0<x<1\\ 1-\frac{1}{2}x^{-c} &  & x\geq1\end{array}\right.\\ G\left(q;c\right) & = & \left\{ \begin{array}{ccc} \left(2q\right)^{1/c} &  & 0\leq q<\frac{1}{2}\\ \left(2-2q\right)^{-1/c} &  & \frac{1}{2}\leq q\leq1\end{array}\right.\end{eqnarray*}

.. math::

     h\left[X\right]=\log\left(\frac{2e}{c}\right)


Implementation: `scipy.stats.loglaplace`


Log Gamma
=========

A single shape parameter :math:`c>0` (Defined for all :math:`x` )

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{\exp\left(cx-e^{x}\right)}{\Gamma\left(c\right)}\\ F\left(x;c\right) & = & \frac{\Gamma\left(c,e^{x}\right)}{\Gamma\left(c\right)}\\ G\left(q;c\right) & = & \log\left[\Gamma^{-1}\left[c,q\Gamma\left(c\right)\right]\right]\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\int_{0}^{\infty}\left[\log y\right]^{n}y^{c-1}\exp\left(-y\right)dy.

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu_{1}^{\prime}\\ \mu_{2} & = & \mu_{2}^{\prime}-\mu^{2}\\ \gamma_{1} & = & \frac{\mu_{3}^{\prime}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\mu_{4}^{\prime}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\end{eqnarray*}

Implementation: `scipy.stats.loggamma`


Log Normal (Cobb-Douglass)
==========================

Has one shape parameter :math:`\sigma` >0. (Notice that the "Regress ":math:`A=\log S` where :math:`S` is the scale parameter and :math:`A` is the mean of the underlying normal distribution). The standard form
is :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\sigma\right) & = & \frac{1}{\sigma x\sqrt{2\pi}}\exp\left[-\frac{1}{2}\left(\frac{\log x}{\sigma}\right)^{2}\right]\\ F\left(x;\sigma\right) & = & \Phi\left(\frac{\log x}{\sigma}\right)\\ G\left(q;\sigma\right) & = & \exp\left\{ \sigma\Phi^{-1}\left(q\right)\right\} \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \exp\left(\sigma^{2}/2\right)\\ \mu_{2} & = & \exp\left(\sigma^{2}\right)\left[\exp\left(\sigma^{2}\right)-1\right]\\ \gamma_{1} & = & \sqrt{p-1}\left(2+p\right)\\ \gamma_{2} & = & p^{4}+2p^{3}+3p^{2}-6\quad\quad p=e^{\sigma^{2}}\end{eqnarray*}

Notice that using JKB notation we have :math:`\theta=L,` :math:`\zeta=\log S` and we have given the so-called antilognormal form of the
distribution. This is more consistent with the location, scale
parameter description of general probability distributions.

.. math::

     h\left[X\right]=\frac{1}{2}\left[1+\log\left(2\pi\right)+2\log\left(\sigma\right)\right].

Also, note that if :math:`X` is a log-normally distributed random-variable with :math:`L=0` and :math:`S` and shape parameter :math:`\sigma.` Then, :math:`\log X` is normally distributed with variance :math:`\sigma^{2}` and mean :math:`\log S.`

Implementation: `scipy.stats.lognorm`


Maxwell
=======

This is a special case of the Chi distribution with :math:`L=0` and :math:`S=S=\frac{1}{\sqrt{a}}` and :math:`\nu=3.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \sqrt{\frac{2}{\pi}}x^{2}e^{-x^{2}/2}I_{\left(0,\infty\right)}\left(x\right)\\ F\left(x\right) & = & \Gamma\left(\frac{3}{2},\frac{x^{2}}{2}\right)\\ G\left(\alpha\right) & = & \sqrt{2\Gamma^{-1}\left(\frac{3}{2},\alpha\right)}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 2\sqrt{\frac{2}{\pi}}\\ \mu_{2} & = & 3-\frac{8}{\pi}\\ \gamma_{1} & = & \sqrt{2}\frac{32-10\pi}{\left(3\pi-8\right)^{3/2}}\\ \gamma_{2} & = & \frac{-12\pi^{2}+160\pi-384}{\left(3\pi-8\right)^{2}}\\ m_{d} & = & \sqrt{2}\\ m_{n} & = & \sqrt{2\Gamma^{-1}\left(\frac{3}{2},\frac{1}{2}\right)}\end{eqnarray*}

.. math::

     h\left[X\right]=\log\left(\sqrt{\frac{2\pi}{e}}\right)+\gamma.

Implementation: `scipy.stats.maxwell`


Mielke's Beta-Kappa
===================

A generalized F distribution. Two shape parameters :math:`\kappa` and :math:`\theta` , and :math:`x>0` . The :math:`\beta` in the DATAPLOT reference is a scale parameter.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\kappa,\theta\right) & = & \frac{\kappa x^{\kappa-1}}{\left(1+x^{\theta}\right)^{1+\frac{\kappa}{\theta}}}\\ F\left(x;\kappa,\theta\right) & = & \frac{x^{\kappa}}{\left(1+x^{\theta}\right)^{\kappa/\theta}}\\ G\left(q;\kappa,\theta\right) & = & \left(\frac{q^{\theta/\kappa}}{1-q^{\theta/\kappa}}\right)^{1/\theta}\end{eqnarray*}

Implementation: `scipy.stats.mielke`


Nakagami
========

Generalization of the chi distribution. Shape parameter is :math:`\nu>0.` Defined for :math:`x>0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\nu\right) & = & \frac{2\nu^{\nu}}{\Gamma\left(\nu\right)}x^{2\nu-1}\exp\left(-\nu x^{2}\right)\\ F\left(x;\nu\right) & = & \Gamma\left(\nu,\nu x^{2}\right)\\ G\left(q;\nu\right) & = & \sqrt{\frac{1}{\nu}\Gamma^{-1}\left(v,q\right)}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\Gamma\left(\nu+\frac{1}{2}\right)}{\sqrt{\nu}\Gamma\left(\nu\right)}\\ \mu_{2} & = & \left[1-\mu^{2}\right]\\ \gamma_{1} & = & \frac{\mu\left(1-4v\mu_{2}\right)}{2\nu\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{-6\mu^{4}\nu+\left(8\nu-2\right)\mu^{2}-2\nu+1}{\nu\mu_{2}^{2}}\end{eqnarray*}

Implementation: `scipy.stats.nakagami`


Noncentral chi-squared
======================

The distribution of :math:`\sum_{i=1}^{\nu}\left(Z_{i}+\delta_{i}\right)^{2}` where :math:`Z_{i}` are independent standard normal variables and :math:`\delta_{i}` are constants. :math:`\lambda=\sum_{i=1}^{\nu}\delta_{i}^{2}>0.` (In communications it is called the Marcum-Q function). Can be thought
of as a Generalized Rayleigh-Rice distribution. For :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\nu,\lambda\right) & = & e^{-\left(\lambda+x\right)/2}\frac{1}{2}\left(\frac{x}{\lambda}\right)^{\left(\nu-2\right)/4}I_{\left(\nu-2\right)/2}\left(\sqrt{\lambda x}\right)\\ F\left(x;\nu,\lambda\right) & = & \sum_{j=0}^{\infty}\left\{ \frac{\left(\lambda/2\right)^{j}}{j!}e^{-\lambda/2}\right\} \mathrm{Pr}\left[\chi_{\nu+2j}^{2}\leq x\right]\\ G\left(q;\nu,\lambda\right) & = & F^{-1}\left(x;\nu,\lambda\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \nu+\lambda\\ \mu_{2} & = & 2\left(\nu+2\lambda\right)\\ \gamma_{1} & = & \frac{\sqrt{8}\left(\nu+3\lambda\right)}{\left(\nu+2\lambda\right)^{3/2}}\\ \gamma_{2} & = & \frac{12\left(\nu+4\lambda\right)}{\left(\nu+2\lambda\right)^{2}}\end{eqnarray*}

Implementation: `scipy.stats.ncx2`


Noncentral F
============

Let :math:`\lambda>0` and :math:`\nu_{1}>0` and :math:`\nu_{2}>0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\lambda,\nu_{1},\nu_{2}\right) & = & \exp\left[\frac{\lambda}{2}+\frac{\left(\lambda\nu_{1}x\right)}{2\left(\nu_{1}x+\nu_{2}\right)}\right]\nu_{1}^{\nu_{1}/2}\nu_{2}^{\nu_{2}/2}x^{\nu_{1}/2-1}\\  &  & \times\left(\nu_{2}+\nu_{1}x\right)^{-\left(\nu_{1}+\nu_{2}\right)/2}\frac{\Gamma\left(\frac{\nu_{1}}{2}\right)\Gamma\left(1+\frac{\nu_{2}}{2}\right)L_{\nu_{2}/2}^{\nu_{1}/2-1}\left(-\frac{\lambda\nu_{1}x}{2\left(\nu_{1}x+\nu_{2}\right)}\right)}{B\left(\frac{\nu_{1}}{2},\frac{\nu_{2}}{2}\right)\Gamma\left(\frac{\nu_{1}+\nu_{2}}{2}\right)}\end{eqnarray*}

Implementation: `scipy.stats.ncf`


Noncentral t
============

The distribution of the ratio

.. math::

     \frac{U+\lambda}{\chi_{\nu}/\sqrt{\nu}}

where :math:`U` and :math:`\chi_{\nu}` are independent and distributed as a standard normal and chi with :math:`\nu` degrees of freedom. Note :math:`\lambda>0` and :math:`\nu>0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\lambda,\nu\right) & = & \frac{\nu^{\nu/2}\Gamma\left(\nu+1\right)}{2^{\nu}e^{\lambda^{2}/2}\left(\nu+x^{2}\right)^{\nu/2}\Gamma\left(\nu/2\right)}\\  &  & \times\left\{ \frac{\sqrt{2}\lambda x\,_{1}F_{1}\left(\frac{\nu}{2}+1;\frac{3}{2};\frac{\lambda^{2}x^{2}}{2\left(\nu+x^{2}\right)}\right)}{\left(\nu+x^{2}\right)\Gamma\left(\frac{\nu+1}{2}\right)}\right.\\  &  & -\left.\frac{\,_{1}F_{1}\left(\frac{\nu+1}{2};\frac{1}{2};\frac{\lambda^{2}x^{2}}{2\left(\nu+x^{2}\right)}\right)}{\sqrt{\nu+x^{2}}\Gamma\left(\frac{\nu}{2}+1\right)}\right\} \\  & = & \frac{\Gamma\left(\nu+1\right)}{2^{\left(\nu-1\right)/2}\sqrt{\pi\nu}\Gamma\left(\nu/2\right)}\exp\left[-\frac{\nu\lambda^{2}}{\nu+x^{2}}\right]\\  &  & \times\left(\frac{\nu}{\nu+x^{2}}\right)^{\left(\nu-1\right)/2}Hh_{\nu}\left(-\frac{\lambda x}{\sqrt{\nu+x^{2}}}\right)\\ F\left(x;\lambda,\nu\right) & =\end{eqnarray*}

Implementation: `scipy.stats.nct`


Normal
======

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{e^{-x^{2}/2}}{\sqrt{2\pi}}\\ F\left(x\right) & = & \Phi\left(x\right)=\frac{1}{2}+\frac{1}{2}\mathrm{erf}\left(\frac{\mathrm{x}}{\sqrt{2}}\right)\\ G\left(q\right) & = & \Phi^{-1}\left(q\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=m_{n}=\mu & = & 0\\ \mu_{2} & = & 1\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & 0\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(\sqrt{2\pi e}\right)\\  & \approx & 1.4189385332046727418\end{eqnarray*}

Implementation: `scipy.stats.norm`


Pareto
======

For :math:`x\geq1` and :math:`b>0` . Standard form is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;b\right) & = & \frac{b}{x^{b+1}}\\ F\left(x;b\right) & = & 1-\frac{1}{x^{b}}\\ G\left(q;b\right) & = & \left(1-q\right)^{-1/b}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{b}{b-1}\quad b>1\\ \mu_{2} & = & \frac{b}{\left(b-2\right)\left(b-1\right)^{2}}\quad b>2\\ \gamma_{1} & = & \frac{2\left(b+1\right)\sqrt{b-2}}{\left(b-3\right)\sqrt{b}}\quad b>3\\ \gamma_{2} & = & \frac{6\left(b^{3}+b^{2}-6b-2\right)}{b\left(b^{2}-7b+12\right)}\quad b>4\end{eqnarray*}

.. math::

     h\left(X\right)=\frac{1}{c}+1-\log\left(c\right)

Implementation: `scipy.stats.pareto`


Pareto Second Kind (Lomax)
==========================

:math:`c>0.` This is Pareto of the first kind with :math:`L=-1.0` so :math:`x\geq0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{c}{\left(1+x\right)^{c+1}}\\ F\left(x;c\right) & = & 1-\frac{1}{\left(1+x\right)^{c}}\\ G\left(q;c\right) & = & \left(1-q\right)^{-1/c}-1\end{eqnarray*}

.. math::

     h\left[X\right]=\frac{1}{c}+1-\log\left(c\right).

Implementation: `scipy.stats.lomax`


Power Log Normal
================

A generalization of the log-normal distribution :math:`\sigma>0` and :math:`c>0` and :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\sigma,c\right) & = & \frac{c}{x\sigma}\phi\left(\frac{\log x}{\sigma}\right)\left(\Phi\left(-\frac{\log x}{\sigma}\right)\right)^{c-1}\\ F\left(x;\sigma,c\right) & = & 1-\left(\Phi\left(-\frac{\log x}{\sigma}\right)\right)^{c}\\ G\left(q;\sigma,c\right) & = & \exp\left[-\sigma\Phi^{-1}\left[\left(1-q\right)^{1/c}\right]\right]\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\int_{0}^{1}\exp\left[-n\sigma\Phi^{-1}\left(y^{1/c}\right)\right]dy

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu_{1}^{\prime}\\ \mu_{2} & = & \mu_{2}^{\prime}-\mu^{2}\\ \gamma_{1} & = & \frac{\mu_{3}^{\prime}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\mu_{4}^{\prime}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\end{eqnarray*}

This distribution reduces to the log-normal distribution when :math:`c=1.`

Implementation: `scipy.stats.powerlognorm`


Power Normal
============

A generalization of the normal distribution, :math:`c>0` for

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & c\phi\left(x\right)\left(\Phi\left(-x\right)\right)^{c-1}\\ F\left(x;c\right) & = & 1-\left(\Phi\left(-x\right)\right)^{c}\\ G\left(q;c\right) & = & -\Phi^{-1}\left[\left(1-q\right)^{1/c}\right]\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\left(-1\right)^{n}\int_{0}^{1}\left[\Phi^{-1}\left(y^{1/c}\right)\right]^{n}dy

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu_{1}^{\prime}\\ \mu_{2} & = & \mu_{2}^{\prime}-\mu^{2}\\ \gamma_{1} & = & \frac{\mu_{3}^{\prime}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\mu_{4}^{\prime}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\end{eqnarray*}

For :math:`c=1` this reduces to the normal distribution.

Implementation: `scipy.stats.powernorm`


Power-function
==============

A special case of the beta distribution with :math:`b=1` : defined for :math:`x\in\left[0,1\right]`

.. math::

     a>0

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a\right) & = & ax^{a-1}\\ F\left(x;a\right) & = & x^{a}\\ G\left(q;a\right) & = & q^{1/a}\\ \mu & = & \frac{a}{a+1}\\ \mu_{2} & = & \frac{a\left(a+2\right)}{\left(a+1\right)^{2}}\\ \gamma_{1} & = & 2\left(1-a\right)\sqrt{\frac{a+2}{a\left(a+3\right)}}\\ \gamma_{2} & = & \frac{6\left(a^{3}-a^{2}-6a+2\right)}{a\left(a+3\right)\left(a+4\right)}\\ m_{d} & = & 1\end{eqnarray*}

.. math::

     h\left[X\right]=1-\frac{1}{a}-\log\left(a\right)

Implementation: `scipy.stats.powerlaw`


R-distribution
==============

A general-purpose distribution with a variety of shapes controlled by :math:`c>0.` Range of standard distribution is :math:`x\in\left[-1,1\right]`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{\left(1-x^{2}\right)^{c/2-1}}{B\left(\frac{1}{2},\frac{c}{2}\right)}\\ F\left(x;c\right) & = & \frac{1}{2}+\frac{x}{B\left(\frac{1}{2},\frac{c}{2}\right)}\,_{2}F_{1}\left(\frac{1}{2},1-\frac{c}{2};\frac{3}{2};x^{2}\right)\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\frac{\left(1+\left(-1\right)^{n}\right)}{2}B\left(\frac{n+1}{2},\frac{c}{2}\right)

The R-distribution with parameter :math:`n` is the distribution of the correlation coefficient of a random sample
of size :math:`n` drawn from a bivariate normal distribution with :math:`\rho=0.` The mean of the standard distribution is always zero and as the sample
size grows, the distribution's mass concentrates more closely about
this mean.

Implementation: `scipy.stats.rdist`


Rayleigh
========

This is Chi distribution with :math:`L=0.0` and :math:`\nu=2` and :math:`S=S` (no location parameter is generally used), the mode of the
distribution is :math:`S.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(r\right) & = & re^{-r^{2}/2}I_{[0,\infty)}\left(x\right)\\ F\left(r\right) & = & 1-e^{-r^{2}/2}I_{[0,\infty)}\left(x\right)\\ G\left(q\right) & = & \sqrt{-2\log\left(1-q\right)}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \sqrt{\frac{\pi}{2}}\\ \mu_{2} & = & \frac{4-\pi}{2}\\ \gamma_{1} & = & \frac{2\left(\pi-3\right)\sqrt{\pi}}{\left(4-\pi\right)^{3/2}}\\ \gamma_{2} & = & \frac{24\pi-6\pi^{2}-16}{\left(4-\pi\right)^{2}}\\ m_{d} & = & 1\\ m_{n} & = & \sqrt{2\log\left(2\right)}\end{eqnarray*}

.. math::

     h\left[X\right]=\frac{\gamma}{2}+\log\left(\frac{e}{\sqrt{2}}\right).

.. math::

     \mu_{n}^{\prime}=\sqrt{2^{n}}\Gamma\left(\frac{n}{2}+1\right)

Implementation: `scipy.stats.rayleigh`


Rice
====

Defined for :math:`x>0` and :math:`b>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;b\right) & = & x\exp\left(-\frac{x^{2}+b^{2}}{2}\right)I_{0}\left(xb\right)\\ F\left(x;b\right) & = & \int_{0}^{x}\alpha\exp\left(-\frac{\alpha^{2}+b^{2}}{2}\right)I_{0}\left(\alpha b\right)d\alpha\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\sqrt{2^{n}}\Gamma\left(1+\frac{n}{2}\right)\,_{1}F_{1}\left(-\frac{n}{2};1;-\frac{b^{2}}{2}\right)

Implementation: `scipy.stats.rice`


Reciprocal
==========

Shape parameters :math:`a,b>0` :math:`x\in\left[a,b\right]`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,b\right) & = & \frac{1}{x\log\left(b/a\right)}\\ F\left(x;a,b\right) & = & \frac{\log\left(x/a\right)}{\log\left(b/a\right)}\\ G\left(q;a,b\right) & = & a\exp\left(q\log\left(b/a\right)\right)=a\left(\frac{b}{a}\right)^{q}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} d & = & \log\left(a/b\right)\\ \mu & = & \frac{a-b}{d}\\ \mu_{2} & = & \mu\frac{a+b}{2}-\mu^{2}=\frac{\left(a-b\right)\left[a\left(d-2\right)+b\left(d+2\right)\right]}{2d^{2}}\\ \gamma_{1} & = & \frac{\sqrt{2}\left[12d\left(a-b\right)^{2}+d^{2}\left(a^{2}\left(2d-9\right)+2abd+b^{2}\left(2d+9\right)\right)\right]}{3d\sqrt{a-b}\left[a\left(d-2\right)+b\left(d+2\right)\right]^{3/2}}\\ \gamma_{2} & = & \frac{-36\left(a-b\right)^{3}+36d\left(a-b\right)^{2}\left(a+b\right)-16d^{2}\left(a^{3}-b^{3}\right)+3d^{3}\left(a^{2}+b^{2}\right)\left(a+b\right)}{3\left(a-b\right)\left[a\left(d-2\right)+b\left(d+2\right)\right]^{2}}-3\\ m_{d} & = & a\\ m_{n} & = & \sqrt{ab}\end{eqnarray*}

.. math::

     h\left[X\right]=\frac{1}{2}\log\left(ab\right)+\log\left[\log\left(\frac{b}{a}\right)\right].

Implementation: `scipy.stats.reciprocal`


Reciprocal Inverse Gaussian
===========================

The pdf is found from the inverse gaussian (IG), :math:`f_{RIG}\left(x;\mu\right)=\frac{1}{x^{2}}f_{IG}\left(\frac{1}{x};\mu\right)` defined for :math:`x\geq0` as

.. math::
   :nowrap:

    \begin{eqnarray*} f_{IG}\left(x;\mu\right) & = & \frac{1}{\sqrt{2\pi x^{3}}}\exp\left(-\frac{\left(x-\mu\right)^{2}}{2x\mu^{2}}\right).\\ F_{IG}\left(x;\mu\right) & = & \Phi\left(\frac{1}{\sqrt{x}}\frac{x-\mu}{\mu}\right)+\exp\left(\frac{2}{\mu}\right)\Phi\left(-\frac{1}{\sqrt{x}}\frac{x+\mu}{\mu}\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} f_{RIG}\left(x;\mu\right) & = & \frac{1}{\sqrt{2\pi x}}\exp\left(-\frac{\left(1-\mu x\right)^{2}}{2x\mu^{2}}\right)\\ F_{RIG}\left(x;\mu\right) & = & 1-F_{IG}\left(\frac{1}{x},\mu\right)\\  & = & 1-\Phi\left(\frac{1}{\sqrt{x}}\frac{1-\mu x}{\mu}\right)-\exp\left(\frac{2}{\mu}\right)\Phi\left(-\frac{1}{\sqrt{x}}\frac{1+\mu x}{\mu}\right)\end{eqnarray*}

Implementation: `scipy.stats.recipinvgauss`


Semicircular
============

Defined on :math:`x\in\left[-1,1\right]`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{2}{\pi}\sqrt{1-x^{2}}\\ F\left(x\right) & = & \frac{1}{2}+\frac{1}{\pi}\left[x\sqrt{1-x^{2}}+\arcsin x\right]\\ G\left(q\right) & = & F^{-1}\left(q\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=m_{n}=\mu & = & 0\\ \mu_{2} & = & \frac{1}{4}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & -1\end{eqnarray*}

.. math::

     h\left[X\right]=0.64472988584940017414.

Implementation: `scipy.stats.semicircular`


Student t
=========

Shape parameter :math:`\nu>0.` :math:`I\left(a,b,x\right)` is the incomplete beta integral and :math:`I^{-1}\left(a,b,I\left(a,b,x\right)\right)=x`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\nu\right) & = & \frac{\Gamma\left(\frac{\nu+1}{2}\right)}{\sqrt{\pi\nu}\Gamma\left(\frac{\nu}{2}\right)\left[1+\frac{x^{2}}{\nu}\right]^{\frac{\nu+1}{2}}}\\ F\left(x;\nu\right) & = & \left\{ \begin{array}{ccc} \frac{1}{2}I\left(\frac{\nu}{2},\frac{1}{2},\frac{\nu}{\nu+x^{2}}\right) &  & x\leq0\\ 1-\frac{1}{2}I\left(\frac{\nu}{2},\frac{1}{2},\frac{\nu}{\nu+x^{2}}\right) &  & x\geq0\end{array}\right.\\ G\left(q;\nu\right) & = & \left\{ \begin{array}{ccc} -\sqrt{\frac{\nu}{I^{-1}\left(\frac{\nu}{2},\frac{1}{2},2q\right)}-\nu} &  & q\leq\frac{1}{2}\\ \sqrt{\frac{\nu}{I^{-1}\left(\frac{\nu}{2},\frac{1}{2},2-2q\right)}-\nu} &  & q\geq\frac{1}{2}\end{array}\right.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{n}=m_{d}=\mu & = & 0\\ \mu_{2} & = & \frac{\nu}{\nu-2}\quad\nu>2\\ \gamma_{1} & = & 0\quad\nu>3\\ \gamma_{2} & = & \frac{6}{\nu-4}\quad\nu>4\end{eqnarray*}

As :math:`\nu\rightarrow\infty,` this distribution approaches the standard normal distribution.

.. math::

     h\left[X\right]=\frac{1}{4}\log\left(\frac{\pi c\Gamma^{2}\left(\frac{c}{2}\right)}{\Gamma^{2}\left(\frac{c+1}{2}\right)}\right)-\frac{\left(c+1\right)}{4}\left[\Psi\left(\frac{c}{2}\right)-cZ\left(c\right)+\pi\tan\left(\frac{\pi c}{2}\right)+\gamma+2\log2\right]

where

.. math::

     Z\left(c\right)=\,_{3}F_{2}\left(1,1,1+\frac{c}{2};\frac{3}{2},2;1\right)=\sum_{k=0}^{\infty}\frac{k!}{k+1}\frac{\Gamma\left(\frac{c}{2}+1+k\right)}{\Gamma\left(\frac{c}{2}+1\right)}\frac{\Gamma\left(\frac{3}{2}\right)}{\Gamma\left(\frac{3}{2}+k\right)}

Implementation: `scipy.stats.t`


Triangular
==========

One shape parameter :math:`c\in[0,1]` giving the distance to the peak as a percentage of the total extent of
the non-zero portion. The location parameter is the start of the non-
zero portion, and the scale-parameter is the width of the non-zero
portion. In standard form we have :math:`x\in\left[0,1\right].`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \left\{ \begin{array}{ccc} 2\frac{x}{c} &  & x<c\\ 2\frac{1-x}{1-c} &  & x\geq c\end{array}\right.\\ F\left(x;c\right) & = & \left\{ \begin{array}{ccc} \frac{x^{2}}{c} &  & x<c\\ \frac{x^{2}-2x+c}{c-1} &  & x\geq c\end{array}\right.\\ G\left(q;c\right) & = & \left\{ \begin{array}{ccc} \sqrt{cq} &  & q<c\\ 1-\sqrt{\left(1-c\right)\left(1-q\right)} &  & q\geq c\end{array}\right.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{c}{3}+\frac{1}{3}\\ \mu_{2} & = & \frac{1-c+c^{2}}{18}\\ \gamma_{1} & = & \frac{\sqrt{2}\left(2c-1\right)\left(c+1\right)\left(c-2\right)}{5\left(1-c+c^{2}\right)^{3/2}}\\ \gamma_{2} & = & -\frac{3}{5}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left(X\right) & = & \log\left(\frac{1}{2}\sqrt{e}\right)\\  & \approx & -0.19314718055994530942.\end{eqnarray*}

Implementation: `scipy.stats.triang`


Truncated Exponential
=====================

This is an exponential distribution defined only over a certain region :math:`0<x<B` . In standard form this is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;B\right) & = & \frac{e^{-x}}{1-e^{-B}}\\ F\left(x;B\right) & = & \frac{1-e^{-x}}{1-e^{-B}}\\ G\left(q;B\right) & = & -\log\left(1-q+qe^{-B}\right)\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\Gamma\left(1+n\right)-\Gamma\left(1+n,B\right)

.. math::

     h\left[X\right]=\log\left(e^{B}-1\right)+\frac{1+e^{B}\left(B-1\right)}{1-e^{B}}.

Implementation: `scipy.stats.truncexpon`


Truncated Normal
================

A normal distribution restricted to lie within a certain range given
by two parameters :math:`A` and :math:`B` . Notice that this :math:`A` and :math:`B` correspond to the bounds on :math:`x` in standard form. For :math:`x\in\left[A,B\right]` we get

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;A,B\right) & = & \frac{\phi\left(x\right)}{\Phi\left(B\right)-\Phi\left(A\right)}\\ F\left(x;A,B\right) & = & \frac{\Phi\left(x\right)-\Phi\left(A\right)}{\Phi\left(B\right)-\Phi\left(A\right)}\\ G\left(q;A,B\right) & = & \Phi^{-1}\left[q\Phi\left(B\right)+\Phi\left(A\right)\left(1-q\right)\right]\end{eqnarray*}

where

.. math::
   :nowrap:

    \begin{eqnarray*} \phi\left(x\right) & = & \frac{1}{\sqrt{2\pi}}e^{-x^{2}/2}\\ \Phi\left(x\right) & = & \int_{-\infty}^{x}\phi\left(u\right)du.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\phi\left(A\right)-\phi\left(B\right)}{\Phi\left(B\right)-\Phi\left(A\right)}\\ \mu_{2} & = & 1+\frac{A\phi\left(A\right)-B\phi\left(B\right)}{\Phi\left(B\right)-\Phi\left(A\right)}-\left(\frac{\phi\left(A\right)-\phi\left(B\right)}{\Phi\left(B\right)-\Phi\left(A\right)}\right)^{2}\end{eqnarray*}

Implementation: `scipy.stats.truncnorm`


Tukey-Lambda
============

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\lambda\right) & = & F^{\prime}\left(x;\lambda\right)=\frac{1}{G^{\prime}\left(F\left(x;\lambda\right);\lambda\right)}=\frac{1}{F^{\lambda-1}\left(x;\lambda\right)+\left[1-F\left(x;\lambda\right)\right]^{\lambda-1}}\\ F\left(x;\lambda\right) & = & G^{-1}\left(x;\lambda\right)\\ G\left(p;\lambda\right) & = & \frac{p^{\lambda}-\left(1-p\right)^{\lambda}}{\lambda}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 0\\ \mu_{2} & = & \int_{0}^{1}G^{2}\left(p;\lambda\right)dp\\  & = & \frac{2\Gamma\left(\lambda+\frac{3}{2}\right)-\lambda4^{-\lambda}\sqrt{\pi}\Gamma\left(\lambda\right)\left(1-2\lambda\right)}{\lambda^{2}\left(1+2\lambda\right)\Gamma\left(\lambda+\frac{3}{2}\right)}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & \frac{\mu_{4}}{\mu_{2}^{2}}-3\\ \mu_{4} & = & \frac{3\Gamma\left(\lambda\right)\Gamma\left(\lambda+\frac{1}{2}\right)2^{-2\lambda}}{\lambda^{3}\Gamma\left(2\lambda+\frac{3}{2}\right)}+\frac{2}{\lambda^{4}\left(1+4\lambda\right)}\\  &  & -\frac{2\sqrt{3}\Gamma\left(\lambda\right)2^{-6\lambda}3^{3\lambda}\Gamma\left(\lambda+\frac{1}{3}\right)\Gamma\left(\lambda+\frac{2}{3}\right)}{\lambda^{3}\Gamma\left(2\lambda+\frac{3}{2}\right)\Gamma\left(\lambda+\frac{1}{2}\right)}.\end{eqnarray*}

Notice that the :math:`\lim_{\lambda\rightarrow0}G\left(p;\lambda\right)=\log\left(p/\left(1-p\right)\right)`

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \int_{0}^{1}\log\left[G^{\prime}\left(p\right)\right]dp\\  & = & \int_{0}^{1}\log\left[p^{\lambda-1}+\left(1-p\right)^{\lambda-1}\right]dp.\end{eqnarray*}

Implementation: `scipy.stats.tukeylambda`


Uniform
=======

Standard form :math:`x\in\left(0,1\right).` In general form, the lower limit is :math:`L,` the upper limit is :math:`S+L.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & 1\\ F\left(x\right) & = & x\\ G\left(q\right) & = & q\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{1}{2}\\ \mu_{2} & = & \frac{1}{12}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & -\frac{6}{5}\end{eqnarray*}

.. math::

     h\left[X\right]=0

Implementation: `scipy.stats.uniform`


Von Mises
=========

Defined for :math:`x\in\left[-\pi,\pi\right]` with shape parameter :math:`\kappa>0` . Note, the PDF and CDF functions are periodic and are always defined
over :math:`x\in\left[-\pi,\pi\right]` regardless of the location parameter. Thus, if an input beyond this
range is given, it is converted to the equivalent angle in this range.
For values of :math:`\kappa<100` the PDF and CDF formulas below are used. Otherwise, a normal
approximation with variance :math:`1/\kappa` is used.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\kappa\right) & = & \frac{e^{\kappa\cos x}}{2\pi I_{0}\left(\kappa\right)}\\ F\left(x;\kappa\right) & = & \frac{1}{2}+\frac{x}{2\pi}+\sum_{k=1}^{\infty}\frac{I_{k}\left(\kappa\right)\sin\left(kx\right)}{I_{0}\left(\kappa\right)\pi k}\\ G\left(q; \kappa\right) & = & F^{-1}\left(x;\kappa\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 0\\ \mu_{2} & = & \int_{-\pi}^{\pi}x^{2}f\left(x;\kappa\right)dx\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & \frac{\int_{-\pi}^{\pi}x^{4}f\left(x;\kappa\right)dx}{\mu_{2}^{2}}-3\end{eqnarray*}

This can be used for defining circular variance.

Implementation: `scipy.stats.vonmises`


Wald
====

Special case of the Inverse Normal with shape parameter set to :math:`1.0` . Defined for :math:`x>0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{\sqrt{2\pi x^{3}}}\exp\left(-\frac{\left(x-1\right)^{2}}{2x}\right).\\ F\left(x\right) & = & \Phi\left(\frac{x-1}{\sqrt{x}}\right)+\exp\left(2\right)\Phi\left(-\frac{x+1}{\sqrt{x}}\right)\\ G\left(q;\mu\right) & = & F^{-1}\left(q;\mu\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 1\\ \mu_{2} & = & 1\\ \gamma_{1} & = & 3\\ \gamma_{2} & = & 15\\ m_{d} & = & \frac{1}{2}\left(\sqrt{13}-3\right)\end{eqnarray*}

Implementation: `scipy.stats.wald`


Wrapped Cauchy
==============

For :math:`x\in\left[0,2\pi\right]` :math:`c\in\left(0,1\right)`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{1-c^{2}}{2\pi\left(1+c^{2}-2c\cos x\right)}\\ g_{c}\left(x\right) & = & \frac{1}{\pi}\arctan\left[\frac{1+c}{1-c}\tan\left(\frac{x}{2}\right)\right]\\ r_{c}\left(q\right) & = & 2\arctan\left[\frac{1-c}{1+c}\tan\left(\pi q\right)\right]\\ F\left(x;c\right) & = & \left\{ \begin{array}{ccc} g_{c}\left(x\right) &  & 0\leq x<\pi\\ 1-g_{c}\left(2\pi-x\right) &  & \pi\leq x\leq2\pi\end{array}\right.\\ G\left(q;c\right) & = & \left\{ \begin{array}{ccc} r_{c}\left(q\right) &  & 0\leq q<\frac{1}{2}\\ 2\pi-r_{c}\left(1-q\right) &  & \frac{1}{2}\leq q\leq1\end{array}\right.\end{eqnarray*}

.. math::

     h\left[X\right]=\log\left(2\pi\left(1-c^{2}\right)\right).

Implementation: `scipy.stats.wrapcauchy`
