.. _discrete-random-variables:


==================================
Discrete Statistical Distributions
==================================

Discrete random variables take on only a countable number of values.
The commonly used distributions are included in SciPy and described in
this document. Each discrete distribution can take one extra integer
parameter: :math:`L.` The relationship between the general distribution
:math:`p` and the standard distribution :math:`p_{0}` is

.. math::
   :nowrap:

    \[ p\left(x\right)=p_{0}\left(x-L\right)\]

which allows for shifting of the input. When a distribution generator
is initialized, the discrete distribution can either specify the
beginning and ending (integer) values :math:`a` and :math:`b` which must be such that

.. math::
   :nowrap:

    \[ p_{0}\left(x\right)=0\quad x<a\textrm{ or }x>b\]

in which case, it is assumed that the pdf function is specified on the
integers :math:`a+mk\leq b` where :math:`k` is a non-negative integer ( :math:`0,1,2,\ldots` ) and :math:`m` is a positive integer multiplier. Alternatively, the two lists :math:`x_{k}` and :math:`p\left(x_{k}\right)` can be provided directly in which case a dictionary is set up
internally to evaulate probabilities and generate random variates.


Probability Mass Function (PMF)
-------------------------------

The probability mass function of a random variable X is defined as the
probability that the random variable takes on a particular value.

.. math::
   :nowrap:

    \[ p\left(x_{k}\right)=P\left[X=x_{k}\right]\]

This is also sometimes called the probability density function,
although technically

.. math::
   :nowrap:

    \[ f\left(x\right)=\sum_{k}p\left(x_{k}\right)\delta\left(x-x_{k}\right)\]

is the probability density function for a discrete distribution [#]_ .



.. [#]
    XXX: Unknown layout Plain Layout: Note that we will be using :math:`p` to represent the probability mass function and a parameter (a
    XXX: probability). The usage should be obvious from context.



Cumulative Distribution Function (CDF)
--------------------------------------

The cumulative distribution function is

.. math::
   :nowrap:

    \[ F\left(x\right)=P\left[X\leq x\right]=\sum_{x_{k}\leq x}p\left(x_{k}\right)\]

and is also useful to be able to compute. Note that

.. math::
   :nowrap:

    \[ F\left(x_{k}\right)-F\left(x_{k-1}\right)=p\left(x_{k}\right)\]




Survival Function
-----------------

The survival function is just

.. math::
   :nowrap:

    \[ S\left(x\right)=1-F\left(x\right)=P\left[X>k\right]\]

the probability that the random variable is strictly larger than :math:`k` .


Percent Point Function (Inverse CDF)
------------------------------------

The percent point function is the inverse of the cumulative
distribution function and is

.. math::
   :nowrap:

    \[ G\left(q\right)=F^{-1}\left(q\right)\]

for discrete distributions, this must be modified for cases where
there is no :math:`x_{k}` such that :math:`F\left(x_{k}\right)=q.` In these cases we choose :math:`G\left(q\right)` to be the smallest value :math:`x_{k}=G\left(q\right)` for which :math:`F\left(x_{k}\right)\geq q` . If :math:`q=0` then we define :math:`G\left(0\right)=a-1` . This definition allows random variates to be defined in the same way
as with continuous rv's using the inverse cdf on a uniform
distribution to generate random variates.


Inverse survival function
-------------------------

The inverse survival function is the inverse of the survival function

.. math::
   :nowrap:

    \[ Z\left(\alpha\right)=S^{-1}\left(\alpha\right)=G\left(1-\alpha\right)\]

and is thus the smallest non-negative integer :math:`k` for which :math:`F\left(k\right)\geq1-\alpha` or the smallest non-negative integer :math:`k` for which :math:`S\left(k\right)\leq\alpha.`


Hazard functions
----------------

If desired, the hazard function and the cumulative hazard function
could be defined as

.. math::
   :nowrap:

    \[ h\left(x_{k}\right)=\frac{p\left(x_{k}\right)}{1-F\left(x_{k}\right)}\]

and

.. math::
   :nowrap:

    \[ H\left(x\right)=\sum_{x_{k}\leq x}h\left(x_{k}\right)=\sum_{x_{k}\leq x}\frac{F\left(x_{k}\right)-F\left(x_{k-1}\right)}{1-F\left(x_{k}\right)}.\]




Moments
-------

Non-central moments are defined using the PDF

.. math::
   :nowrap:

    \[ \mu_{m}^{\prime}=E\left[X^{m}\right]=\sum_{k}x_{k}^{m}p\left(x_{k}\right).\]

Central moments are computed similarly :math:`\mu=\mu_{1}^{\prime}`

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{m}=E\left[\left(X-\mu\right)^{m}\right] & = & \sum_{k}\left(x_{k}-\mu\right)^{m}p\left(x_{k}\right)\\  & = & \sum_{k=0}^{m}\left(-1\right)^{m-k}\left(\begin{array}{c} m\\ k\end{array}\right)\mu^{m-k}\mu_{k}^{\prime}\end{eqnarray*}

The mean is the first moment

.. math::
   :nowrap:

    \[ \mu=\mu_{1}^{\prime}=E\left[X\right]=\sum_{k}x_{k}p\left(x_{k}\right)\]

the variance is the second central moment

.. math::
   :nowrap:

    \[ \mu_{2}=E\left[\left(X-\mu\right)^{2}\right]=\sum_{x_{k}}x_{k}^{2}p\left(x_{k}\right)-\mu^{2}.\]

Skewness is defined as

.. math::
   :nowrap:

    \[ \gamma_{1}=\frac{\mu_{3}}{\mu_{2}^{3/2}}\]

while (Fisher) kurtosis is

.. math::
   :nowrap:

    \[ \gamma_{2}=\frac{\mu_{4}}{\mu_{2}^{2}}-3,\]

so that a normal distribution has a kurtosis of zero.


Moment generating function
--------------------------

The moment generating funtion is defined as

.. math::
   :nowrap:

    \[ M_{X}\left(t\right)=E\left[e^{Xt}\right]=\sum_{x_{k}}e^{x_{k}t}p\left(x_{k}\right)\]

Moments are found as the derivatives of the moment generating function
evaluated at :math:`0.`


Fitting data
------------

To fit data to a distribution, maximizing the likelihood function is
common. Alternatively, some distributions have well-known minimum
variance unbiased estimators. These will be chosen by default, but the
likelihood function will always be available for minimizing.

If :math:`f_{i}\left(k;\boldsymbol{\theta}\right)` is the PDF of a random-variable where :math:`\boldsymbol{\theta}` is a vector of parameters ( *e.g.* :math:`L` and :math:`S` ), then for a collection of :math:`N` independent samples from this distribution, the joint distribution the
random vector :math:`\mathbf{k}` is

.. math::
   :nowrap:

    \[ f\left(\mathbf{k};\boldsymbol{\theta}\right)=\prod_{i=1}^{N}f_{i}\left(k_{i};\boldsymbol{\theta}\right).\]

The maximum likelihood estimate of the parameters :math:`\boldsymbol{\theta}` are the parameters which maximize this function with :math:`\mathbf{x}` fixed and given by the data:

.. math::
   :nowrap:

    \begin{eqnarray*} \hat{\boldsymbol{\theta}} & = & \arg\max_{\boldsymbol{\theta}}f\left(\mathbf{k};\boldsymbol{\theta}\right)\\  & = & \arg\min_{\boldsymbol{\theta}}l_{\mathbf{k}}\left(\boldsymbol{\theta}\right).\end{eqnarray*}

Where

.. math::
   :nowrap:

    \begin{eqnarray*} l_{\mathbf{k}}\left(\boldsymbol{\theta}\right) & = & -\sum_{i=1}^{N}\log f\left(k_{i};\boldsymbol{\theta}\right)\\  & = & -N\overline{\log f\left(k_{i};\boldsymbol{\theta}\right)}\end{eqnarray*}




Standard notation for mean
--------------------------

We will use

.. math::
   :nowrap:

    \[ \overline{y\left(\mathbf{x}\right)}=\frac{1}{N}\sum_{i=1}^{N}y\left(x_{i}\right)\]

where :math:`N` should be clear from context.


Combinations
------------

Note that

.. math::
   :nowrap:

    \[ k!=k\cdot\left(k-1\right)\cdot\left(k-2\right)\cdot\cdots\cdot1=\Gamma\left(k+1\right)\]

and has special cases of

.. math::
   :nowrap:

    \begin{eqnarray*} 0! & \equiv & 1\\ k! & \equiv & 0\quad k<0\end{eqnarray*}

and

.. math::
   :nowrap:

    \[ \left(\begin{array}{c} n\\ k\end{array}\right)=\frac{n!}{\left(n-k\right)!k!}.\]

If :math:`n<0` or :math:`k<0` or :math:`k>n` we define :math:`\left(\begin{array}{c} n\\ k\end{array}\right)=0`


Bernoulli
=========

A Bernoulli random variable of parameter :math:`p` takes one of only two values :math:`X=0` or :math:`X=1` . The probability of success ( :math:`X=1` ) is :math:`p` , and the probability of failure ( :math:`X=0` ) is :math:`1-p.` It can be thought of as a binomial random variable with :math:`n=1` . The PMF is :math:`p\left(k\right)=0` for :math:`k\neq0,1` and

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;p\right) & = & \begin{cases} 1-p & k=0\\ p & k=1\end{cases}\\ F\left(x;p\right) & = & \begin{cases} 0 & x<0\\ 1-p & 0\le x<1\\ 1 & 1\leq x\end{cases}\\ G\left(q;p\right) & = & \begin{cases} 0 & 0\leq q<1-p\\ 1 & 1-p\leq q\leq1\end{cases}\\ \mu & = & p\\ \mu_{2} & = & p\left(1-p\right)\\ \gamma_{3} & = & \frac{1-2p}{\sqrt{p\left(1-p\right)}}\\ \gamma_{4} & = & \frac{1-6p\left(1-p\right)}{p\left(1-p\right)}\end{eqnarray*}





.. math::
   :nowrap:

    \[ M\left(t\right)=1-p\left(1-e^{t}\right)\]





.. math::
   :nowrap:

    \[ \mu_{m}^{\prime}=p\]





.. math::
   :nowrap:

    \[ h\left[X\right]=p\log p+\left(1-p\right)\log\left(1-p\right)\]




Binomial
========

A binomial random variable with parameters :math:`\left(n,p\right)` can be described as the sum of :math:`n` independent Bernoulli random variables of parameter :math:`p;`

.. math::
   :nowrap:

    \[ Y=\sum_{i=1}^{n}X_{i}.\]

Therefore, this random variable counts the number of successes in :math:`n` independent trials of a random experiment where the probability of
success is :math:`p.`

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;n,p\right) & = & \left(\begin{array}{c} n\\ k\end{array}\right)p^{k}\left(1-p\right)^{n-k}\,\, k\in\left\{ 0,1,\ldots n\right\} ,\\ F\left(x;n,p\right) & = & \sum_{k\leq x}\left(\begin{array}{c} n\\ k\end{array}\right)p^{k}\left(1-p\right)^{n-k}=I_{1-p}\left(n-\left\lfloor x\right\rfloor ,\left\lfloor x\right\rfloor +1\right)\quad x\geq0\end{eqnarray*}

where the incomplete beta integral is

.. math::
   :nowrap:

    \[ I_{x}\left(a,b\right)=\frac{\Gamma\left(a+b\right)}{\Gamma\left(a\right)\Gamma\left(b\right)}\int_{0}^{x}t^{a-1}\left(1-t\right)^{b-1}dt.\]

Now

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & np\\ \mu_{2} & = & np\left(1-p\right)\\ \gamma_{1} & = & \frac{1-2p}{\sqrt{np\left(1-p\right)}}\\ \gamma_{2} & = & \frac{1-6p\left(1-p\right)}{np\left(1-p\right)}.\end{eqnarray*}



.. math::
   :nowrap:

    \[ M\left(t\right)=\left[1-p\left(1-e^{t}\right)\right]^{n}\]




Boltzmann (truncated Planck)
============================



.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;N,\lambda\right) & = & \frac{1-e^{-\lambda}}{1-e^{-\lambda N}}\exp\left(-\lambda k\right)\quad k\in\left\{ 0,1,\ldots,N-1\right\} \\ F\left(x;N,\lambda\right) & = & \left\{ \begin{array}{cc} 0 & x<0\\ \frac{1-\exp\left[-\lambda\left(\left\lfloor x\right\rfloor +1\right)\right]}{1-\exp\left(-\lambda N\right)} & 0\leq x\leq N-1\\ 1 & x\geq N-1\end{array}\right.\\ G\left(q,\lambda\right) & = & \left\lceil -\frac{1}{\lambda}\log\left[1-q\left(1-e^{-\lambda N}\right)\right]-1\right\rceil \end{eqnarray*}

Define :math:`z=e^{-\lambda}`

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{z}{1-z}-\frac{Nz^{N}}{1-z^{N}}\\ \mu_{2} & = & \frac{z}{\left(1-z\right)^{2}}-\frac{N^{2}z^{N}}{\left(1-z^{N}\right)^{2}}\\ \gamma_{1} & = & \frac{z\left(1+z\right)\left(\frac{1-z^{N}}{1-z}\right)^{3}-N^{3}z^{N}\left(1+z^{N}\right)}{\left[z\left(\frac{1-z^{N}}{1-z}\right)^{2}-N^{2}z^{N}\right]^{3/2}}\\ \gamma_{2} & = & \frac{z\left(1+4z+z^{2}\right)\left(\frac{1-z^{N}}{1-z}\right)^{4}-N^{4}z^{N}\left(1+4z^{N}+z^{2N}\right)}{\left[z\left(\frac{1-z^{N}}{1-z}\right)^{2}-N^{2}z^{N}\right]^{2}}\end{eqnarray*}



.. math::
   :nowrap:

    \[ M\left(t\right)=\frac{1-e^{N\left(t-\lambda\right)}}{1-e^{t-\lambda}}\frac{1-e^{-\lambda}}{1-e^{-\lambda N}}\]




Planck (discrete exponential)
=============================

Named Planck because of its relationship to the black-body problem he
solved.



.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;\lambda\right) & = & \left(1-e^{-\lambda}\right)e^{-\lambda k}\quad k\lambda\geq0\\ F\left(x;\lambda\right) & = & 1-e^{-\lambda\left(\left\lfloor x\right\rfloor +1\right)}\quad x\lambda\geq0\\ G\left(q;\lambda\right) & = & \left\lceil -\frac{1}{\lambda}\log\left[1-q\right]-1\right\rceil .\end{eqnarray*}



.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{1}{e^{\lambda}-1}\\ \mu_{2} & = & \frac{e^{-\lambda}}{\left(1-e^{-\lambda}\right)^{2}}\\ \gamma_{1} & = & 2\cosh\left(\frac{\lambda}{2}\right)\\ \gamma_{2} & = & 4+2\cosh\left(\lambda\right)\end{eqnarray*}





.. math::
   :nowrap:

    \[ M\left(t\right)=\frac{1-e^{-\lambda}}{1-e^{t-\lambda}}\]



.. math::
   :nowrap:

    \[ h\left[X\right]=\frac{\lambda e^{-\lambda}}{1-e^{-\lambda}}-\log\left(1-e^{-\lambda}\right)\]




Poisson
=======

The Poisson random variable counts the number of successes in :math:`n` independent Bernoulli trials in the limit as :math:`n\rightarrow\infty` and :math:`p\rightarrow0` where the probability of success in each trial is :math:`p` and :math:`np=\lambda\geq0` is a constant. It can be used to approximate the Binomial random
variable or in it's own right to count the number of events that occur
in the interval :math:`\left[0,t\right]` for a process satisfying certain "sparsity "constraints. The functions are

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;\lambda\right) & = & e^{-\lambda}\frac{\lambda^{k}}{k!}\quad k\geq0,\\ F\left(x;\lambda\right) & = & \sum_{n=0}^{\left\lfloor x\right\rfloor }e^{-\lambda}\frac{\lambda^{n}}{n!}=\frac{1}{\Gamma\left(\left\lfloor x\right\rfloor +1\right)}\int_{\lambda}^{\infty}t^{\left\lfloor x\right\rfloor }e^{-t}dt,\\ \mu & = & \lambda\\ \mu_{2} & = & \lambda\\ \gamma_{1} & = & \frac{1}{\sqrt{\lambda}}\\ \gamma_{2} & = & \frac{1}{\lambda}.\end{eqnarray*}





.. math::
   :nowrap:

    \[ M\left(t\right)=\exp\left[\lambda\left(e^{t}-1\right)\right].\]




Geometric
=========

The geometric random variable with parameter :math:`p\in\left(0,1\right)` can be defined as the number of trials required to obtain a success
where the probability of success on each trial is :math:`p` . Thus,

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;p\right) & = & \left(1-p\right)^{k-1}p\quad k\geq1\\ F\left(x;p\right) & = & 1-\left(1-p\right)^{\left\lfloor x\right\rfloor }\quad x\geq1\\ G\left(q;p\right) & = & \left\lceil \frac{\log\left(1-q\right)}{\log\left(1-p\right)}\right\rceil \\ \mu & = & \frac{1}{p}\\ \mu_{2} & = & \frac{1-p}{p^{2}}\\ \gamma_{1} & = & \frac{2-p}{\sqrt{1-p}}\\ \gamma_{2} & = & \frac{p^{2}-6p+6}{1-p}.\end{eqnarray*}





.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & \frac{p}{e^{-t}-\left(1-p\right)}\end{eqnarray*}




Negative Binomial
=================

The negative binomial random variable with parameters :math:`n` and :math:`p\in\left(0,1\right)` can be defined as the number of *extra* independent trials (beyond :math:`n` ) required to accumulate a total of :math:`n` successes where the probability of a success on each trial is :math:`p.` Equivalently, this random variable is the number of failures
encoutered while accumulating :math:`n` successes during independent trials of an experiment that succeeds
with probability :math:`p.` Thus,

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;n,p\right) & = & \left(\begin{array}{c} k+n-1\\ n-1\end{array}\right)p^{n}\left(1-p\right)^{k}\quad k\geq0\\ F\left(x;n,p\right) & = & \sum_{i=0}^{\left\lfloor x\right\rfloor }\left(\begin{array}{c} i+n-1\\ i\end{array}\right)p^{n}\left(1-p\right)^{i}\quad x\geq0\\  & = & I_{p}\left(n,\left\lfloor x\right\rfloor +1\right)\quad x\geq0\\ \mu & = & n\frac{1-p}{p}\\ \mu_{2} & = & n\frac{1-p}{p^{2}}\\ \gamma_{1} & = & \frac{2-p}{\sqrt{n\left(1-p\right)}}\\ \gamma_{2} & = & \frac{p^{2}+6\left(1-p\right)}{n\left(1-p\right)}.\end{eqnarray*}

Recall that :math:`I_{p}\left(a,b\right)` is the incomplete beta integral.


Hypergeometric
==============

The hypergeometric random variable with parameters :math:`\left(M,n,N\right)` counts the number of "good "objects in a sample of size :math:`N` chosen without replacement from a population of :math:`M` objects where :math:`n` is the number of "good "objects in the total population.

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;N,n,M\right) & = & \frac{\left(\begin{array}{c} n\\ k\end{array}\right)\left(\begin{array}{c} M-n\\ N-k\end{array}\right)}{\left(\begin{array}{c} M\\ N\end{array}\right)}\quad N-\left(M-n\right)\leq k\leq\min\left(n,N\right)\\ F\left(x;N,n,M\right) & = & \sum_{k=0}^{\left\lfloor x\right\rfloor }\frac{\left(\begin{array}{c} m\\ k\end{array}\right)\left(\begin{array}{c} N-m\\ n-k\end{array}\right)}{\left(\begin{array}{c} N\\ n\end{array}\right)},\\ \mu & = & \frac{nN}{M}\\ \mu_{2} & = & \frac{nN\left(M-n\right)\left(M-N\right)}{M^{2}\left(M-1\right)}\\ \gamma_{1} & = & \frac{\left(M-2n\right)\left(M-2N\right)}{M-2}\sqrt{\frac{M-1}{nN\left(M-m\right)\left(M-n\right)}}\\ \gamma_{2} & = & \frac{g\left(N,n,M\right)}{nN\left(M-n\right)\left(M-3\right)\left(M-2\right)\left(N-M\right)}\end{eqnarray*}

where (defining :math:`m=M-n` )

.. math::
   :nowrap:

    \begin{eqnarray*} g\left(N,n,M\right) & = & m^{3}-m^{5}+3m^{2}n-6m^{3}n+m^{4}n+3mn^{2}\\  &  & -12m^{2}n^{2}+8m^{3}n^{2}+n^{3}-6mn^{3}+8m^{2}n^{3}\\  &  & +mn^{4}-n^{5}-6m^{3}N+6m^{4}N+18m^{2}nN\\  &  & -6m^{3}nN+18mn^{2}N-24m^{2}n^{2}N-6n^{3}N\\  &  & -6mn^{3}N+6n^{4}N+6m^{2}N^{2}-6m^{3}N^{2}-24mnN^{2}\\  &  & +12m^{2}nN^{2}+6n^{2}N^{2}+12mn^{2}N^{2}-6n^{3}N^{2}.\end{eqnarray*}




Zipf (Zeta)
===========

A random variable has the zeta distribution (also called the zipf
distribution) with parameter :math:`\alpha>1` if it's probability mass function is given by

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;\alpha\right) & = & \frac{1}{\zeta\left(\alpha\right)k^{\alpha}}\quad k\geq1\end{eqnarray*}

where

.. math::
   :nowrap:

    \[ \zeta\left(\alpha\right)=\sum_{n=1}^{\infty}\frac{1}{n^{\alpha}}\]

is the Riemann zeta function. Other functions of this distribution are

.. math::
   :nowrap:

    \begin{eqnarray*} F\left(x;\alpha\right) & = & \frac{1}{\zeta\left(\alpha\right)}\sum_{k=1}^{\left\lfloor x\right\rfloor }\frac{1}{k^{\alpha}}\\ \mu & = & \frac{\zeta_{1}}{\zeta_{0}}\quad\alpha>2\\ \mu_{2} & = & \frac{\zeta_{2}\zeta_{0}-\zeta_{1}^{2}}{\zeta_{0}^{2}}\quad\alpha>3\\ \gamma_{1} & = & \frac{\zeta_{3}\zeta_{0}^{2}-3\zeta_{0}\zeta_{1}\zeta_{2}+2\zeta_{1}^{3}}{\left[\zeta_{2}\zeta_{0}-\zeta_{1}^{2}\right]^{3/2}}\quad\alpha>4\\ \gamma_{2} & = & \frac{\zeta_{4}\zeta_{0}^{3}-4\zeta_{3}\zeta_{1}\zeta_{0}^{2}+12\zeta_{2}\zeta_{1}^{2}\zeta_{0}-6\zeta_{1}^{4}-3\zeta_{2}^{2}\zeta_{0}^{2}}{\left(\zeta_{2}\zeta_{0}-\zeta_{1}^{2}\right)^{2}}.\end{eqnarray*}





.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & \frac{\textrm{Li}_{\alpha}\left(e^{t}\right)}{\zeta\left(\alpha\right)}\end{eqnarray*}

where :math:`\zeta_{i}=\zeta\left(\alpha-i\right)` and :math:`\textrm{Li}_{n}\left(z\right)` is the :math:`n^{\textrm{th}}` polylogarithm function of :math:`z` defined as

.. math::
   :nowrap:

    \[ \textrm{Li}_{n}\left(z\right)\equiv\sum_{k=1}^{\infty}\frac{z^{k}}{k^{n}}\]



.. math::
   :nowrap:

    \[ \mu_{n}^{\prime}=\left.M^{\left(n\right)}\left(t\right)\right|_{t=0}=\left.\frac{\textrm{Li}_{\alpha-n}\left(e^{t}\right)}{\zeta\left(a\right)}\right|_{t=0}=\frac{\zeta\left(\alpha-n\right)}{\zeta\left(\alpha\right)}\]




Logarithmic (Log-Series, Series)
================================

The logarimthic distribution with parameter :math:`p` has a probability mass function with terms proportional to the Taylor
series expansion of :math:`\log\left(1-p\right)`

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;p\right) & = & -\frac{p^{k}}{k\log\left(1-p\right)}\quad k\geq1\\ F\left(x;p\right) & = & -\frac{1}{\log\left(1-p\right)}\sum_{k=1}^{\left\lfloor x\right\rfloor }\frac{p^{k}}{k}=1+\frac{p^{1+\left\lfloor x\right\rfloor }\Phi\left(p,1,1+\left\lfloor x\right\rfloor \right)}{\log\left(1-p\right)}\end{eqnarray*}

where

.. math::
   :nowrap:

    \[ \Phi\left(z,s,a\right)=\sum_{k=0}^{\infty}\frac{z^{k}}{\left(a+k\right)^{s}}\]

is the Lerch Transcendent. Also define :math:`r=\log\left(1-p\right)`

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & -\frac{p}{\left(1-p\right)r}\\ \mu_{2} & = & -\frac{p\left[p+r\right]}{\left(1-p\right)^{2}r^{2}}\\ \gamma_{1} & = & -\frac{2p^{2}+3pr+\left(1+p\right)r^{2}}{r\left(p+r\right)\sqrt{-p\left(p+r\right)}}r\\ \gamma_{2} & = & -\frac{6p^{3}+12p^{2}r+p\left(4p+7\right)r^{2}+\left(p^{2}+4p+1\right)r^{3}}{p\left(p+r\right)^{2}}.\end{eqnarray*}



.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & -\frac{1}{\log\left(1-p\right)}\sum_{k=1}^{\infty}\frac{e^{tk}p^{k}}{k}\\  & = & \frac{\log\left(1-pe^{t}\right)}{\log\left(1-p\right)}\end{eqnarray*}

Thus,

.. math::
   :nowrap:

    \[ \mu_{n}^{\prime}=\left.M^{\left(n\right)}\left(t\right)\right|_{t=0}=\left.\frac{\textrm{Li}_{1-n}\left(pe^{t}\right)}{\log\left(1-p\right)}\right|_{t=0}=-\frac{\textrm{Li}_{1-n}\left(p\right)}{\log\left(1-p\right)}.\]




Discrete Uniform (randint)
==========================

The discrete uniform distribution with parameters :math:`\left(a,b\right)` constructs a random variable that has an equal probability of being
any one of the integers in the half-open range :math:`[a,b).` If :math:`a` is not given it is assumed to be zero and the only parameter is :math:`b.` Therefore,

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;a,b\right) & = & \frac{1}{b-a}\quad a\leq k<b\\ F\left(x;a,b\right) & = & \frac{\left\lfloor x\right\rfloor -a}{b-a}\quad a\leq x\leq b\\ G\left(q;a,b\right) & = & \left\lceil q\left(b-a\right)+a\right\rceil \\ \mu & = & \frac{b+a-1}{2}\\ \mu_{2} & = & \frac{\left(b-a-1\right)\left(b-a+1\right)}{12}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & -\frac{6}{5}\frac{\left(b-a\right)^{2}+1}{\left(b-a-1\right)\left(b-a+1\right)}.\end{eqnarray*}





.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & \frac{1}{b-a}\sum_{k=a}^{b-1}e^{tk}\\  & = & \frac{e^{bt}-e^{at}}{\left(b-a\right)\left(e^{t}-1\right)}\end{eqnarray*}




Discrete Laplacian
==================

Defined over all integers for :math:`a>0`

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k\right) & = & \tanh\left(\frac{a}{2}\right)e^{-a\left|k\right|},\\ F\left(x\right) & = & \left\{ \begin{array}{cc} \frac{e^{a\left(\left\lfloor x\right\rfloor +1\right)}}{e^{a}+1} & \left\lfloor x\right\rfloor <0,\\ 1-\frac{e^{-a\left\lfloor x\right\rfloor }}{e^{a}+1} & \left\lfloor x\right\rfloor \geq0.\end{array}\right.\\ G\left(q\right) & = & \left\{ \begin{array}{cc} \left\lceil \frac{1}{a}\log\left[q\left(e^{a}+1\right)\right]-1\right\rceil  & q<\frac{1}{1+e^{-a}},\\ \left\lceil -\frac{1}{a}\log\left[\left(1-q\right)\left(1+e^{a}\right)\right]\right\rceil  & q\geq\frac{1}{1+e^{-a}}.\end{array}\right.\end{eqnarray*}



.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & \tanh\left(\frac{a}{2}\right)\sum_{k=-\infty}^{\infty}e^{tk}e^{-a\left|k\right|}\\  & = & C\left(1+\sum_{k=1}^{\infty}e^{-\left(t+a\right)k}+\sum_{1}^{\infty}e^{\left(t-a\right)k}\right)\\  & = & \tanh\left(\frac{a}{2}\right)\left(1+\frac{e^{-\left(t+a\right)}}{1-e^{-\left(t+a\right)}}+\frac{e^{t-a}}{1-e^{t-a}}\right)\\  & = & \frac{\tanh\left(\frac{a}{2}\right)\sinh a}{\cosh a-\cosh t}.\end{eqnarray*}

Thus,

.. math::
   :nowrap:

    \[ \mu_{n}^{\prime}=M^{\left(n\right)}\left(0\right)=\left[1+\left(-1\right)^{n}\right]\textrm{Li}_{-n}\left(e^{-a}\right)\]

where :math:`\textrm{Li}_{-n}\left(z\right)` is the polylogarithm function of order :math:`-n` evaluated at :math:`z.`

.. math::
   :nowrap:

    \[ h\left[X\right]=-\log\left(\tanh\left(\frac{a}{2}\right)\right)+\frac{a}{\sinh a}\]




Discrete Gaussian*
==================

Defined for all :math:`\mu` and :math:`\lambda>0` and :math:`k`

.. math::
   :nowrap:

    \[ p\left(k;\mu,\lambda\right)=\frac{1}{Z\left(\lambda\right)}\exp\left[-\lambda\left(k-\mu\right)^{2}\right]\]

where

.. math::
   :nowrap:

    \[ Z\left(\lambda\right)=\sum_{k=-\infty}^{\infty}\exp\left[-\lambda k^{2}\right]\]



.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu\\ \mu_{2} & = & -\frac{\partial}{\partial\lambda}\log Z\left(\lambda\right)\\  & = & G\left(\lambda\right)e^{-\lambda}\end{eqnarray*}

where :math:`G\left(0\right)\rightarrow\infty` and :math:`G\left(\infty\right)\rightarrow2` with a minimum less than 2 near :math:`\lambda=1`

.. math::
   :nowrap:

    \[ G\left(\lambda\right)=\frac{1}{Z\left(\lambda\right)}\sum_{k=-\infty}^{\infty}k^{2}\exp\left[-\lambda\left(k+1\right)\left(k-1\right)\right]\]
