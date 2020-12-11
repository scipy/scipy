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

    p\left(x\right) = p_{0}\left(x-L\right)

which allows for shifting of the input. When a distribution generator
is initialized, the discrete distribution can either specify the
beginning and ending (integer) values :math:`a` and :math:`b` which must be such that

.. math::

    p_{0}\left(x\right) = 0\quad x < a \textrm{ or } x > b

in which case, it is assumed that the pdf function is specified on the
integers :math:`a+mk\leq b` where :math:`k` is a non-negative integer ( :math:`0,1,2,\ldots` ) and :math:`m` is a positive integer multiplier. Alternatively, the two lists :math:`x_{k}` and :math:`p\left(x_{k}\right)` can be provided directly in which case a dictionary is set up
internally to evaluate probabilities and generate random variates.


Probability Mass Function (PMF)
-------------------------------

The probability mass function of a random variable X is defined as the
probability that the random variable takes on a particular value.

.. math::

    p\left(x_{k}\right)=P\left[X=x_{k}\right]

This is also sometimes called the probability density function,
although technically

.. math::

    f\left(x\right)=\sum_{k}p\left(x_{k}\right)\delta\left(x-x_{k}\right)

is the probability density function for a discrete distribution [#]_ .

.. [#]
    XXX: Unknown layout Plain Layout: Note that we will be using :math:`p` to represent the probability mass function and a parameter (a
    XXX: probability). The usage should be obvious from context.


Cumulative Distribution Function (CDF)
--------------------------------------

The cumulative distribution function is

.. math::

    F\left(x\right)=P\left[X\leq x\right]=\sum_{x_{k}\leq x}p\left(x_{k}\right)

and is also useful to be able to compute. Note that

.. math::

    F\left(x_{k}\right)-F\left(x_{k-1}\right)=p\left(x_{k}\right)


Survival Function
-----------------

The survival function is just

.. math::

    S\left(x\right)=1-F\left(x\right)=P\left[X>k\right]

the probability that the random variable is strictly larger than :math:`k` .

.. _discrete-ppf:

Percent Point Function (Inverse CDF)
------------------------------------

The percent point function is the inverse of the cumulative
distribution function and is

.. math::

    G\left(q\right)=F^{-1}\left(q\right)

for discrete distributions, this must be modified for cases where
there is no :math:`x_{k}` such that :math:`F\left(x_{k}\right)=q.` In these cases we choose :math:`G\left(q\right)` to be the smallest value :math:`x_{k}=G\left(q\right)` for which :math:`F\left(x_{k}\right)\geq q` . If :math:`q=0` then we define :math:`G\left(0\right)=a-1` . This definition allows random variates to be defined in the same way
as with continuous rv's using the inverse cdf on a uniform
distribution to generate random variates.


Inverse survival function
-------------------------

The inverse survival function is the inverse of the survival function

.. math::

    Z\left(\alpha\right)=S^{-1}\left(\alpha\right)=G\left(1-\alpha\right)

and is thus the smallest non-negative integer :math:`k` for which :math:`F\left(k\right)\geq1-\alpha` or the smallest non-negative integer :math:`k` for which :math:`S\left(k\right)\leq\alpha.`


Hazard functions
----------------

If desired, the hazard function and the cumulative hazard function
could be defined as

.. math::

    h\left(x_{k}\right)=\frac{p\left(x_{k}\right)}{1-F\left(x_{k}\right)}

and

.. math::

    H\left(x\right)=\sum_{x_{k}\leq x}h\left(x_{k}\right)=\sum_{x_{k}\leq x}\frac{F\left(x_{k}\right)-F\left(x_{k-1}\right)}{1-F\left(x_{k}\right)}.


Moments
-------

Non-central moments are defined using the PDF

.. math::

    \mu_{m}^{\prime}=E\left[X^{m}\right]=\sum_{k}x_{k}^{m}p\left(x_{k}\right).

Central moments are computed similarly :math:`\mu=\mu_{1}^{\prime}`

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{m}=E\left[\left(X-\mu\right)^{m}\right] & = & \sum_{k}\left(x_{k}-\mu\right)^{m}p\left(x_{k}\right)\\  & = & \sum_{k=0}^{m}\left(-1\right)^{m-k}\left(\begin{array}{c} m\\ k\end{array}\right)\mu^{m-k}\mu_{k}^{\prime}\end{eqnarray*}

The mean is the first moment

.. math::

    \mu=\mu_{1}^{\prime}=E\left[X\right]=\sum_{k}x_{k}p\left(x_{k}\right)

the variance is the second central moment

.. math::

    \mu_{2}=E\left[\left(X-\mu\right)^{2}\right]=\sum_{x_{k}}x_{k}^{2}p\left(x_{k}\right)-\mu^{2}.

Skewness is defined as

.. math::

    \gamma_{1}=\frac{\mu_{3}}{\mu_{2}^{3/2}}

while (Fisher) kurtosis is

.. math::

    \gamma_{2}=\frac{\mu_{4}}{\mu_{2}^{2}}-3,

so that a normal distribution has a kurtosis of zero.


Moment generating function
--------------------------

The moment generating function is defined as

.. math::

    M_{X}\left(t\right)=E\left[e^{Xt}\right]=\sum_{x_{k}}e^{x_{k}t}p\left(x_{k}\right)

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

    f\left(\mathbf{k};\boldsymbol{\theta}\right)=\prod_{i=1}^{N}f_{i}\left(k_{i};\boldsymbol{\theta}\right).

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

    \overline{y\left(\mathbf{x}\right)}=\frac{1}{N}\sum_{i=1}^{N}y\left(x_{i}\right)

where :math:`N` should be clear from context.


Combinations
------------

Note that

.. math::

    k!=k\cdot\left(k-1\right)\cdot\left(k-2\right)\cdot\cdots\cdot1=\Gamma\left(k+1\right)

and has special cases of

.. math::
   :nowrap:

    \begin{eqnarray*} 0! & \equiv & 1\\ k! & \equiv & 0\quad k<0\end{eqnarray*}

and

.. math::

    \left(\begin{array}{c} n\\ k\end{array}\right)=\frac{n!}{\left(n-k\right)!k!}.

If :math:`n<0` or :math:`k<0` or :math:`k>n` we define :math:`\left(\begin{array}{c} n\\ k\end{array}\right)=0`



Discrete Distributions in `scipy.stats`
---------------------------------------
.. toctree::
   :maxdepth: 1

   discrete_bernoulli
   discrete_betabinom
   discrete_binom
   discrete_boltzmann
   discrete_planck
   discrete_poisson
   discrete_geom
   discrete_nbinom
   discrete_hypergeom
   discrete_nhypergeom
   discrete_zipf
   discrete_logser
   discrete_randint
   discrete_dlaplace
   discrete_yulesimon
