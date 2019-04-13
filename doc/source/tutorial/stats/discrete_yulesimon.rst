
.. _discrete-yulesimon:

Yule-Simon Distribution
========================

A random variable has the Yule-Simon distribution
with parameter :math:`\alpha>0` if it's probability mass function is given by

.. math::
   :nowrap:

    p\left(k;\alpha\right) = \alpha \B\left(k, \alpha+1\right),

for :math:`k = 1,2,...` where

.. math::

    \B\left(a, b\right) = \frac{\Gamma\left(a\right)\Gamma\left(b\right)}{\Gamma\left(a+b\right)},

is the Beta function. The Yule-Simon can be represented as a mixture of 
exponential random variates. To see this write :math:`W` as an exponential 
random variate with rate :math:`\rho` and a Geometric random variate :math:`K` 
with probability :math:`1-exp(-W)` then :math:`K` marginally has a Yule-Simon
distribution. The latent variable representation described above is used for
random variate generation. 

The mean is 

.. math::

    \mathbb{E}(K) = \frac{\alpha}{\alpha-1},

for :math:`\alpha>1` otherwise the mean does not exist.

Similarly, the variance is 

.. math::

    \mathbb{V}ar(K) = \frac{\alpha^2}{(\alpha-1)^2(\alpha-2)}.

Implementation: `scipy.stats.yulesimon`
