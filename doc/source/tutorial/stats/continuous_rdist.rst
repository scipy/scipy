
.. _continuous-rdist:

R-distribution Distribution
===========================

A general-purpose distribution with a variety of shapes controlled by one shape parameter :math:`c>0.`
The support of the standard distribution is :math:`x\in\left[-1,1\right]`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{\left(1-x^{2}\right)^{c/2-1}}{B\left(\frac{1}{2},\frac{c}{2}\right)}\\ F\left(x;c\right) & = & \frac{1}{2}+\frac{x}{B\left(\frac{1}{2},\frac{c}{2}\right)}\,_{2}F_{1}\left(\frac{1}{2},1-\frac{c}{2};\frac{3}{2};x^{2}\right)\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\frac{\left(1+\left(-1\right)^{n}\right)}{2}B\left(\frac{n+1}{2},\frac{c}{2}\right)

The R-distribution with parameter :math:`n` is the distribution of the
correlation coefficient of a random sample of size :math:`n` drawn from a
bivariate normal distribution with :math:`\rho=0.` The mean of the standard
distribution is always zero and as the sample size grows, the distribution's
mass concentrates more closely about this mean.

Implementation: `scipy.stats.rdist`
