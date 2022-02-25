
.. _discrete-yulesimon:

Yule-Simon Distribution
========================

A Yule-Simon random variable with parameter :math:`\alpha>0`
can be represented as a mixture of
exponential random variates. To see this write :math:`W` as an exponential
random variate with rate :math:`\rho` and a Geometric random variate :math:`K`
with probability :math:`1-exp(-W)` then :math:`K` marginally has a Yule-Simon
distribution. The latent variable representation described above is used for
random variate generation.

.. math::
   :nowrap:

    \begin{eqnarray*}
    p \left( k; \alpha \right) & = & \alpha \frac{\Gamma\left(k\right)\Gamma\left(\alpha + 1\right)}{\Gamma\left(k+\alpha+1\right)} \\
    F \left( k; \alpha \right) & = &  1 - \frac{ k \Gamma\left(k\right)\Gamma\left(\alpha + 1\right)}{\Gamma\left(k+\alpha+1\right)}
    \end{eqnarray*}

for :math:`k = 1,2,...`.

Now

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\alpha}{\alpha-1}\\
    \mu_{2} & = &  \frac{\alpha^2}{\left(\alpha-1\right)^2\left( \alpha - 2 \right)}\\
    \gamma_{1} & = & \frac{ \sqrt{\left( \alpha - 2 \right)} \left( \alpha + 1 \right)^2}{ \alpha  \left( \alpha - 3 \right)}\\
    \gamma_{2} & = & \frac{ \left(\alpha + 3\right) + \left(\alpha^3 - 49\alpha - 22\right)}{\alpha \left(\alpha - 4\right)\left(\alpha - 3 \right) }
    \end{eqnarray*}

for :math:`\alpha>1` otherwise the mean is infinite and the variance does not exist.
For the variance, :math:`\alpha>2` otherwise the variance does not exist.
Similarly, for the skewness and
kurtosis to be finite, :math:`\alpha>3` and :math:`\alpha>4` respectively.


Implementation: `scipy.stats.yulesimon`
