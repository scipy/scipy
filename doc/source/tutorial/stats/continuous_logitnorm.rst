
.. _continuous-logitnorm:

Logit-Normal Distribution
=========================

Has two shape parameters: real :math:`\mu` and :math:`\sigma > 0`. If
:math:`X` is normally distributed with mean :math:`\mu` and standard deviation
:math:`\sigma`, then :math:`Y = \operatorname{expit}(X) = 1 / (1 + e^{-X})`
follows a logit-normal distribution. The support is :math:`0 < x < 1`.

With :math:`\operatorname{logit}(x) = \log\left(x / (1 - x)\right)`,

.. math::
   :nowrap:

    \begin{eqnarray*}
        f\left(x;\mu,\sigma\right) & = & \frac{1}{\sigma x (1-x)}
            \phi\left(\frac{\operatorname{logit}(x) - \mu}{\sigma}\right)\\
        F\left(x;\mu,\sigma\right) & = &
            \Phi\left(\frac{\operatorname{logit}(x) - \mu}{\sigma}\right)\\
        G\left(q;\mu,\sigma\right) & = &
            \operatorname{expit}\left(\mu + \sigma\Phi^{-1}\left(q\right)\right)
    \end{eqnarray*}

where :math:`\phi` and :math:`\Phi` are the PDF and CDF of the standard normal
distribution. The moments have no closed form.

When the location and scale are fixed, the maximum likelihood estimates of
:math:`\mu` and :math:`\sigma` are the sample mean and (biased) sample
standard deviation of :math:`\operatorname{logit}(x)`.

Implementation: `scipy.stats.logitnorm`
