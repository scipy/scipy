
.. _continuous-norminvgauss:

Normal Inverse Gaussian Distribution
==============================================

The probability density function is given by:

.. math::
	:nowrap:

	\begin{eqnarray*}
	        f(x; a, b) = \frac{a \exp\left(\sqrt{a^2 - b^2} + b x \right)}{\pi \sqrt{1 + x^2}} \, K_1\left(a * \sqrt{1 + x^2}\right),
	\end{eqnarray*}

where :math:`x` is a real number, the parameter :math:`a` is the tail heaviness and :math:`b` is the asymmetry parameter satisfying :math:`a > 0` and :math:`|b| \leq a`. :math:`K_1` is the modified Bessel function of second kind (`scipy.special.k1`).

A normal inverse Gaussian random variable with parameters :math:`a` and :math:`b` can be expressed  as :math:`X = b V + \sqrt(V) X` where :math:`X` is `norm(0,1)` and :math:`V` is `invgauss(mu=1/sqrt(a**2 - b**2))`. Hence, the normal inverse Gaussian distribution is a special case of normal variance-mean mixtures.

Another common parametrization of the distribution is given by the following expression of the pdf:

.. math::
	:nowrap:

	\begin{eqnarray*}
        g(x, \alpha, \beta, \delta, \mu) = \frac{\alpha\delta K_1 \left(\alpha\sqrt{\delta^2 + (x - \mu)^2}\right)}{\pi \sqrt{\delta^2 + (x - \mu)^2}} \,
        e^{\delta \sqrt{\alpha^2 - \beta^2} + \beta (x - \mu)}
	\end{eqnarray*}

In SciPy, this corresponds to :math:`a = \alpha \delta, b = \beta  \delta, \text{loc} = \mu, \text{scale}=\delta`.


Implementation: `scipy.stats.norminvgauss`
