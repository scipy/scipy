
.. _continuous-geninvgauss:

Generalized Inverse Gaussian Distribution
=========================================

The probability density function is given by:

.. math::
	:nowrap:
	
	\begin{eqnarray*}
	        f(x; p, b) = x^{p-1} \exp(-b(x + 1/x)/2) / (2 K_p(b)),
	\end{eqnarray*}

where :math:`x > 0` is a real number and the parameters :math:`p, b` satisfy :math:`b > 0`. :math:`K_v` is the modified Bessel function of second kind of order :math:`v` (`scipy.special.kv`).

If `X` is ``geninvgauss(p, b)``, then the distribution of `1/X` is ``geninvgauss(-p, b)``. The inverse Gaussian distribution (`scipy.stats.invgauss`) is a special case with p=-1/2.

Implementation: `scipy.stats.geninvgauss`
