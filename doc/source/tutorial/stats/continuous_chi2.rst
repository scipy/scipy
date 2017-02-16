
.. _continuous-chi2:

Chi-squared Distribution
========================

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
