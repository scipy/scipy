
.. _continuous-ncf:

Noncentral F Distribution
=========================

The distribution of :math:`\left(X_{1}/X_{2}\right)\left(\nu_{2}/\nu_{1}\right)`
if :math:`X_{1}` is non-central chi-squared with :math:`v_{1}` degrees of freedom
and parameter :math:`\lambda`, and :math:`X_{2}` is chi-squared with :math:`v_{2}` degrees of freedom.

There are 3 shape parameters: the degrees of freedom :math:`\nu_{1}>0` and :math:`\nu_{2}>0`; and :math:`\lambda>0`.


.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\lambda,\nu_{1},\nu_{2}\right) & = & \exp\left[\frac{\lambda}{2}+\frac{\left(\lambda\nu_{1}x\right)}{2\left(\nu_{1}x+\nu_{2}\right)}\right]\nu_{1}^{\nu_{1}/2}\nu_{2}^{\nu_{2}/2}x^{\nu_{1}/2-1}\\
    &  & \times\left(\nu_{2}+\nu_{1}x\right)^{-\left(\nu_{1}+\nu_{2}\right)/2}\frac{\Gamma\left(\frac{\nu_{1}}{2}\right)\Gamma\left(1+\frac{\nu_{2}}{2}\right)L_{\nu_{2}/2}^{\nu_{1}/2-1}\left(-\frac{\lambda\nu_{1}x}{2\left(\nu_{1}x+\nu_{2}\right)}\right)}{B\left(\frac{\nu_{1}}{2},\frac{\nu_{2}}{2}\right)\Gamma\left(\frac{\nu_{1}+\nu_{2}}{2}\right)}\end{eqnarray*}

where :math:`L_{\nu_{2}/2}^{\nu_{1}/2-1}(x)` is an associated Laguerre polynomial.

Implementation: `scipy.stats.ncf`
