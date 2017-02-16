
.. _continuous-ncf:

Noncentral F Distribution
=========================

Let :math:`\lambda>0` and :math:`\nu_{1}>0` and :math:`\nu_{2}>0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\lambda,\nu_{1},\nu_{2}\right) & = & \exp\left[\frac{\lambda}{2}+\frac{\left(\lambda\nu_{1}x\right)}{2\left(\nu_{1}x+\nu_{2}\right)}\right]\nu_{1}^{\nu_{1}/2}\nu_{2}^{\nu_{2}/2}x^{\nu_{1}/2-1}\\  &  & \times\left(\nu_{2}+\nu_{1}x\right)^{-\left(\nu_{1}+\nu_{2}\right)/2}\frac{\Gamma\left(\frac{\nu_{1}}{2}\right)\Gamma\left(1+\frac{\nu_{2}}{2}\right)L_{\nu_{2}/2}^{\nu_{1}/2-1}\left(-\frac{\lambda\nu_{1}x}{2\left(\nu_{1}x+\nu_{2}\right)}\right)}{B\left(\frac{\nu_{1}}{2},\frac{\nu_{2}}{2}\right)\Gamma\left(\frac{\nu_{1}+\nu_{2}}{2}\right)}\end{eqnarray*}

Implementation: `scipy.stats.ncf`
