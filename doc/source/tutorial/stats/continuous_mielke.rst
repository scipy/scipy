
.. _continuous-mielke:

Mielke's Beta-Kappa Distribution
================================

A generalized F distribution. Two shape parameters
:math:`\kappa` and :math:`\theta`, with support :math:`x\geq0`.
The :math:`\beta` in the DATAPLOT reference is a scale parameter.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\kappa,\theta\right) & = & \frac{\kappa x^{\kappa-1}}{\left(1+x^{\theta}\right)^{1+\frac{\kappa}{\theta}}}\\ F\left(x;\kappa,\theta\right) & = & \frac{x^{\kappa}}{\left(1+x^{\theta}\right)^{\kappa/\theta}}\\ G\left(q;\kappa,\theta\right) & = & \left(\frac{q^{\theta/\kappa}}{1-q^{\theta/\kappa}}\right)^{1/\theta}\end{eqnarray*}

Implementation: `scipy.stats.mielke`
