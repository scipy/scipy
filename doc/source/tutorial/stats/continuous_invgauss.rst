
.. _continuous-invgauss:

Inverse Normal (Inverse Gaussian) Distribution
==============================================

The standard form involves the shape parameter :math:`\mu` (in most
definitions, :math:`L=0.0` is used). (In terms of the regress
documentation :math:`\mu=A/B` ) and :math:`B=S` and :math:`L` is not
a parameter in that distribution. A standard form is :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\mu\right) & = & \frac{1}{\sqrt{2\pi x^{3}}}\exp\left(-\frac{\left(x-\mu\right)^{2}}{2x\mu^{2}}\right).\\
    F\left(x;\mu\right) & = & \Phi\left(\frac{1}{\sqrt{x}}\frac{x-\mu}{\mu}\right)+\exp\left(\frac{2}{\mu}\right)\Phi\left(-\frac{1}{\sqrt{x}}\frac{x+\mu}{\mu}\right)\\
    G\left(q;\mu\right) & = & F^{-1}\left(q;\mu\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu\\
    \mu_{2} & = & \mu^{3}\\
    \gamma_{1} & = & 3\sqrt{\mu}\\
    \gamma_{2} & = & 15\mu\\
    m_{d} & = & \frac{\mu}{2}\left(\sqrt{9\mu^{2}+4}-3\mu\right)\end{eqnarray*}

This is related to the canonical form or JKB "two-parameter" inverse
Gaussian when written in it's full form with scale parameter
:math:`S` and location parameter :math:`L` by taking
:math:`L=0` and :math:`S\equiv\lambda,` then :math:`\mu S` is equal to
:math:`\mu_{2}` where :math:`\mu_{2}` is the parameter used by JKB.
We prefer this form because of it's consistent use of the scale parameter.
Notice that in JKB the skew :math:`\left(\sqrt{\beta_{1}}\right)` and the
kurtosis ( :math:`\beta_{2}-3` ) are both functions only of
:math:`\mu_{2}/\lambda=\mu S/S=\mu` as shown here, while the variance
and mean of the standard form here are transformed appropriately.

Implementation: `scipy.stats.invgauss`
