
.. _continuous-wald:

Wald Distribution
=================

Special case of the Inverse Normal with shape parameter set to :math:`1.0`. It has support :math:`x\geq0`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{\sqrt{2\pi x^{3}}}\exp\left(-\frac{\left(x-1\right)^{2}}{2x}\right).\\ F\left(x\right) & = & \Phi\left(\frac{x-1}{\sqrt{x}}\right)+\exp\left(2\right)\Phi\left(-\frac{x+1}{\sqrt{x}}\right)\\ G\left(q;\mu\right) & = & F^{-1}\left(q;\mu\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 1\\ \mu_{2} & = & 1\\ \gamma_{1} & = & 3\\ \gamma_{2} & = & 15\\ m_{d} & = & \frac{1}{2}\left(\sqrt{13}-3\right)\end{eqnarray*}

Implementation: `scipy.stats.wald`
