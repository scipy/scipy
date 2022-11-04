
.. _continuous-recipinvgauss:

Reciprocal Inverse Gaussian Distribution
========================================

The pdf is found from the inverse gaussian (IG), :math:`f_{RIG}\left(x;\mu\right)=\frac{1}{x^{2}}f_{IG}\left(\frac{1}{x};\mu\right)` defined for :math:`x\geq0` as

.. math::
   :nowrap:

    \begin{eqnarray*} f_{IG}\left(x;\mu\right) & = & \frac{1}{\sqrt{2\pi x^{3}}}\exp\left(-\frac{\left(x-\mu\right)^{2}}{2x\mu^{2}}\right).\\ F_{IG}\left(x;\mu\right) & = & \Phi\left(\frac{1}{\sqrt{x}}\frac{x-\mu}{\mu}\right)+\exp\left(\frac{2}{\mu}\right)\Phi\left(-\frac{1}{\sqrt{x}}\frac{x+\mu}{\mu}\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} f_{RIG}\left(x;\mu\right) & = & \frac{1}{\sqrt{2\pi x}}\exp\left(-\frac{\left(1-\mu x\right)^{2}}{2x\mu^{2}}\right)\\ F_{RIG}\left(x;\mu\right) & = & 1-F_{IG}\left(\frac{1}{x},\mu\right)\\  & = & 1-\Phi\left(\frac{1}{\sqrt{x}}\frac{1-\mu x}{\mu}\right)-\exp\left(\frac{2}{\mu}\right)\Phi\left(-\frac{1}{\sqrt{x}}\frac{1+\mu x}{\mu}\right)\end{eqnarray*}

Implementation: `scipy.stats.recipinvgauss`
