
.. _continuous-levy_l:

Left-skewed Lévy Distribution
==============================

Special case of Lévy-stable distribution with :math:`\alpha=\frac{1}{2}` and :math:`\beta=-1` the support is :math:`x<0` . In standard form

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{\left|x\right|\sqrt{2\pi\left|x\right|}}\exp\left(-\frac{1}{2\left|x\right|}\right)\\ F\left(x\right) & = & 2\Phi\left(\frac{1}{\sqrt{\left|x\right|}}\right)-1\\ G\left(q\right) & = & -\left[\Phi^{-1}\left(\frac{q+1}{2}\right)\right]^{-2}.\end{eqnarray*}

No moments.

Implementation: `scipy.stats.levy_l`
