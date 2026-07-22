
.. _continuous-rice:

Rice Distribution
=================

There is one shape parameter :math:`b\geq0` (the "distance from the origin") and the support is :math:`x\geq0`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;b\right) & = & x\exp\left(-\frac{x^{2}+b^{2}}{2}\right)I_{0}\left(xb\right)\\
    F\left(x;b\right) & = & \int_{0}^{x}\alpha\exp\left(-\frac{\alpha^{2}+b^{2}}{2}\right)I_{0}\left(\alpha b\right)d\alpha\end{eqnarray*}

were  :math:`I_{0}(y)` is the modified Bessel function of the first kind of order 0.

.. math::

     \mu_{n}^{\prime}=\sqrt{2^{n}}\Gamma\left(1+\frac{n}{2}\right)\,_{1}F_{1}\left(-\frac{n}{2};1;-\frac{b^{2}}{2}\right)

Implementation: `scipy.stats.rice`
