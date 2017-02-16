
.. _continuous-rice:

Rice Distribution
=================

Defined for :math:`x>0` and :math:`b>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;b\right) & = & x\exp\left(-\frac{x^{2}+b^{2}}{2}\right)I_{0}\left(xb\right)\\ F\left(x;b\right) & = & \int_{0}^{x}\alpha\exp\left(-\frac{\alpha^{2}+b^{2}}{2}\right)I_{0}\left(\alpha b\right)d\alpha\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\sqrt{2^{n}}\Gamma\left(1+\frac{n}{2}\right)\,_{1}F_{1}\left(-\frac{n}{2};1;-\frac{b^{2}}{2}\right)

Implementation: `scipy.stats.rice`
