
.. _continuous-johnsonsb:

Johnson SB Distribution
=======================

There are two shape parameters :math:`a\in\mathbb{R}` and :math:`b>0`, and the support is :math:`x\in\left[0,1\right]`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,b\right) & = & \frac{b}{x\left(1-x\right)}\phi\left(a+b\log\frac{x}{1-x}\right)\\
    F\left(x;a,b\right) & = & \Phi\left(a+b\log\frac{x}{1-x}\right)\\
    G\left(q;a,b\right) & = & \frac{1}{1+\exp\left(-\frac{1}{b}\left(\Phi^{-1}\left(q\right)-a\right)\right)}\end{eqnarray*}

Implementation: `scipy.stats.johnsonsb`
