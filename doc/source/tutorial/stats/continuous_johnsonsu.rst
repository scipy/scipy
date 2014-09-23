
.. _continuous-johnsonsu:

Johnson SU Distribution
=======================

Defined for all :math:`x` with two shape parameters :math:`a` and :math:`b>0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,b\right) & = & \frac{b}{\sqrt{x^{2}+1}}\phi\left(a+b\log\left(x+\sqrt{x^{2}+1}\right)\right)\\ F\left(x;a,b\right) & = & \Phi\left(a+b\log\left(x+\sqrt{x^{2}+1}\right)\right)\\ G\left(q;a,b\right) & = & \sinh\left[\frac{\Phi^{-1}\left(q\right)-a}{b}\right]\end{eqnarray*}

Implementation: `scipy.stats.johnsonsu`
