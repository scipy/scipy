
.. _continuous-genhalflogistic:

Generalized Half-Logistic Distribution
======================================

For :math:`x\in\left[0,1/c\right]` and :math:`c>0` we have

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{2\left(1-cx\right)^{\frac{1}{c}-1}}{\left(1+\left(1-cx\right)^{1/c}\right)^{2}}\\ F\left(x;c\right) & = & \frac{1-\left(1-cx\right)^{1/c}}{1+\left(1-cx\right)^{1/c}}\\ G\left(q;c\right) & = & \frac{1}{c}\left[1-\left(\frac{1-q}{1+q}\right)^{c}\right]\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & 2-\left(2c+1\right)\log2.\end{eqnarray*}

Implementation: `scipy.stats.genhalflogistic`
