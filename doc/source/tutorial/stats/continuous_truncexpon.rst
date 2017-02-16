
.. _continuous-truncexpon:

Truncated Exponential Distribution
==================================

This is an exponential distribution defined only over a certain region :math:`0<x<B` . In standard form this is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;B\right) & = & \frac{e^{-x}}{1-e^{-B}}\\ F\left(x;B\right) & = & \frac{1-e^{-x}}{1-e^{-B}}\\ G\left(q;B\right) & = & -\log\left(1-q+qe^{-B}\right)\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\Gamma\left(1+n\right)-\Gamma\left(1+n,B\right)

.. math::

     h\left[X\right]=\log\left(e^{B}-1\right)+\frac{1+e^{B}\left(B-1\right)}{1-e^{B}}.

Implementation: `scipy.stats.truncexpon`
