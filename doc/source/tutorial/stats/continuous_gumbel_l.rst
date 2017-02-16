
.. _continuous-gumbel_l:

Gumbel Left-skewed (for minimum order statistic) Distribution
=============================================================

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \exp\left(x-e^{x}\right)\\ F\left(x\right) & = & 1-\exp\left(-e^{x}\right)\\ G\left(q\right) & = & \log\left(-\log\left(1-q\right)\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\Gamma\left(1+t\right)

Note, that :math:`\mu` is negative the mean for the right-skewed distribution. Similar for
median and mode. All other moments are the same.

.. math::

     h\left[X\right]\approx1.0608407169541684911.

Implementation: `scipy.stats.gumbel_l`
