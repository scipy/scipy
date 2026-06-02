
.. _continuous-gumbel_r:

Gumbel (LogWeibull, Fisher-Tippetts, Type I Extreme Value) Distribution
=======================================================================

One of a class of extreme value distributions (right-skewed).

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \exp\left(-\left(x+e^{-x}\right)\right)\\ F\left(x\right) & = & \exp\left(-e^{-x}\right)\\ G\left(q\right) & = & -\log\left(-\log\left(q\right)\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\Gamma\left(1-t\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \gamma=-\psi_{0}\left(1\right)\\ \mu_{2} & = & \frac{\pi^{2}}{6}\\ \gamma_{1} & = & \frac{12\sqrt{6}}{\pi^{3}}\zeta\left(3\right)\\ \gamma_{2} & = & \frac{12}{5}\\ m_{d} & = & 0\\ m_{n} & = & -\log\left(\log2\right)\end{eqnarray*}

.. math::

     h\left[X\right]\approx1.0608407169541684911

Implementation: `scipy.stats.gumbel_r`
