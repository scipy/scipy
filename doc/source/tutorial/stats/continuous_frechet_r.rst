
.. _continuous-frechet_r:

Fréchet (ExtremeLB, Extreme Value II, Weibull minimum) Distribution
====================================================================

A type of extreme-value distribution with a lower bound. Defined for :math:`x>0` and :math:`c>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & cx^{c-1}\exp\left(-x^{c}\right)\\ F\left(x;c\right) & = & 1-\exp\left(-x^{c}\right)\\ G\left(q;c\right) & = & \left[-\log\left(1-q\right)\right]^{1/c}\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\Gamma\left(1+\frac{n}{c}\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \Gamma\left(1+\frac{1}{c}\right)\\ \mu_{2} & = & \Gamma\left(1+\frac{2}{c}\right)-\Gamma^{2}\left(1-\frac{1}{c}\right)\\ \gamma_{1} & = & \frac{\Gamma\left(1+\frac{3}{c}\right)-3\Gamma\left(1+\frac{2}{c}\right)\Gamma\left(1+\frac{1}{c}\right)+2\Gamma^{3}\left(1+\frac{1}{c}\right)}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\Gamma\left(1+\frac{4}{c}\right)-4\Gamma\left(1+\frac{1}{c}\right)\Gamma\left(1+\frac{3}{c}\right)+6\Gamma^{2}\left(1+\frac{1}{c}\right)\Gamma\left(1+\frac{2}{c}\right)-\Gamma^{4}\left(1+\frac{1}{c}\right)}{\mu_{2}^{2}}-3\\ m_{d} & = & \left(\frac{c}{1+c}\right)^{1/c}\\ m_{n} & = & G\left(\frac{1}{2};c\right)\end{eqnarray*}

.. math::

     h\left[X\right]=-\frac{\gamma}{c}-\log\left(c\right)+\gamma+1

where :math:`\gamma` is Euler's constant and equal to

.. math::

     \gamma\approx0.57721566490153286061.

Implementation: `scipy.stats.frechet_r`
