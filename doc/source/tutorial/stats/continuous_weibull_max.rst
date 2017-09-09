
.. _continuous-weibull_max:

Weibull Maximum Extreme Value Distribution
==========================================

Defined for :math:`x<0` and :math:`c>0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & c\left(-x\right)^{c-1}\exp\left(-\left(-x\right)^{c}\right)\\ F\left(x;c\right) & = & \exp\left(-\left(-x\right)^{c}\right)\\ G\left(q;c\right) & = & -\left(-\log q\right)^{1/c}\end{eqnarray*}

The mean is the negative of the right-skewed Frechet distribution
given above, and the other statistical parameters can be computed from

.. math::

     \mu_{n}^{\prime}=\left(-1\right)^{n}\Gamma\left(1+\frac{n}{c}\right).

.. math::
   :nowrap:

    \begin{eqnarray*}
        \mu & = & -\Gamma\left(1+\frac{1}{c}\right) \\
        \mu_{2} & = & \Gamma\left(1+\frac{2}{c}\right) -
                      \Gamma^{2}\left(1+\frac{1}{c}\right) \\
        \gamma_{1} & = & -\frac{\Gamma\left(1+\frac{3}{c}\right) -
                                3\Gamma\left(1+\frac{2}{c}\right)\Gamma\left(1+\frac{1}{c}\right) +
                                2\Gamma^{3}\left(1+\frac{1}{c}\right)}
                               {\mu_{2}^{3/2}} \\
        \gamma_{2} & = & \frac{\Gamma\left(1+\frac{4}{c}\right) -
                               4\Gamma\left(1+\frac{1}{c}\right)\Gamma\left(1+\frac{3}{c}\right) +
                               6\Gamma^{2}\left(1+\frac{1}{c}\right)\Gamma\left(1+\frac{2}{c}\right) -
                               3\Gamma^{4}\left(1+\frac{1}{c}\right)}
                              {\mu_{2}^{2}} - 3 \\
        m_{d} & = & \begin{cases}
                        -\left(\frac{c-1}{c}\right)^{\frac{1}{c}} & \text{if}\; c > 1 \\
                        0 & \text{if}\; c <= 1
                    \end{cases} \\
        m_{n} & = & -\ln\left(2\right)^{\frac{1}{c}}
    \end{eqnarray*}

.. math::

     h\left[X\right]=-\frac{\gamma}{c}-\log\left(c\right)+\gamma+1

where :math:`\gamma` is Euler's constant and equal to

.. math::

     \gamma\approx0.57721566490153286061.

Implementation: `scipy.stats.weibull_max`
