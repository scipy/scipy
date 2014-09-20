
.. _continuous-genextreme:

Generalized Extreme Value Distribution
======================================

Extreme value distributions with shape parameter :math:`c` .

For :math:`c>0` defined on :math:`-\infty<x\leq1/c.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \exp\left[-\left(1-cx\right)^{1/c}\right]\left(1-cx\right)^{1/c-1}\\ F\left(x;c\right) & = & \exp\left[-\left(1-cx\right)^{1/c}\right]\\ G\left(q;c\right) & = & \frac{1}{c}\left[1-\left(-\log q\right)^{c}\right]\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\frac{1}{c^{n}}\sum_{k=0}^{n}\left(\begin{array}{c} n\\ k\end{array}\right)\left(-1\right)^{k}\Gamma\left(ck+1\right)\quad cn>-1

So,

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{1}^{\prime} & = & \frac{1}{c}\left(1-\Gamma\left(1+c\right)\right)\quad c>-1\\ \mu_{2}^{\prime} & = & \frac{1}{c^{2}}\left(1-2\Gamma\left(1+c\right)+\Gamma\left(1+2c\right)\right)\quad c>-\frac{1}{2}\\ \mu_{3}^{\prime} & = & \frac{1}{c^{3}}\left(1-3\Gamma\left(1+c\right)+3\Gamma\left(1+2c\right)-\Gamma\left(1+3c\right)\right)\quad c>-\frac{1}{3}\\ \mu_{4}^{\prime} & = & \frac{1}{c^{4}}\left(1-4\Gamma\left(1+c\right)+6\Gamma\left(1+2c\right)-4\Gamma\left(1+3c\right)+\Gamma\left(1+4c\right)\right)\quad c>-\frac{1}{4}\end{eqnarray*}

For :math:`c<0` defined on :math:`\frac{1}{c}\leq x<\infty.` For :math:`c=0` defined over all space

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;0\right) & = & \exp\left[-e^{-x}\right]e^{-x}\\ F\left(x;0\right) & = & \exp\left[-e^{-x}\right]\\ G\left(q;0\right) & = & -\log\left(-\log q\right)\end{eqnarray*}

This is just the (left-skewed) Gumbel distribution for c=0.

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \gamma=-\psi_{0}\left(1\right)\\ \mu_{2} & = & \frac{\pi^{2}}{6}\\ \gamma_{1} & = & \frac{12\sqrt{6}}{\pi^{3}}\zeta\left(3\right)\\ \gamma_{2} & = & \frac{12}{5}\end{eqnarray*}

Implementation: `scipy.stats.genextreme`
