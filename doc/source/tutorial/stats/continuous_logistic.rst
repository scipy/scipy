
.. _continuous-logistic:

Logistic (Sech-squared) Distribution
====================================

A special case of the Generalized Logistic distribution with :math:`c=1.` Defined for :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{\exp\left(-x\right)}{\left[1+\exp\left(-x\right)\right]^{2}}\\ F\left(x\right) & = & \frac{1}{1+\exp\left(-x\right)}\\ G\left(q\right) & = & -\log\left(1/q-1\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \gamma+\psi_{0}\left(1\right)=0\\ \mu_{2} & = & \frac{\pi^{2}}{6}+\psi_{1}\left(1\right)=\frac{\pi^{2}}{3}\\ \gamma_{1} & = & \frac{\psi_{2}\left(c\right)+2\zeta\left(3\right)}{\mu_{2}^{3/2}}=0\\ \gamma_{2} & = & \frac{\left(\frac{\pi^{4}}{15}+\psi_{3}\left(c\right)\right)}{\mu_{2}^{2}}=\frac{6}{5}\\ m_{d} & = & \log1=0\\ m_{n} & = & -\log\left(2-1\right)=0\end{eqnarray*}

.. math::

     h\left[X\right]=1.

Implementation: `scipy.stats.logistic`
