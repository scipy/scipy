
.. _continuous-semicircular:

Semicircular Distribution
=========================

Defined on :math:`x\in\left[-1,1\right]`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{2}{\pi}\sqrt{1-x^{2}}\\ F\left(x\right) & = & \frac{1}{2}+\frac{1}{\pi}\left[x\sqrt{1-x^{2}}+\arcsin x\right]\\ G\left(q\right) & = & F^{-1}\left(q\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=m_{n}=\mu & = & 0\\ \mu_{2} & = & \frac{1}{4}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & -1\end{eqnarray*}

.. math::

     h\left[X\right]=0.64472988584940017414.

Implementation: `scipy.stats.semicircular`
