
.. _continuous-cauchy:

Cauchy Distribution
===================

The support is :math:`x\in\mathbb{R}`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{\pi\left(1+x^{2}\right)}\\
    F\left(x\right) & = & \frac{1}{2}+\frac{1}{\pi}\tan^{-1}x\\
    G\left(q\right) & = & \tan\left(\pi q-\frac{\pi}{2}\right)\\
    m_{d} & = & 0\\
    m_{n} & = & 0\end{eqnarray*}

No finite moments. This is the :math:`t` distribution with one degree of
freedom.

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(4\pi\right)\\  & \approx & 2.5310242469692907930.\end{eqnarray*}

Implementation: `scipy.stats.cauchy`
