
.. _continuous-cosine:

Cosine Distribution
===================

Approximation to the normal distribution.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{2\pi}\left[1+\cos x\right]I_{\left[-\pi,\pi\right]}\left(x\right)\\ F\left(x\right) & = & \frac{1}{2\pi}\left[\pi+x+\sin x\right]I_{\left[-\pi,\pi\right]}\left(x\right)+I_{\left(\pi,\infty\right)}\left(x\right)\\ G\left(\alpha\right) & = & F^{-1}\left(\alpha\right)\\ M\left(t\right) & = & \frac{\sinh\left(\pi t\right)}{\pi t\left(1+t^{2}\right)}\\ \mu=m_{d}=m_{n} & = & 0\\ \mu_{2} & = & \frac{\pi^{2}}{3}-2\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & \frac{-6\left(\pi^{4}-90\right)}{5\left(\pi^{2}-6\right)^{2}}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(4\pi\right)-1\\  & \approx & 1.5310242469692907930.\end{eqnarray*}

Implementation: `scipy.stats.cosine`
