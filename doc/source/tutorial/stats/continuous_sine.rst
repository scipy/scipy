
.. _continuous-sine:

Sine Distribution
===================

Polar-angle distribution required for ensuring uniform distribution over a sphere.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{2}\sin \left(x\right) I_{\left[0,\pi\right]}\left(x\right)\\ F\left(x\right) & = & \frac{1}{2}\left[1-\cos(x)\right]I_{\left[0,\pi\right]}\left(x\right)+I_{\left(\pi,\infty\right)}\left(x\right)\\ G\left(\alpha\right) & = & 2\sin^{-1}\left(\sqrt{x}\right)\\ M\left(t\right) & = & \frac{1+e^{\pi t}}{2\left(1+t^2\right)}\\ \mu=m_{d}=m_{n} & = & \frac{\pi}{2}\\ \mu_{2} & = & \frac{\pi^2}{4}-2\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & -2\frac{\pi^4-96}{\left(\pi^2-8\right)^2}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & 1\end{eqnarray*}

Implementation: `scipy.stats.sine`
