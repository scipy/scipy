
.. _continuous-rayleigh:

Rayleigh Distribution
=====================

This is Chi distribution with :math:`L=0.0` and :math:`\nu=2` and :math:`S=S` (no location parameter is generally used), the mode of the
distribution is :math:`S.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(r\right) & = & re^{-r^{2}/2}I_{[0,\infty)}\left(x\right)\\ F\left(r\right) & = & 1-e^{-r^{2}/2}I_{[0,\infty)}\left(x\right)\\ G\left(q\right) & = & \sqrt{-2\log\left(1-q\right)}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \sqrt{\frac{\pi}{2}}\\ \mu_{2} & = & \frac{4-\pi}{2}\\ \gamma_{1} & = & \frac{2\left(\pi-3\right)\sqrt{\pi}}{\left(4-\pi\right)^{3/2}}\\ \gamma_{2} & = & \frac{24\pi-6\pi^{2}-16}{\left(4-\pi\right)^{2}}\\ m_{d} & = & 1\\ m_{n} & = & \sqrt{2\log\left(2\right)}\end{eqnarray*}

.. math::

     h\left[X\right]=\frac{\gamma}{2}+\log\left(\frac{e}{\sqrt{2}}\right).

.. math::

     \mu_{n}^{\prime}=\sqrt{2^{n}}\Gamma\left(\frac{n}{2}+1\right)

Implementation: `scipy.stats.rayleigh`
