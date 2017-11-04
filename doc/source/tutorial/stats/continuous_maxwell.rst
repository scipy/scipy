
.. _continuous-maxwell:

Maxwell Distribution
====================

This is a special case of the Chi distribution with :math:`L=0` and :math:`S=S=\frac{1}{\sqrt{a}}` and :math:`\nu=3.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \sqrt{\frac{2}{\pi}}x^{2}e^{-x^{2}/2}I_{\left(0,\infty\right)}\left(x\right)\\ F\left(x\right) & = & \Gamma\left(\frac{3}{2},\frac{x^{2}}{2}\right)\\ G\left(\alpha\right) & = & \sqrt{2\Gamma^{-1}\left(\frac{3}{2},\alpha\right)}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 2\sqrt{\frac{2}{\pi}}\\ \mu_{2} & = & 3-\frac{8}{\pi}\\ \gamma_{1} & = & \sqrt{2}\frac{32-10\pi}{\left(3\pi-8\right)^{3/2}}\\ \gamma_{2} & = & \frac{-12\pi^{2}+160\pi-384}{\left(3\pi-8\right)^{2}}\\ m_{d} & = & \sqrt{2}\\ m_{n} & = & \sqrt{2\Gamma^{-1}\left(\frac{3}{2},\frac{1}{2}\right)}\end{eqnarray*}

.. math::

     h\left[X\right]=\log\left(\sqrt{\frac{2\pi}{e}}\right)+\gamma.

Implementation: `scipy.stats.maxwell`
