
.. _continuous-maxwell:

Maxwell Distribution
====================

This is a special case of the Chi distribution with :math:`L=0` and :math:`S=\frac{1}{\sqrt{a}}` and :math:`\nu=3.`
The support is :math:`x\geq0`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \sqrt{\frac{2}{\pi}}x^{2}e^{-x^{2}/2}\\
    F\left(x\right) & = & \frac{\gamma\left(\frac{3}{2},\frac{x^2}{2}\right)}{\Gamma(\frac{3}{2})}\\
    G\left(q\right) & = & \sqrt{2\gamma^{-1}\left(\frac{3}{2},q\Gamma(\frac{3}{2})\right)}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 2\sqrt{\frac{2}{\pi}}\\
    \mu_{2} & = & 3-\frac{8}{\pi}\\
    \gamma_{1} & = & \sqrt{2}\frac{32-10\pi}{\left(3\pi-8\right)^{3/2}}\\
    \gamma_{2} & = & \frac{-12\pi^{2}+160\pi-384}{\left(3\pi-8\right)^{2}}\\
    m_{d} & = & \sqrt{2}\\
    m_{n} & = & \sqrt{2\gamma^{-1}\left(\frac{3}{2},\frac{1}{2}\Gamma(\frac{3}{2})\right)}\end{eqnarray*}

.. math::

     h\left[X\right]=\log\left(\sqrt{\frac{2\pi}{e}}\right)+\gamma.

Implementation: `scipy.stats.maxwell`
