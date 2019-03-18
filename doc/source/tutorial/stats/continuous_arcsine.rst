
.. _continuous-arcsine:

Arcsine Distribution
====================

Defined over :math:`x\in\left[0,1\right]`.  To get the JKB definition put :math:`x=\frac{u+1}{2}.` i.e. :math:`L=-1` and :math:`S=2.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{\pi\sqrt{x\left(1-x\right)}}\\ F\left(x\right) & = & \frac{2}{\pi}\arcsin\left(\sqrt{x}\right)\\ G\left(q\right) & = & \sin^{2}\left(\frac{\pi}{2}q\right)\end{eqnarray*}

.. math::

     M\left(t\right)=E^{t/2}I_{0}\left(\frac{t}{2}\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{n}^{\prime} & = & \frac{1}{\pi}\int_{0}^{1} x^{n-1/2}\left(1-x\right)^{-1/2} dx\\
     & = & \frac{1}{\pi}B\left(\frac{1}{2},n+\frac{1}{2}\right)=\frac{\left(2n-1\right)!!}{2^{n}n!}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{1}{2}\\ \mu_{2} & = & \frac{1}{8}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & -\frac{3}{2}\end{eqnarray*}

.. math::

     h\left[X\right] = \log(\frac{\pi}{4}) \approx-0.24156447527049044468

.. math::

     l_{\mathbf{x}}\left(\cdot\right)=N\log\pi+\frac{N}{2}\overline{\log\mathbf{x}}+\frac{N}{2}\overline{\log\left(1-\mathbf{x}\right)}

Implementation: `scipy.stats.arcsine`
