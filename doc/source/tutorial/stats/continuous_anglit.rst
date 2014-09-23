
.. _continuous-anglit:

Anglit Distribution
===================

Defined over :math:`x\in\left[-\frac{\pi}{4},\frac{\pi}{4}\right]`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \sin\left(2x+\frac{\pi}{2}\right)=\cos\left(2x\right)\\ F\left(x\right) & = & \sin^{2}\left(x+\frac{\pi}{4}\right)\\ G\left(q\right) & = & \arcsin\left(\sqrt{q}\right)-\frac{\pi}{4}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 0\\ \mu_{2} & = & \frac{\pi^{2}}{16}-\frac{1}{2}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & -2\frac{\pi^{4}-96}{\left(\pi^{2}-8\right)^{2}}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & 1-\log2\\  & \approx & 0.30685281944005469058\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & \int_{-\frac{\pi}{4}}^{\frac{\pi}{4}}\cos\left(2x\right)e^{xt}dx\\  & = & \frac{4\cosh\left(\frac{\pi t}{4}\right)}{t^{2}+4}\end{eqnarray*}

.. math::

     l_{\mathbf{x}}\left(\cdot\right)=-N\overline{\log\left[\cos\left(2\mathbf{x}\right)\right]}

Implementation: `scipy.stats.anglit`
