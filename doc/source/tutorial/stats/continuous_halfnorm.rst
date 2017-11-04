
.. _continuous-halfnorm:

HalfNormal Distribution
=======================

This is a special case of the chi distribution with :math:`L=a` and :math:`S=b` and :math:`\nu=1.` This is also a special case of the folded normal with shape parameter :math:`c=0` and :math:`S=S.` If :math:`Z` is (standard) normally distributed then, :math:`\left|Z\right|` is half-normal. The standard form is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \sqrt{\frac{2}{\pi}}e^{-x^{2}/2}I_{\left(0,\infty\right)}\left(x\right)\\ F\left(x\right) & = & 2\Phi\left(x\right)-1\\ G\left(q\right) & = & \Phi^{-1}\left(\frac{1+q}{2}\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\sqrt{2\pi}e^{t^{2}/2}\Phi\left(t\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \sqrt{\frac{2}{\pi}}\\ \mu_{2} & = & 1-\frac{2}{\pi}\\ \gamma_{1} & = & \frac{\sqrt{2}\left(4-\pi\right)}{\left(\pi-2\right)^{3/2}}\\ \gamma_{2} & = & \frac{8\left(\pi-3\right)}{\left(\pi-2\right)^{2}}\\ m_{d} & = & 0\\ m_{n} & = & \Phi^{-1}\left(\frac{3}{4}\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(\sqrt{\frac{\pi e}{2}}\right)\\  & \approx & 0.72579135264472743239.\end{eqnarray*}

Implementation: `scipy.stats.halfnorm`
