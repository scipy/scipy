
.. _continuous-vonmises:

Von Mises Distribution
======================

There is one shape parameter :math:`\kappa>0`, with support :math:`x\in\left[-\pi,\pi\right]`.
For values of :math:`\kappa<100` the PDF and CDF formulas below are used. Otherwise, a normal
approximation with variance :math:`1/\kappa` is used.
[Note that the PDF and CDF functions below are periodic with period :math:`2\pi`.
If an input outside :math:`x\in\left[-\pi,\pi\right]` is given, it is converted
to the equivalent angle in this range.]


.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\kappa\right) & = & \frac{e^{\kappa\cos x}}{2\pi I_{0}\left(\kappa\right)}\\
    F\left(x;\kappa\right) & = & \frac{1}{2} + \frac{x}{2\pi} + \sum_{k=1}^{\infty}\frac{I_{k}\left(\kappa\right)\sin\left(kx\right)}{I_{0}\left(\kappa\right)\pi k}\\
    G\left(q; \kappa\right) & = & F^{-1}\left(x;\kappa\right)\end{eqnarray*}


where  :math:`I_{k}(\kappa)` is a modified Bessel function of the first kind.

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 0\\
    \mu_{2} & = & \int_{-\pi}^{\pi}x^{2}f\left(x;\kappa\right)dx\\
    \gamma_{1} & = & 0\\
    \gamma_{2} & = & \frac{\int_{-\pi}^{\pi}x^{4}f\left(x;\kappa\right)dx}{\mu_{2}^{2}}-3\end{eqnarray*}

This can be used for defining circular variance.

Implementation: `scipy.stats.vonmises`
