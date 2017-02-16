
.. _continuous-halfcauchy:

HalfCauchy Distribution
=======================

If :math:`Z` is Hyperbolic Secant distributed then :math:`e^{Z}` is Half-Cauchy distributed. Also, if :math:`W` is (standard) Cauchy distributed, then :math:`\left|W\right|` is Half-Cauchy distributed. Special case of the Folded Cauchy
distribution with :math:`c=0.` The standard form is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{2}{\pi\left(1+x^{2}\right)}I_{[0,\infty)}\left(x\right)\\ F\left(x\right) & = & \frac{2}{\pi}\arctan\left(x\right)I_{\left[0,\infty\right]}\left(x\right)\\ G\left(q\right) & = & \tan\left(\frac{\pi}{2}q\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\cos t+\frac{2}{\pi}\left[\mathrm{Si}\left(t\right)\cos t-\mathrm{Ci}\left(\mathrm{-}t\right)\sin t\right]

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d} & = & 0\\ m_{n} & = & \tan\left(\frac{\pi}{4}\right)\end{eqnarray*}

No moments, as the integrals diverge.

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(2\pi\right)\\  & \approx & 1.8378770664093454836.\end{eqnarray*}

Implementation: `scipy.stats.halfcauchy`
