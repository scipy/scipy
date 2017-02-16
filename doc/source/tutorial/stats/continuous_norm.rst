
.. _continuous-norm:

Normal Distribution
===================

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{e^{-x^{2}/2}}{\sqrt{2\pi}}\\ F\left(x\right) & = & \Phi\left(x\right)=\frac{1}{2}+\frac{1}{2}\mathrm{erf}\left(\frac{\mathrm{x}}{\sqrt{2}}\right)\\ G\left(q\right) & = & \Phi^{-1}\left(q\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=m_{n}=\mu & = & 0\\ \mu_{2} & = & 1\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & 0\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(\sqrt{2\pi e}\right)\\  & \approx & 1.4189385332046727418\end{eqnarray*}

Implementation: `scipy.stats.norm`
