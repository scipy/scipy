
.. _continuous-gilbrat:

Gilbrat Distribution
====================

Special case of the log-normal with :math:`\sigma=1` and :math:`S=1.0`, typically also :math:`L=0.0`.)

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\sigma\right) & = & \frac{1}{x\sqrt{2\pi}}\exp\left(-\frac{1}{2}\left(\log x\right)^{2}\right)\\
    F\left(x;\sigma\right) & = & \Phi\left(\log x\right)=\frac{1}{2}\left(1+\mathrm{erf}\left(\frac{\log x}{\sqrt{2}}\right)\right)\\
    G\left(q;\sigma\right) & = & \exp\left( \Phi^{-1}\left(q\right)\right) \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \sqrt{e}\\
    \mu_{2} & = & e\left[e-1\right]\\
    \gamma_{1} & = & \sqrt{e-1}\left(2+e\right)\\
    \gamma_{2} & = & e^{4}+2e^{3}+3e^{2}-6\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(\sqrt{2\pi e}\right)\\
     & \approx & 1.4189385332046727418\end{eqnarray*}

Implementation: `scipy.stats.gilbrat`
