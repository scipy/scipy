
.. _continuous-tukeylambda:

Tukey-Lambda Distribution
=========================

There is one shape parameter :math:`\lambda`.  The support is :math:`x\in\mathbb{R}`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\lambda\right) & = & F^{\prime}\left(x;\lambda\right)=\frac{1}{G^{\prime}\left(F\left(x;\lambda\right);\lambda\right)}=\frac{1}{F^{\lambda-1}\left(x;\lambda\right)+\left[1-F\left(x;\lambda\right)\right]^{\lambda-1}}\\ F\left(x;\lambda\right) & = & G^{-1}\left(x;\lambda\right)\\ G\left(p;\lambda\right) & = & \frac{p^{\lambda}-\left(1-p\right)^{\lambda}}{\lambda}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 0\\ \mu_{2} & = & \int_{0}^{1}G^{2}\left(p;\lambda\right)dp\\  & = & \frac{2\Gamma\left(\lambda+\frac{3}{2}\right)-\lambda4^{-\lambda}\sqrt{\pi}\Gamma\left(\lambda\right)\left(1-2\lambda\right)}{\lambda^{2}\left(1+2\lambda\right)\Gamma\left(\lambda+\frac{3}{2}\right)}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & \frac{\mu_{4}}{\mu_{2}^{2}}-3\\ \mu_{4} & = & \frac{3\Gamma\left(\lambda\right)\Gamma\left(\lambda+\frac{1}{2}\right)2^{-2\lambda}}{\lambda^{3}\Gamma\left(2\lambda+\frac{3}{2}\right)}+\frac{2}{\lambda^{4}\left(1+4\lambda\right)}\\  &  & -\frac{2\sqrt{3}\Gamma\left(\lambda\right)2^{-6\lambda}3^{3\lambda}\Gamma\left(\lambda+\frac{1}{3}\right)\Gamma\left(\lambda+\frac{2}{3}\right)}{\lambda^{3}\Gamma\left(2\lambda+\frac{3}{2}\right)\Gamma\left(\lambda+\frac{1}{2}\right)}.\end{eqnarray*}

Notice that the :math:`\lim_{\lambda\rightarrow0}G\left(p;\lambda\right)=\log\left(p/\left(1-p\right)\right)`

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \int_{0}^{1}\log\left[G^{\prime}\left(p\right)\right]dp\\  & = & \int_{0}^{1}\log\left[p^{\lambda-1}+\left(1-p\right)^{\lambda-1}\right]dp.\end{eqnarray*}

Implementation: `scipy.stats.tukeylambda`
