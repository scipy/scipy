
.. _continuous-ncx2:

Noncentral chi-squared Distribution
===================================

The distribution of :math:`\sum_{i=1}^{\nu}\left(Z_{i}+\delta_{i}\right)^{2}` where :math:`Z_{i}` are independent standard normal variables and :math:`\delta_{i}` are constants. :math:`\lambda=\sum_{i=1}^{\nu}\delta_{i}^{2}>0.` (In communications it is called the Marcum-Q function). Can be thought
of as a Generalized Rayleigh-Rice distribution. For :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\nu,\lambda\right) & = & e^{-\left(\lambda+x\right)/2}\frac{1}{2}\left(\frac{x}{\lambda}\right)^{\left(\nu-2\right)/4}I_{\left(\nu-2\right)/2}\left(\sqrt{\lambda x}\right)\\ F\left(x;\nu,\lambda\right) & = & \sum_{j=0}^{\infty}\left\{ \frac{\left(\lambda/2\right)^{j}}{j!}e^{-\lambda/2}\right\} \mathrm{Pr}\left[\chi_{\nu+2j}^{2}\leq x\right]\\ G\left(q;\nu,\lambda\right) & = & F^{-1}\left(x;\nu,\lambda\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \nu+\lambda\\ \mu_{2} & = & 2\left(\nu+2\lambda\right)\\ \gamma_{1} & = & \frac{\sqrt{8}\left(\nu+3\lambda\right)}{\left(\nu+2\lambda\right)^{3/2}}\\ \gamma_{2} & = & \frac{12\left(\nu+4\lambda\right)}{\left(\nu+2\lambda\right)^{2}}\end{eqnarray*}

Implementation: `scipy.stats.ncx2`
