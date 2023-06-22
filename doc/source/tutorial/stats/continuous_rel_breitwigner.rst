
.. _continuous-rel_breitwigner:

Relativistic Breit-Wigner Distribution
======================================

There is a single shape parameter :math:`\rho` which takes values in :math:`(0, \infty)`.
The support is :math:`0 \leq x < \infty`.

.. math::
   :nowrap:

   \begin{eqnarray*}
       f\left(x, \rho\right) & = & \frac{k}{\left(x^2 - \rho^2\right)^2 + \rho^2}\\
       F\left(x, \rho\right) & = & -\frac{i k\left(\frac{\tan^{-1}\left(\frac{x}{c}\right)}{c} -
                                    \frac{\tan^{-1}\left(\frac{x}{\bar{c}}\right)}{\bar{c}}\right)}{2\rho}
   \end{eqnarray*}

.. math::
   :nowrap:

   \begin{eqnarray*}
       \mu & = & \frac{k}{2\rho} \left[\frac{\pi}{2} + \tan^{-1}\left(\rho\right)\right]\\
       \mu_2 & = & \frac{k\pi}{4} \left[\frac{1 - \rho i}{\sqrt{-1 - \rho i}} + \frac{1 + \rho i}{\sqrt{-1 + \rho i}}\right]\\
       \mu_3 & = & \infty\\
       \mu_4 & = & \infty\\
   \end{eqnarray*}

where

.. math::
   :nowrap:

   \begin{eqnarray*}
   c & = & \sqrt{-\rho (\rho + i)}\\
   \bar{c} & = & \sqrt{-\rho (\rho - i)}\text{ is its complex conjugate}\\
   k & = & \frac{2\sqrt{2}\rho^2\sqrt{\rho^2 + 1}}{\pi\sqrt{\rho^2 + \rho\sqrt{\rho^2 + 1}}}
   \end{eqnarray*}

Implementation `scipy.stats.rel_breitwigner`
