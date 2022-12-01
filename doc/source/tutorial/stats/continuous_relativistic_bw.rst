
.. _continuous-relativistic_bw:

Relativistic Breit-Wigner Distribution
======================================

There is a single shape parameter :math:`\rho`. The support is
:math:`0 \leq x < \infty`

.. math::
   :nowrap:


   \begin{eqnarray*}
       f\left(x, \rho\right) & = & \frac{k}{\left(x^2 - \rho^2\right)^2 + \rho^2}\\
       F\left(x, \rho\right) & = & -\frac{i k\left(\frac{\tan^{-1}\left(\frac{x}{c}\right)}{c} -
                                    \frac{\tan^{-1}\left(\frac{x}{\bar{c}}\right)}{\bar{c}}\right)}{2\rho}
   \end{eqnarray*}

where

.. math::
   :nowrap:


   \begin{eqnarray*}
   c & = & \sqrt{-\rho (\rho + i)}\\
   \bar{c} & = & \sqrt{-\rho (\rho - i)}\text{ is its complex conjugate}\\
   k & = & \frac{2\sqrt{2}\rho^2\sqrt{\rho^2 + 1}}{\pi\sqrt{\rho^2 + \rho\sqrt{\rho^2 + 1}}}
   \end{eqnarray*}

Implementation `scipy.stats.relativistic_bw`
