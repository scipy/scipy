
.. _continuous-f:

Fratio (or F) Distribution
==========================

Defined for :math:`x>0` . The distribution of :math:`\left(X_{1}/X_{2}\right)\left(\nu_{2}/\nu_{1}\right)` if :math:`X_{1}` is chi-squared with :math:`v_{1}` degrees of freedom and :math:`X_{2}` is chi-squared with :math:`v_{2}` degrees of freedom.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\nu_{1},\nu_{2}\right) & = & \frac{\nu_{2}^{\nu_{2}/2}\nu_{1}^{\nu_{1}/2}x^{\nu_{1}/2-1}}{\left(\nu_{2}+\nu_{1}x\right)^{\left(\nu_{1}+\nu_{2}\right)/2}B\left(\frac{\nu_{1}}{2},\frac{\nu_{2}}{2}\right)}\\ F\left(x;v_{1},v_{2}\right) & = & I\left(\frac{\nu_{1}}{2},\frac{\nu_{2}}{2},\frac{\nu_{2}x}{\nu_{2}+\nu_{1}x}\right)\\ G\left(q;\nu_{1},\nu_{2}\right) & = & \left[\frac{\nu_{2}}{I^{-1}\left(\nu_{1}/2,\nu_{2}/2,q\right)}-\frac{\nu_{1}}{\nu_{2}}\right]^{-1}.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\nu_{2}}{\nu_{2}-2}\quad\nu_{2}>2\\ \mu_{2} & = & \frac{2\nu_{2}^{2}\left(\nu_{1}+\nu_{2}-2\right)}{\nu_{1}\left(\nu_{2}-2\right)^{2}\left(\nu_{2}-4\right)}\quad v_{2}>4\\ \gamma_{1} & = & \frac{2\left(2\nu_{1}+\nu_{2}-2\right)}{\nu_{2}-6}\sqrt{\frac{2\left(\nu_{2}-4\right)}{\nu_{1}\left(\nu_{1}+\nu_{2}-2\right)}}\quad\nu_{2}>6\\ \gamma_{2} & = & \frac{3\left[8+\left(\nu_{2}-6\right)\gamma_{1}^{2}\right]}{2\nu-16}\quad\nu_{2}>8\end{eqnarray*}

Implementation: `scipy.stats.f`
