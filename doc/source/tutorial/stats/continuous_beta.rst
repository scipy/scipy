
.. _continuous-beta:

Beta Distribution
=================

There are two shape parameters :math:`a,b > 0` and the support is :math:`x\in[0,1]`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,b\right) & = & \frac{\Gamma\left(a+b\right)}{\Gamma\left(a\right)\Gamma\left(b\right)}x^{a-1}\left(1-x\right)^{b-1} \\
    F\left(x;a,b\right) & = & \int_{0}^{x}f\left(y;a,b\right)dy=I\left(x;a,b\right)\\
    G\left(q;a,b\right) & = & I^{-1}\left(q;a,b\right)\\
    M\left(t\right) & = & \frac{\Gamma\left(a\right)\Gamma\left(b\right)}{\Gamma\left(a+b\right)}\,_{1}F_{1}\left(a;a+b;t\right)\\
    \mu & = & \frac{a}{a+b}\\
    \mu_{2} & = & \frac{ab\left(a+b+1\right)}{\left(a+b\right)^{2}}\\
    \gamma_{1} & = & 2\frac{b-a}{a+b+2}\sqrt{\frac{a+b+1}{ab}}\\
    \gamma_{2} & = & \frac{6\left(a^{3}+a^{2}\left(1-2b\right)+b^{2}\left(b+1\right)-2ab\left(b+2\right)\right)}{ab\left(a+b+2\right)\left(a+b+3\right)}\\
    m_{d} & = & \frac{\left(a-1\right)}{\left(a+b-2\right)}\, a+b\neq2\end{eqnarray*}

where :math:`I\left(x;a,b\right)` is the regularized incomplete Beta function.  :math:`f\left(x;a,1\right)` is also called the Power-function distribution.

.. math::

     l_{\mathbf{x}}\left(a,b\right)=-N\log\Gamma\left(a+b\right)+N\log\Gamma\left(a\right)+N\log\Gamma\left(b\right)-N\left(a-1\right)\overline{\log\mathbf{x}}-N\left(b-1\right)\overline{\log\left(1-\mathbf{x}\right)}

Implementation: `scipy.stats.beta`
