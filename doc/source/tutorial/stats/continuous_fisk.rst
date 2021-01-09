
.. _continuous-fisk:

Fisk (Log Logistic) Distribution
================================

Special case of the Burr distribution with :math:`d=1`.
There is are one shape parameter :math:`c > 0` and the support is :math:`x \in [0,\infty)`.

.. math::
   :nowrap:

    \begin{eqnarray*}\textrm{Let }k & = & \Gamma\left(1-\frac{2}{c}\right)\Gamma\left(\frac{2}{c}+1\right)-\Gamma^{2}\left(1-\frac{1}{c}\right)\Gamma^{2}\left(\frac{1}{c}+1\right)\\
    f\left(x;c,d\right) & = & \frac{cx^{c-1}}{\left(1+x^{c}\right)^{2}} \\
    F\left(x;c,d\right) & = & \left(1+x^{-c}\right)^{-1}\\
    G\left(q;c,d\right) & = & \left(q^{-1}-1\right)^{-1/c}\\
    \mu & = & \Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+1\right)\\
    \mu_{2} & = & k\\
    \gamma_{1} & = & \frac{1}{\sqrt{k^{3}}}\left[2\Gamma^{3}\left(1-\frac{1}{c}\right)\Gamma^{3}\left(\frac{1}{c}+1\right)+\Gamma\left(1-\frac{3}{c}\right)\Gamma\left(\frac{3}{c}+1\right)\right.\\  &  & \left.-3\Gamma\left(1-\frac{2}{c}\right)\Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+1\right)\Gamma\left(\frac{2}{c}+1\right)\right]\\
    \gamma_{2} & = & -3+\frac{1}{k^{2}}\left[6\Gamma\left(1-\frac{2}{c}\right)\Gamma^{2}\left(1-\frac{1}{c}\right)\Gamma^{2}\left(\frac{1}{c}+1\right)\Gamma\left(\frac{2}{c}+1\right)\right.\\  &  & -3\Gamma^{4}\left(1-\frac{1}{c}\right)\Gamma^{4}\left(\frac{1}{c}+1\right)+\Gamma\left(1-\frac{4}{c}\right)\Gamma\left(\frac{4}{c}+1\right)\\  &  & \left.-4\Gamma\left(1-\frac{3}{c}\right)\Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+1\right)\Gamma\left(\frac{3}{c}+1\right)\right]\\
    m_{d} & = & \left(\frac{c-1}{c+1}\right)^{1/c}\, \text{if }c>1, \text{otherwise } 0\\
    m_{n} & = & 1\\
    h\left[X\right] & = & 2-\log c\end{eqnarray*}


Implementation: `scipy.stats.fisk`
