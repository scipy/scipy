
.. _continuous-triang:

Triangular Distribution
=======================

One shape parameter :math:`c\in[0,1]` giving the distance to the peak as a percentage of the total extent of
the non-zero portion. The location parameter is the start of the non-
zero portion, and the scale-parameter is the width of the non-zero
portion. In standard form we have :math:`x\in\left[0,1\right].`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \left\{ \begin{array}{ccc} 2\frac{x}{c} &  & x<c\\ 2\frac{1-x}{1-c} &  & x\geq c\end{array}\right.\\ F\left(x;c\right) & = & \left\{ \begin{array}{ccc} \frac{x^{2}}{c} &  & x<c\\ \frac{x^{2}-2x+c}{c-1} &  & x\geq c\end{array}\right.\\ G\left(q;c\right) & = & \left\{ \begin{array}{ccc} \sqrt{cq} &  & q<c\\ 1-\sqrt{\left(1-c\right)\left(1-q\right)} &  & q\geq c\end{array}\right.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{c}{3}+\frac{1}{3}\\ \mu_{2} & = & \frac{1-c+c^{2}}{18}\\ \gamma_{1} & = & \frac{\sqrt{2}\left(2c-1\right)\left(c+1\right)\left(c-2\right)}{5\left(1-c+c^{2}\right)^{3/2}}\\ \gamma_{2} & = & -\frac{3}{5}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left(X\right) & = & \log\left(\frac{1}{2}\sqrt{e}\right)\\  & \approx & -0.19314718055994530942.\end{eqnarray*}

Implementation: `scipy.stats.triang`
