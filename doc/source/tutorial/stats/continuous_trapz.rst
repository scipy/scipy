
.. _continuous-trapz:

Trapezoidal Distribution
========================

Two shape parameters :math:`c\in[0,1], d\in[0, 1]` giving the distances to the
first and second modes as a percentage of the total extent of
the non-zero portion. The location parameter is the start of the non-
zero portion, and the scale-parameter is the width of the non-zero
portion. In standard form we have :math:`x\in\left[0,1\right].`

.. math::
   :nowrap:

    \begin{eqnarray*}
        u(c, d) & = & \frac{2}{d - c + 1} \\
        f\left(x;c, d\right) & = & \left\{
                                    \begin{array}{ccc}
                                        \frac{ux}{c} &  & x < c \\
                                        u & & c\leq x \leq d \\
                                        u\frac{1-x}{1-d} &  & x > d \\
                                    \end{array}
                                \right.\\
        F\left(x;c, d\right) & = & \left\{
                                    \begin{array}{ccc}
                                        \frac{ux^{2}}{2c} &  & x < c \\
                                        \frac{uc}{2} + u(x-c) &  & c\leq x \leq d \\
                                        1 - \frac{u(1 - x)^2}{2(1 - d)} &  & x > d \\
                                    \end{array}
                                \right.\\
        G\left(q;c, d\right) & = & \left\{
                                    \begin{array}{ccc}
                                        \sqrt{qc(d-c+1)} &  & q < c \\
                                        \frac{q}{u}+ \frac{c}{2} &  & q \leq d \\
                                        1 - \sqrt{\frac{2(1 - q) (1 - d)}{u}} &  & q > d \\
                                    \end{array}
                                \right.
    \end{eqnarray*}


Implementation: `scipy.stats.trapz`
