
.. _continuous-loglaplace:

Log Double Exponential (Log-Laplace) Distribution
=================================================

Defined over :math:`x>0` with :math:`c>0`

.. math::
   :nowrap:

    \begin{eqnarray*}
        f\left(x;c\right) & = & \left\{
                                    \begin{array}{ccc}
                                        \frac{c}{2}x^{c-1} &  & 0 < x < 1 \\
                                        \frac{c}{2}x^{-c-1} &  & x \geq 1
                                    \end{array}
                                \right. \\
        F\left(x;c\right) & = & \left\{
                                    \begin{array}{ccc}
                                        \frac{1}{2}x^{c} &  & 0 < x < 1 \\
                                        1-\frac{1}{2}x^{-c} &  & x \geq 1
                                    \end{array}
                                \right. \\
        G\left(q;c\right) & = & \left\{
                                    \begin{array}{ccc}
                                        \left(2q\right)^{1/c} &  & 0 \leq q < \frac{1}{2} \\
                                        \left(2-2q\right)^{-1/c} &  & \frac{1}{2} \leq q \leq 1
                                    \end{array}
                                \right.
    \end{eqnarray*}

.. math::

     h\left[X\right]=\log\left(\frac{2e}{c}\right)


Implementation: `scipy.stats.loglaplace`
