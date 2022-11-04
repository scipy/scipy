
.. _continuous-burr12:

Burr12 Distribution
===================

There are two shape parameters :math:`c,d > 0` and the support is :math:`x \in [0,\infty)`.
The Burr12 distribution is also known as the Singh-Maddala distribution.



.. math::
   :nowrap:

    \begin{eqnarray*}
    f\left(x;c,d\right) & = & {cd} \frac{x^{c-1}} {\left(1+x^{c}\right)^{d+1}} \\
    F\left(x;c,d\right) & = & 1 - \left(1+x^{c}\right)^{-d}\\
    G\left(q;c,d\right) & = & \left((1-q)^{-1/d}-1\right)^{-1/c}\\
    S\left(x;c,d\right) & = & \left(1+x^{c}\right)^{-d}\\
    \mu & = & d \,  B\left(d-\frac{1}{c}, 1+\frac{1}{c}\right)\\
    \mu_{n} & = & d \, B\left(d-\frac{n}{c}, 1+\frac{n}{c}\right)\\
    m_{d} & = & \left(\frac{c-1}{c d + 1}\right)^{1/c} \,\text{if }\quad c>1 \text{, otherwise }\quad 0\\
    m_{n} & = & \left(2^{1/d}-1\right)^{-1/c}
    \end{eqnarray*}

where :math:`B(a, b) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}` is the Beta function.

Implementation: `scipy.stats.burr12`
