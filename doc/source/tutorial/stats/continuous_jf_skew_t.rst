
.. _continuous-jf_skew_t:

Jones and Faddy Skew-T Distribution
===================================

A skew extension of the `t` distribution, defined for :math:`a>0` and
:math:`b>0`.

.. math::
   :nowrap:

    \begin{eqnarray*}
    f(x;a,b) & = & C_{a,b}^{-1} \left(1+\frac{x}{\left(a+b+x^2\right)^{1/2}}\right)^{a+1/2} \left(1-\frac{x}{\left(a+b+x^2\right)^{1/2}}\right)^{b+1/2} \\
    F(x;a,b) & = & I\left(\frac{1+x(a+b+x^2)^{-1/2}}{2};a,b\right) \\
    \mu_{n}^{\prime} & = & \frac{(a+b)^{n/2}}{2^nB(a,b)}\sum_{i=0}^{n}{n \choose i}(-1)^iB\left(a+\frac{n}{2}-i, b-\frac{n}{2}+i\right)
    \end{eqnarray*}

where :math:`C_{a,b}=2^{a+b-1}B(a,b)(a+b)^{1/2}`, :math:`B` is the beta
function `scipy.special.beta` and the formula for the moments
:math:`\mu_{n}^{\prime}` holds provided that :math:`a>n/2` and :math:`b>n/2`.

When :math:`a<b`, the distribution is negatively skewed, and when :math:`a>b`,
the distribution is positively skewed. If :math:`a=b`, then we recover the `t`
distribution with :math:`2a` degrees of freedom.

References
----------

-  M.C. Jones and M.J. Faddy. "A skew extension of the t distribution, with
   applications" *Journal of the Royal Statistical Society*, Series B
   (Statistical Methodology) 65, no. 1 (2003): 159-174.
   :doi:`10.1111/1467-9868.00378`

Implementation: `scipy.stats.jf_skew_t`
