
.. _discrete-randint:

Discrete Uniform (randint) Distribution
=======================================

The discrete uniform distribution with parameters :math:`\left(a,b\right)` constructs a random variable that has an equal probability of being
any one of the integers in the half-open range :math:`[a,b)`. If :math:`a` is not given it is assumed to be zero and the only parameter is :math:`b`. Therefore,

.. math::
   :nowrap:

    \begin{eqnarray*}
        p\left(k,a,b\right) & = & \frac{1}{b-a} \quad a \leq k < b \\
        F\left(x;a,b\right) & = & \frac{\left\lfloor x\right\rfloor -a}{b-a} \quad a \leq x \leq b \\
        G\left(q;a,b\right) & = & \left\lceil q\left(b-a\right)+a\right\rceil \\
        \mu & = & \frac{b+a-1}{2}\\
        \mu_{2} & = & \frac{\left(b-a-1\right)\left(b-a+1\right)}{12}\\
        \gamma_{1} & = & 0 \\
        \gamma_{2} & = & -\frac{6}{5}\frac{\left(b-a\right)^{2}+1}{\left(b-a-1\right)\left(b-a+1\right)}.
    \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*}
        M\left(t\right) & = & \frac{1}{b-a}\sum_{k=a}^{b-1}e^{tk}\\
                        & = & \frac{e^{bt}-e^{at}}{\left(b-a\right)\left(e^{t}-1\right)}
    \end{eqnarray*}

Implementation: `scipy.stats.randint`
