
.. _discrete-bernoulli:

Bernoulli Distribution
======================

A Bernoulli random variable of parameter :math:`p` takes one of only two values :math:`X=0` or :math:`X=1` . The probability of success ( :math:`X=1` ) is :math:`p` , and the probability of failure ( :math:`X=0` ) is :math:`1-p.` It can be thought of as a binomial random variable with :math:`n=1` . The PMF is :math:`p\left(k\right)=0` for :math:`k\neq0,1` and

.. math::
   :nowrap:

    \begin{eqnarray*}
        p\left(k;p\right) & = & \begin{cases} 1-p & k=0\\ p & k=1\end{cases}\\
        F\left(x;p\right) & = & \begin{cases} 0 & x<0\\ 1-p & 0\le x<1\\ 1 & 1\leq x\end{cases}\\
        G\left(q;p\right) & = & \begin{cases} 0 & 0\leq q<1-p\\ 1 & 1-p\leq q\leq1\end{cases}\\
        \mu & = & p\\ \mu_{2} & = & p\left(1-p\right)\\
        \gamma_{3} & = & \frac{1-2p}{\sqrt{p\left(1-p\right)}}\\
        \gamma_{4} & = & \frac{1-6p\left(1-p\right)}{p\left(1-p\right)}
    \end{eqnarray*}

.. math::

    M\left(t\right) = 1-p\left(1-e^{t}\right)

.. math::

    \mu_{m}^{\prime}=p

.. math::

    h\left[X\right]=p\log p+\left(1-p\right)\log\left(1-p\right)

Implementation: `scipy.stats.bernoulli`
