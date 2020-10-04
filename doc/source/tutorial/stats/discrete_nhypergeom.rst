
.. _discrete-nhypergeom:

Negative Hypergeometric Distribution
====================================

Consider a box containing :math:`M` balls: :math:`n` red and :math:`M-n` blue. We randomly sample balls from the box, one at a time and *without* replacement, until we have picked :math:`r` blue balls. `nhypergeom` is the distribution of the number of red balls :math:`k` we have picked.

.. math::
   :nowrap:

    \begin{eqnarray*}
    p(k;M,n,r) & = & \frac{\left(\begin{array}{c} k+r-1\\ k\end{array}\right)\left(\begin{array}{c} M-r-k\\ n-k\end{array}\right)}{\left(\begin{array}{c} M\\ n\end{array}\right)}\quad 0 \leq k \leq M-n,\\
    F(x;M,n,r) & = & \sum_{k=0}^{\left\lfloor x\right\rfloor }p\left(k;M,n,r\right),\\
    \mu & = & \frac{rn}{M-n+1},\\
    \mu_{2} & = & \frac{rn(M+1)}{(M-n+1)(M-n+2)}\left(1-\frac{r}{M-n+1}\right)
    \end{eqnarray*}

for :math:`k \in 0, 1, 2, ..., n`, where the binomial coefficients are defined as,

.. math::
   :nowrap:

    \begin{eqnarray*} \binom{n}{k} \equiv \frac{n!}{k! (n - k)!} \end{eqnarray*}

The cumulative distribution, survivor function, hazard function, cumulative hazard function, inverse distribution function, moment generating function, and characteristic function on the support of :math:`k` are mathematically intractable.

Implementation: `scipy.stats.nhypergeom`
