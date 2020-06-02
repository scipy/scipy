
.. _discrete-nhypergeom:

Negative Hypergeometric Distribution
====================================

We choose one element from the box without replacement until we end up with r elements which we consider failures (like getting r blue balls while targeting for red balls) and then we count how many successes k we encountered (like getting k red balls (successes) in a sample with r blue balls (failures)). In the end, there are k+r samples collected without replacement and nhypergeom is the distribution over this k.

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;N,K,r\right) & = & \frac{\left(\begin{array}{c} k+r-1\\ k\end{array}\right)\left(\begin{array}{c} N-r-k\\ K-k\end{array}\right)}{\left(\begin{array}{c} N\\ K\end{array}\right)}\quad 0 \leq k \leq N-K\\ F\left(x;N,K,r\right) & = & \sum_{k=0}^{\left\lfloor x\right\rfloor }\frac{\left(\begin{array}{c} k+r-1\\ k\end{array}\right)\left(\begin{array}{c} N-r-k\\ K-k\end{array}\right)}{\left(\begin{array}{c} N\\ K\end{array}\right)},\\ \mu & = & \frac{rK}{N-K+1}\\ \mu_{2} & = & \frac{rK\left(N+1\right)}{\left(N-K+1\right)\left(N-K+2\right)}\left(1-\frac{r}{N-K+1}\right)\end{eqnarray*}

for :math:`k \in 0, 1, 2, ..., K`, where the binomial coefficients are defined as,

.. math::
   :nowrap:

    \begin{eqnarray*} \binom{n}{k} \equiv \frac{n!}{k! (n - k)!} \end{eqnarray*}

The cumulative distribution, survivor function, hazard function, cumulative hazard function, and inverse distribution function, moment generating function, and characteristic function on the support of :math:`k` are mathematically intractable.

Implementation: `scipy.stats.nhypergeom`
