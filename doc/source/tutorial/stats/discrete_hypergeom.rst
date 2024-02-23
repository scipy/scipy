
.. _discrete-hypergeom:

Hypergeometric Distribution
===========================

The hypergeometric random variable with parameters :math:`\left(M,n,N\right)` counts the number of "good "objects in a sample of size :math:`N` chosen without replacement from a population of :math:`M` objects where :math:`n` is the number of "good "objects in the total population.

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;N,n,M\right) & = & \frac{\left(\begin{array}{c} n\\ k\end{array}\right)\left(\begin{array}{c} M-n\\ N-k\end{array}\right)}{\left(\begin{array}{c} M\\ N\end{array}\right)}\quad N-\left(M-n\right)\leq k\leq\min\left(n,N\right)\\ F\left(x;N,n,M\right) & = & \sum_{k=0}^{\left\lfloor x\right\rfloor }\frac{\left(\begin{array}{c} m\\ k\end{array}\right)\left(\begin{array}{c} N-m\\ n-k\end{array}\right)}{\left(\begin{array}{c} N\\ n\end{array}\right)},\\ \mu & = & \frac{nN}{M}\\ \mu_{2} & = & \frac{nN\left(M-n\right)\left(M-N\right)}{M^{2}\left(M-1\right)}\\ \gamma_{1} & = & \frac{\left(M-2n\right)\left(M-2N\right)}{M-2}\sqrt{\frac{M-1}{nN\left(M-m\right)\left(M-n\right)}}\\ \gamma_{2} & = & \frac{g\left(N,n,M\right)}{nN\left(M-n\right)\left(M-3\right)\left(M-2\right)\left(N-M\right)}\end{eqnarray*}

where (defining :math:`m=M-n` )

.. math::
   :nowrap:

    \begin{eqnarray*} g\left(N,n,M\right) & = & m^{3}-m^{5}+3m^{2}n-6m^{3}n+m^{4}n+3mn^{2}\\  &  & -12m^{2}n^{2}+8m^{3}n^{2}+n^{3}-6mn^{3}+8m^{2}n^{3}\\  &  & +mn^{4}-n^{5}-6m^{3}N+6m^{4}N+18m^{2}nN\\  &  & -6m^{3}nN+18mn^{2}N-24m^{2}n^{2}N-6n^{3}N\\  &  & -6mn^{3}N+6n^{4}N+6m^{2}N^{2}-6m^{3}N^{2}-24mnN^{2}\\  &  & +12m^{2}nN^{2}+6n^{2}N^{2}+12mn^{2}N^{2}-6n^{3}N^{2}.\end{eqnarray*}

Implementation: `scipy.stats.hypergeom`
