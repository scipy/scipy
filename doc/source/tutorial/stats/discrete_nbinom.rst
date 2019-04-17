
.. _discrete-nbinom:

Negative Binomial Distribution
==============================

The negative binomial random variable with parameters :math:`n` and :math:`p\in\left(0,1\right)` can be defined as the number of *extra* independent trials (beyond :math:`n` ) required to accumulate a total of :math:`n` successes where the probability of a success on each trial is :math:`p.` Equivalently, this random variable is the number of failures
encountered while accumulating :math:`n` successes during independent trials of an experiment that succeeds
with probability :math:`p.` Thus,

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;n,p\right) & = & \left(\begin{array}{c} k+n-1\\ n-1\end{array}\right)p^{n}\left(1-p\right)^{k}\quad k\geq0\\ F\left(x;n,p\right) & = & \sum_{i=0}^{\left\lfloor x\right\rfloor }\left(\begin{array}{c} i+n-1\\ i\end{array}\right)p^{n}\left(1-p\right)^{i}\quad x\geq0\\  & = & I_{p}\left(n,\left\lfloor x\right\rfloor +1\right)\quad x\geq0\\ \mu & = & n\frac{1-p}{p}\\ \mu_{2} & = & n\frac{1-p}{p^{2}}\\ \gamma_{1} & = & \frac{2-p}{\sqrt{n\left(1-p\right)}}\\ \gamma_{2} & = & \frac{p^{2}+6\left(1-p\right)}{n\left(1-p\right)}.\end{eqnarray*}

Recall that :math:`I_{p}\left(a,b\right)` is the incomplete beta integral.

Implementation: `scipy.stats.nbinom`
