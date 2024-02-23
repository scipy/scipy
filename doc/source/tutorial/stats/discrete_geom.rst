
.. _discrete-geom:

Geometric Distribution
======================

The geometric random variable with parameter :math:`p\in\left(0,1\right)` can be defined as the number of trials required to obtain a success
where the probability of success on each trial is :math:`p` . Thus,

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;p\right) & = & \left(1-p\right)^{k-1}p\quad k\geq1\\ F\left(x;p\right) & = & 1-\left(1-p\right)^{\left\lfloor x\right\rfloor }\quad x\geq1\\ G\left(q;p\right) & = & \left\lceil \frac{\log\left(1-q\right)}{\log\left(1-p\right)}\right\rceil \\ \mu & = & \frac{1}{p}\\ \mu_{2} & = & \frac{1-p}{p^{2}}\\ \gamma_{1} & = & \frac{2-p}{\sqrt{1-p}}\\ \gamma_{2} & = & \frac{p^{2}-6p+6}{1-p}.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & \frac{p}{e^{-t}-\left(1-p\right)}\end{eqnarray*}

Implementation: `scipy.stats.geom`
