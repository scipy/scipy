
.. _discrete-poisson:

Poisson Distribution
====================

The Poisson random variable counts the number of successes in :math:`n` independent Bernoulli trials in the limit as :math:`n\rightarrow\infty` and :math:`p\rightarrow0` where the probability of success in each trial is :math:`p` and :math:`np=\lambda\geq0` is a constant. It can be used to approximate the Binomial random
variable or in it's own right to count the number of events that occur
in the interval :math:`\left[0,t\right]` for a process satisfying certain "sparsity "constraints. The functions are

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;\lambda\right) & = & e^{-\lambda}\frac{\lambda^{k}}{k!}\quad k\geq0,\\ F\left(x;\lambda\right) & = & \sum_{n=0}^{\left\lfloor x\right\rfloor }e^{-\lambda}\frac{\lambda^{n}}{n!}=\frac{1}{\Gamma\left(\left\lfloor x\right\rfloor +1\right)}\int_{\lambda}^{\infty}t^{\left\lfloor x\right\rfloor }e^{-t}dt,\\ \mu & = & \lambda\\ \mu_{2} & = & \lambda\\ \gamma_{1} & = & \frac{1}{\sqrt{\lambda}}\\ \gamma_{2} & = & \frac{1}{\lambda}.\end{eqnarray*}

.. math::
   :nowrap:

    \[ M\left(t\right)=\exp\left[\lambda\left(e^{t}-1\right)\right].\]

Implementation: `scipy.stats.poisson`
