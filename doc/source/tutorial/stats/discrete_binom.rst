
.. _discrete-binom:

Binomial Distribution
=====================

A binomial random variable with parameters :math:`\left(n,p\right)` can be described as the sum of :math:`n` independent Bernoulli random variables of parameter :math:`p;`

.. math::
   :nowrap:

    Y=\sum_{i=1}^{n}X_{i}.

Therefore, this random variable counts the number of successes in :math:`n` independent trials of a random experiment where the probability of
success is :math:`p.`

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;n,p\right) & = & \left(\begin{array}{c} n\\ k\end{array}\right)p^{k}\left(1-p\right)^{n-k}\,\, k\in\left\{ 0,1,\ldots n\right\} ,\\ F\left(x;n,p\right) & = & \sum_{k\leq x}\left(\begin{array}{c} n\\ k\end{array}\right)p^{k}\left(1-p\right)^{n-k}=I_{1-p}\left(n-\left\lfloor x\right\rfloor ,\left\lfloor x\right\rfloor +1\right)\quad x\geq0\end{eqnarray*}

where the incomplete beta integral is

.. math::
   :nowrap:

    I_{x}\left(a,b\right)=\frac{\Gamma\left(a+b\right)}{\Gamma\left(a\right)\Gamma\left(b\right)}\int_{0}^{x}t^{a-1}\left(1-t\right)^{b-1}dt.

Now

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & np\\ \mu_{2} & = & np\left(1-p\right)\\ \gamma_{1} & = & \frac{1-2p}{\sqrt{np\left(1-p\right)}}\\ \gamma_{2} & = & \frac{1-6p\left(1-p\right)}{np\left(1-p\right)}.\end{eqnarray*}

.. math::
   :nowrap:

    M\left(t\right)=\left[1-p\left(1-e^{t}\right)\right]^{n}

Implementation: `scipy.stats.binom`
