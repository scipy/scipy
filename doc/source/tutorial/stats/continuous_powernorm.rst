
.. _continuous-powernorm:

Power Normal Distribution
=========================

A generalization of the normal distribution, :math:`c>0` for

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & c\phi\left(x\right)\left(\Phi\left(-x\right)\right)^{c-1}\\ F\left(x;c\right) & = & 1-\left(\Phi\left(-x\right)\right)^{c}\\ G\left(q;c\right) & = & -\Phi^{-1}\left[\left(1-q\right)^{1/c}\right]\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\left(-1\right)^{n}\int_{0}^{1}\left[\Phi^{-1}\left(y^{1/c}\right)\right]^{n}dy

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu_{1}^{\prime}\\ \mu_{2} & = & \mu_{2}^{\prime}-\mu^{2}\\ \gamma_{1} & = & \frac{\mu_{3}^{\prime}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\mu_{4}^{\prime}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\end{eqnarray*}

For :math:`c=1` this reduces to the normal distribution.

Implementation: `scipy.stats.powernorm`
