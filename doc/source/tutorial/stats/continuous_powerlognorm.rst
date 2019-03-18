
.. _continuous-powerlognorm:

Power Log Normal Distribution
=============================

A generalization of the log-normal distribution with shape parameters :math:`\sigma>0`, :math:`c>0` and support :math:`x\geq0`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\sigma,c\right) & = & \frac{c}{x\sigma}\phi\left(\frac{\log x}{\sigma}\right)\left(\Phi\left(-\frac{\log x}{\sigma}\right)\right)^{c-1}\\
    F\left(x;\sigma,c\right) & = & 1-\left(\Phi\left(-\frac{\log x}{\sigma}\right)\right)^{c}\\
    G\left(q;\sigma,c\right) & = & \exp\left(-\sigma\Phi^{-1}\left(\left(1-q\right)^{1/c}\right)\right)\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\int_{0}^{1}\exp\left(-n\sigma\Phi^{-1}\left(y^{1/c}\right)\right)dy

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu_{1}^{\prime}\\
    \mu_{2} & = & \mu_{2}^{\prime}-\mu^{2}\\
    \gamma_{1} & = & \frac{\mu_{3}^{\prime}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\
    \gamma_{2} & = & \frac{\mu_{4}^{\prime}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\end{eqnarray*}

This distribution reduces to the log-normal distribution when :math:`c=1.`

Implementation: `scipy.stats.powerlognorm`
