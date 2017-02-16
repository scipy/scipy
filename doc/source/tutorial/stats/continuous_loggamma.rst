
.. _continuous-loggamma:

Log Gamma Distribution
======================

A single shape parameter :math:`c>0` (Defined for all :math:`x` )

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{\exp\left(cx-e^{x}\right)}{\Gamma\left(c\right)}\\ F\left(x;c\right) & = & \frac{\Gamma\left(c,e^{x}\right)}{\Gamma\left(c\right)}\\ G\left(q;c\right) & = & \log\left[\Gamma^{-1}\left[c,q\Gamma\left(c\right)\right]\right]\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\int_{0}^{\infty}\left[\log y\right]^{n}y^{c-1}\exp\left(-y\right)dy.

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu_{1}^{\prime}\\ \mu_{2} & = & \mu_{2}^{\prime}-\mu^{2}\\ \gamma_{1} & = & \frac{\mu_{3}^{\prime}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\mu_{4}^{\prime}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\end{eqnarray*}

Implementation: `scipy.stats.loggamma`
