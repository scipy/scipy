
.. _continuous-irwinhall:

Irwin-Hall ( sum of `n` independent uniforms on (0,1) ) Distribution
====================================================================

.. math::
   :nowrap:

    \begin{eqnarray*}
    F(x) &=& \frac{1}{n!}\sum_{k=0}^{\lfloor x\rfloor}(-1)^k \binom{n}{k}(x-k)^n
    \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=m_{n}=\mu & = & \frac{n}{2}\\ \mu_{2} & = & \frac{n}{12}\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & \frac{-6}{5n}\end{eqnarray*}

The location parameter of `irwinhall` is

.. math::

     \hat{L}=0.0

where :math:`X_{i}` is a sequence of :math:`N` mutually independent Irwin-Hall  R.V.'s.
Also, the ML estimator of the scale parameter is

.. math::

     \hat{S}=\mathrm{max}_i\left(X_{i}\right).

Implementation: `scipy.stats.irwinhall`
