
.. _continuous-expon:

Exponential Distribution
========================

This is a special case of the Gamma (and Erlang) distributions with
shape parameter :math:`\left(\alpha=1\right)` and the same location and scale parameters. The standard form is
therefore ( :math:`x\geq0` )

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & e^{-x}\\ F\left(x\right) & = & \Gamma\left(1,x\right)=1-e^{-x}\\ G\left(q\right) & = & -\log\left(1-q\right)\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=n!

.. math::

     M\left(t\right)=\frac{1}{1-t}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & 1\\ \mu_{2} & = & 1\\ \gamma_{1} & = & 2\\ \gamma_{2} & = & 6\\ m_{d} & = & 0\end{eqnarray*}

.. math::

     h\left[X\right]=1.

Implementation: `scipy.stats.expon`
