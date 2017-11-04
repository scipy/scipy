
.. _continuous-pareto:

Pareto Distribution
===================

For :math:`x\geq1` and :math:`b>0` . Standard form is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;b\right) & = & \frac{b}{x^{b+1}}\\ F\left(x;b\right) & = & 1-\frac{1}{x^{b}}\\ G\left(q;b\right) & = & \left(1-q\right)^{-1/b}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{b}{b-1}\quad b>1\\ \mu_{2} & = & \frac{b}{\left(b-2\right)\left(b-1\right)^{2}}\quad b>2\\ \gamma_{1} & = & \frac{2\left(b+1\right)\sqrt{b-2}}{\left(b-3\right)\sqrt{b}}\quad b>3\\ \gamma_{2} & = & \frac{6\left(b^{3}+b^{2}-6b-2\right)}{b\left(b^{2}-7b+12\right)}\quad b>4\end{eqnarray*}

.. math::

     h\left(X\right)=\frac{1}{c}+1-\log\left(c\right)

Implementation: `scipy.stats.pareto`
