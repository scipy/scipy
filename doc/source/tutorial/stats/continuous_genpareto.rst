
.. _continuous-genpareto:

Generalized Pareto Distribution
===============================

There is one shape parameter :math:`c\neq0`.  The support is :math:`x\geq0` if :math:`c>0`,
and :math:`0\leq x<\frac{1}{\left|c\right|}` if :math:`c` is negative.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \left(1+cx\right)^{-1-\frac{1}{c}}\\
    F\left(x;c\right) & = & 1-\frac{1}{\left(1+cx\right)^{1/c}}\\
    G\left(q;c\right) & = & \frac{1}{c}\left[\left(\frac{1}{1-q}\right)^{c}-1\right]\end{eqnarray*}

.. Comment Is it \gamma\left(-\frac{1}{c},-\frac{t}{c}\right) or \Gamma\left(-\frac{1}{c},-\frac{t}{c}\right) below?


.. math::

    M\left(t\right) = \left\{
      \begin{array}{cc}
        \left(-\frac{t}{c}\right)^{\frac{1}{c}}
        e^{-\frac{t}{c}}
        \left[
        \Gamma\left(1-\frac{1}{c}\right)
        + \left(\gamma\left(-\frac{1}{c},-\frac{t}{c}\right) / \Gamma\left(\frac{1}{-c}\right)\right)
          - \pi\csc\left(\frac{\pi}{c}\right)/\Gamma\left(\frac{1}{c}\right)
          \right] & c>0\\
        \left(
        \frac{\left|c\right|}{t}\right)^{1/\left|c\right|}
        \Gamma\left(\frac{1}{\left|c\right|}, \frac{t}{\left|c\right|}\right)
        \frac{1}{\Gamma\left(\frac{1}{|c|}\right)}
         & c<0
      \end{array}
      \right.

.. math::

     \mu_{n}^{\prime}=\frac{\left(-1\right)^{n}}{c^{n}}\sum_{k=0}^{n}\binom{n}{k}\frac{\left(-1\right)^{k}}{1-ck}\quad \text{ if }cn<1

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{1}^{\prime} & = & \frac{1}{1-c}\quad c<1\\
    \mu_{2}^{\prime} & = & \frac{2}{\left(1-2c\right)\left(1-c\right)}\quad c<\frac{1}{2}\\
    \mu_{3}^{\prime} & = & \frac{6}{\left(1-c\right)\left(1-2c\right)\left(1-3c\right)}\quad c<\frac{1}{3}\\
    \mu_{4}^{\prime} & = & \frac{24}{\left(1-c\right)\left(1-2c\right)\left(1-3c\right)\left(1-4c\right)}\quad c<\frac{1}{4}\end{eqnarray*}

Thus,

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \mu_{1}^{\prime}\\
    \mu_{2} & = & \mu_{2}^{\prime}-\mu^{2}\\
    \gamma_{1} & = & \frac{\mu_{3}^{\prime}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\
    \gamma_{2} & = & \frac{\mu_{4}^{\prime}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\end{eqnarray*}

.. math::

     h\left[X\right]=1+c\quad c>0

Implementation: `scipy.stats.genpareto`
