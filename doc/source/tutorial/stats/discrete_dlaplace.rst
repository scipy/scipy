
.. _discrete-dlaplace:

Discrete Laplacian Distribution
===============================

Defined over all integers for :math:`a>0`

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k\right) & = & \tanh\left(\frac{a}{2}\right)e^{-a\left|k\right|},\\ F\left(x\right) & = & \left\{ \begin{array}{cc} \frac{e^{a\left(\left\lfloor x\right\rfloor +1\right)}}{e^{a}+1} & \left\lfloor x\right\rfloor <0,\\ 1-\frac{e^{-a\left\lfloor x\right\rfloor }}{e^{a}+1} & \left\lfloor x\right\rfloor \geq0.\end{array}\right.\\ G\left(q\right) & = & \left\{ \begin{array}{cc} \left\lceil \frac{1}{a}\log\left[q\left(e^{a}+1\right)\right]-1\right\rceil  & q<\frac{1}{1+e^{-a}},\\ \left\lceil -\frac{1}{a}\log\left[\left(1-q\right)\left(1+e^{a}\right)\right]\right\rceil  & q\geq\frac{1}{1+e^{-a}}.\end{array}\right.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & \tanh\left(\frac{a}{2}\right)\sum_{k=-\infty}^{\infty}e^{tk}e^{-a\left|k\right|}\\  & = & C\left(1+\sum_{k=1}^{\infty}e^{-\left(t+a\right)k}+\sum_{1}^{\infty}e^{\left(t-a\right)k}\right)\\  & = & \tanh\left(\frac{a}{2}\right)\left(1+\frac{e^{-\left(t+a\right)}}{1-e^{-\left(t+a\right)}}+\frac{e^{t-a}}{1-e^{t-a}}\right)\\  & = & \frac{\tanh\left(\frac{a}{2}\right)\sinh a}{\cosh a-\cosh t}.\end{eqnarray*}

Thus,

.. math::
   :nowrap:

    \mu_{n}^{\prime}=M^{\left(n\right)}\left(0\right)=\left[1+\left(-1\right)^{n}\right]\textrm{Li}_{-n}\left(e^{-a}\right)

where :math:`\textrm{Li}_{-n}\left(z\right)` is the polylogarithm function of order :math:`-n` evaluated at :math:`z.`

.. math::
   :nowrap:

    h\left[X\right]=-\log\left(\tanh\left(\frac{a}{2}\right)\right)+\frac{a}{\sinh a}

Implementation: `scipy.stats.dlaplace`
