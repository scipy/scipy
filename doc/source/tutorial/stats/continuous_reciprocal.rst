
.. _continuous-reciprocal:

Reciprocal Distribution
=======================
There are two shape parameters :math:`a,b>0` and the support is :math:`x\in\left[a,b\right]`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,b\right) & = & \frac{1}{x\log\left(b/a\right)}\\ F\left(x;a,b\right) & = & \frac{\log\left(x/a\right)}{\log\left(b/a\right)}\\ G\left(q;a,b\right) & = & a\exp\left(q\log\left(b/a\right)\right)=a\left(\frac{b}{a}\right)^{q}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} d & = & \log\left(a/b\right)\\ \mu & = & \frac{a-b}{d}\\ \mu_{2} & = & \mu\frac{a+b}{2}-\mu^{2}=\frac{\left(a-b\right)\left[a\left(d-2\right)+b\left(d+2\right)\right]}{2d^{2}}\\ \gamma_{1} & = & \frac{\sqrt{2}\left[12d\left(a-b\right)^{2}+d^{2}\left(a^{2}\left(2d-9\right)+2abd+b^{2}\left(2d+9\right)\right)\right]}{3d\sqrt{a-b}\left[a\left(d-2\right)+b\left(d+2\right)\right]^{3/2}}\\ \gamma_{2} & = & \frac{-36\left(a-b\right)^{3}+36d\left(a-b\right)^{2}\left(a+b\right)-16d^{2}\left(a^{3}-b^{3}\right)+3d^{3}\left(a^{2}+b^{2}\right)\left(a+b\right)}{3\left(a-b\right)\left[a\left(d-2\right)+b\left(d+2\right)\right]^{2}}-3\\ m_{d} & = & a\\ m_{n} & = & \sqrt{ab}\end{eqnarray*}

.. math::

     h\left[X\right]=\frac{1}{2}\log\left(ab\right)+\log\left[\log\left(\frac{b}{a}\right)\right].

Implementation: `scipy.stats.reciprocal`
