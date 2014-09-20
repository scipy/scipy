
.. _continuous-gengamma:

Generalized Gamma Distribution
==============================

A general probability form that reduces to many common distributions: :math:`x>0` :math:`a>0` and :math:`c\neq0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,c\right) & = & \frac{\left|c\right|x^{ca-1}}{\Gamma\left(a\right)}\exp\left(-x^{c}\right)\\ F\left(x;a,c\right) & = & \begin{array}{cc} \frac{\Gamma\left(a,x^{c}\right)}{\Gamma\left(a\right)} & c>0\\ 1-\frac{\Gamma\left(a,x^{c}\right)}{\Gamma\left(a\right)} & c<0\end{array}\\ G\left(q;a,c\right) & = & \left\{ \Gamma^{-1}\left[a,\Gamma\left(a\right)q\right]\right\} ^{1/c}\quad c>0\\  &  & \left\{ \Gamma^{-1}\left[a,\Gamma\left(a\right)\left(1-q\right)\right]\right\} ^{1/c}\quad c<0\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\frac{\Gamma\left(a+\frac{n}{c}\right)}{\Gamma\left(a\right)}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\Gamma\left(a+\frac{1}{c}\right)}{\Gamma\left(a\right)}\\ \mu_{2} & = & \frac{\Gamma\left(a+\frac{2}{c}\right)}{\Gamma\left(a\right)}-\mu^{2}\\ \gamma_{1} & = & \frac{\Gamma\left(a+\frac{3}{c}\right)/\Gamma\left(a\right)-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\Gamma\left(a+\frac{4}{c}\right)/\Gamma\left(a\right)-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\\ m_{d} & = & \left(\frac{ac-1}{c}\right)^{1/c}.\end{eqnarray*}

Special cases are Weibull :math:`\left(a=1\right)` , half-normal :math:`\left(a=1/2,c=2\right)` and ordinary gamma distributions :math:`c=1.` If :math:`c=-1` then it is the inverted gamma distribution.

.. math::

     h\left[X\right]=a-a\Psi\left(a\right)+\frac{1}{c}\Psi\left(a\right)+\log\Gamma\left(a\right)-\log\left|c\right|.

Implementation: `scipy.stats.gengamma`
