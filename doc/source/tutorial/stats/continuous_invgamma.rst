
.. _continuous-invgamma:

Inverted Gamma Distribution
===========================

Special case of the generalized Gamma distribution with :math:`c=-1` and :math:`a>0` , :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a\right) & = & \frac{x^{-a-1}}{\Gamma\left(a\right)}\exp\left(-\frac{1}{x}\right)\\ F\left(x;a\right) & = & \frac{\Gamma\left(a,\frac{1}{x}\right)}{\Gamma\left(a\right)}\\ G\left(q;a\right) & = & \left\{ \Gamma^{-1}\left[a,\Gamma\left(a\right)q\right]\right\} ^{-1}\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\frac{\Gamma\left(a-n\right)}{\Gamma\left(a\right)}\quad a>n

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{1}{a-1}\quad a>1\\ \mu_{2} & = & \frac{1}{\left(a-2\right)\left(a-1\right)}-\mu^{2}\quad a>2\\ \gamma_{1} & = & \frac{\frac{1}{\left(a-3\right)\left(a-2\right)\left(a-1\right)}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{\frac{1}{\left(a-4\right)\left(a-3\right)\left(a-2\right)\left(a-1\right)}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}}{\mu_{2}^{2}}-3\end{eqnarray*}

.. math::

     m_{d}=\frac{1}{a+1}

.. math::

     h\left[X\right]=a-\left(a+1\right)\Psi\left(a\right)+\log\Gamma\left(a\right).

Implementation: `scipy.stats.invgamma`
