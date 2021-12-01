
.. _continuous-gausshyper:

Gauss Hypergeometric Distribution
=================================

The four shape parameters are :math:`\alpha>0`, :math:`\beta>0`,
:math:`-\infty < \gamma < \infty`, and :math:`z > -1`.
The support is :math:`x\in\left[0,1\right]`.

.. math::

     \text{Let }C=\frac{1}{B\left(\alpha,\beta\right)\,_{2}F_{1}\left(\gamma,\alpha;\alpha+\beta;-z\right)}

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha,\beta,\gamma,z\right) & = & Cx^{\alpha-1}\frac{\left(1-x\right)^{\beta-1}}{\left(1+zx\right)^{\gamma}}\\ \mu_{n}^{\prime} & = & \frac{B\left(n+\alpha,\beta\right)}{B\left(\alpha,\beta\right)}\frac{\,_{2}F_{1}\left(\gamma,\alpha+n;\alpha+\beta+n;-z\right)}{\,_{2}F_{1}\left(\gamma,\alpha;\alpha+\beta;-z\right)}\end{eqnarray*}

Implementation: `scipy.stats.gausshyper`
