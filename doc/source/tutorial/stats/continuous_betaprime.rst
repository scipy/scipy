
.. _continuous-betaprime:

Beta Prime Distribution
=======================

Defined over :math:`0<x<\infty.` :math:`\alpha,\beta>0.` (Note the CDF evaluation uses Eq. 3.194.1 on pg. 313 of Gradshteyn &
Ryzhik (sixth edition).

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha,\beta\right) & = & \frac{\Gamma\left(\alpha+\beta\right)}{\Gamma\left(\alpha\right)\Gamma\left(\beta\right)}x^{\alpha-1}\left(1+x\right)^{-\alpha-\beta}\\ F\left(x;\alpha,\beta\right) & = & \frac{\Gamma\left(\alpha+\beta\right)}{\alpha\Gamma\left(\alpha\right)\Gamma\left(\beta\right)}x^{\alpha}\,_{2}F_{1}\left(\alpha+\beta,\alpha;1+\alpha;-x\right)\\ G\left(q;\alpha,\beta\right) & = & F^{-1}\left(x;\alpha,\beta\right)\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\left\{ \begin{array}{ccc} \frac{\Gamma\left(n+\alpha\right)\Gamma\left(\beta-n\right)}{\Gamma\left(\alpha\right)\Gamma\left(\beta\right)}=\frac{\left(\alpha\right)_{n}}{\left(\beta-n\right)_{n}} &  & \beta>n\\ \infty &  & \mathrm{otherwise}\end{array}\right.

Therefore,

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\alpha}{\beta-1}\quad\beta>1\\ \mu_{2} & = & \frac{\alpha\left(\alpha+1\right)}{\left(\beta-2\right)\left(\beta-1\right)}-\frac{\alpha^{2}}{\left(\beta-1\right)^{2}}\quad\beta>2\\ \gamma_{1} & = & \frac{\frac{\alpha\left(\alpha+1\right)\left(\alpha+2\right)}{\left(\beta-3\right)\left(\beta-2\right)\left(\beta-1\right)}-3\mu\mu_{2}-\mu^{3}}{\mu_{2}^{3/2}}\quad\beta>3\\ \gamma_{2} & = & \frac{\mu_{4}}{\mu_{2}^{2}}-3\\ \mu_{4} & = & \frac{\alpha\left(\alpha+1\right)\left(\alpha+2\right)\left(\alpha+3\right)}{\left(\beta-4\right)\left(\beta-3\right)\left(\beta-2\right)\left(\beta-1\right)}-4\mu\mu_{3}-6\mu^{2}\mu_{2}-\mu^{4}\quad\beta>4\end{eqnarray*}

Implementation: `scipy.stats.betaprime`
