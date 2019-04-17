
.. _continuous-genlogistic:

Generalized Logistic Distribution
=================================

Has been used in the analysis of extreme values. There is one shape
parameter :math:`c>0.` The support is :math:`x\geq0`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{c\exp\left(-x\right)}{\left[1+\exp\left(-x\right)\right]^{c+1}}\\
    F\left(x;c\right) & = & \frac{1}{\left[1+\exp\left(-x\right)\right]^{c}}\\
    G\left(q;c\right) & = & -\log\left(q^{-1/c}-1\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\frac{c}{1-t}\,_{2}F_{1}\left(1+c,\,1-t\,;\,2-t\,;-1\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \gamma+\psi_{0}\left(c\right)\\
    \mu_{2} & = & \frac{\pi^{2}}{6}+\psi_{1}\left(c\right)\\
    \gamma_{1} & = & \frac{\psi_{2}\left(c\right)+2\zeta\left(3\right)}{\mu_{2}^{3/2}}\\
    \gamma_{2} & = & \frac{\left(\frac{\pi^{4}}{15}+\psi_{3}\left(c\right)\right)}{\mu_{2}^{2}}\\
    m_{d} & = & \log c\\
    m_{n} & = & -\log\left(2^{1/c}-1\right)\end{eqnarray*}

Note that the polygamma function is

.. math::
   :nowrap:

    \begin{eqnarray*} \psi_{n}\left(z\right) & = & \frac{d^{n+1}}{dz^{n+1}}\log\Gamma\left(z\right)\\
    & = & \left(-1\right)^{n+1}n!\sum_{k=0}^{\infty}\frac{1}{\left(z+k\right)^{n+1}}\\
    & = & \left(-1\right)^{n+1}n!\zeta\left(n+1,z\right)\end{eqnarray*}

where :math:`\zeta\left(k,x\right)` is a generalization of the Riemann zeta function called the Hurwitz
zeta function. Note that :math:`\zeta\left(n\right)\equiv\zeta\left(n,1\right)`.

Implementation: `scipy.stats.genlogistic`
