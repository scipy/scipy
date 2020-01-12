
.. _continuous-dgamma:

Double Gamma Distribution
=========================

The double gamma is the signed version of the Gamma distribution. For :math:`\alpha>0:`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha\right) & = & \frac{1}{2\Gamma\left(\alpha\right)}\left|x\right|^{\alpha-1}e^{-\left|x\right|}\\
    F\left(x;\alpha\right) & = & \left\{
      \begin{array}{ccc}
        \frac{1}{2}-\frac{\gamma\left(\alpha,\left|x\right|\right)}{2\Gamma\left(\alpha\right)} &  & x\leq0\\
        \frac{1}{2}+\frac{\gamma\left(\alpha,\left|x\right|\right)}{2\Gamma\left(\alpha\right)} &  & x>0
      \end{array}
    \right.\\
    G\left(q;\alpha\right) & = & \left\{
      \begin{array}{ccc}
        -\gamma^{-1}\left(\alpha,\left|2q-1\right|\Gamma\left(\alpha\right)\right) &  & q\leq\frac{1}{2}\\
        \gamma^{-1}\left(\alpha,\left|2q-1\right|\Gamma\left(\alpha\right)\right) &  & q>\frac{1}{2}
      \end{array}
    \right.\end{eqnarray*}

.. math::

     M\left(t\right)=\frac{1}{2\left(1-t\right)^{a}}+\frac{1}{2\left(1+t\right)^{a}}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu=m_{n} & = & 0\\
    \mu_{2} & = & \alpha\left(\alpha+1\right)\\
    \gamma_{1} & = & 0\\
    \gamma_{2} & = & \frac{\left(\alpha+2\right)\left(\alpha+3\right)}{\alpha\left(\alpha+1\right)}-3\\
    m_{d} & = & \mathrm{NA}\end{eqnarray*}

Implementation: `scipy.stats.dgamma`
