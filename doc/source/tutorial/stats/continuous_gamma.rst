
.. _continuous-gamma:

Gamma Distribution
==================

The standard form for the gamma distribution is :math:`\left(\alpha>0\right)` valid for :math:`x\geq0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha\right) & = & \frac{1}{\Gamma\left(\alpha\right)}x^{\alpha-1}e^{-x}\\
    F\left(x;\alpha\right) & = & \frac{\gamma\left(\alpha,x\right)}{\Gamma(\alpha)}\\
    G\left(q;\alpha\right) & = & \gamma^{-1}\left(\alpha,q\Gamma(\alpha)\right)\end{eqnarray*}

where :math:`\gamma` is the lower incomplete gamma function, :math:`\gamma\left(s, x\right) = \int_0^x t^{s-1} e^{-t} dt`.

.. math::

     M\left(t\right)=\frac{1}{\left(1-t\right)^{\alpha}}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \alpha\\
    \mu_{2} & = & \alpha\\
    \gamma_{1} & = & \frac{2}{\sqrt{\alpha}}\\
    \gamma_{2} & = & \frac{6}{\alpha}\\
    m_{d} & = & \alpha-1\end{eqnarray*}

.. math::

     h\left[X\right]=\Psi\left(a\right)\left[1-a\right]+a+\log\Gamma\left(a\right)

where

.. math::

     \Psi\left(a\right)=\frac{\Gamma^{\prime}\left(a\right)}{\Gamma\left(a\right)}.

Implementation: `scipy.stats.gamma`
