
.. _discrete-planck:

Planck (discrete exponential) Distribution
==========================================

Named Planck because of its relationship to the black-body problem he
solved.

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;\lambda\right) & = & \left(1-e^{-\lambda}\right)e^{-\lambda k}\quad k\lambda\geq0\\ F\left(x;\lambda\right) & = & 1-e^{-\lambda\left(\left\lfloor x\right\rfloor +1\right)}\quad x\lambda\geq0\\ G\left(q;\lambda\right) & = & \left\lceil -\frac{1}{\lambda}\log\left[1-q\right]-1\right\rceil .\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{1}{e^{\lambda}-1}\\ \mu_{2} & = & \frac{e^{-\lambda}}{\left(1-e^{-\lambda}\right)^{2}}\\ \gamma_{1} & = & 2\cosh\left(\frac{\lambda}{2}\right)\\ \gamma_{2} & = & 4+2\cosh\left(\lambda\right)\end{eqnarray*}

.. math::
   :nowrap:

    \[ M\left(t\right)=\frac{1-e^{-\lambda}}{1-e^{t-\lambda}}\]

.. math::
   :nowrap:

    \[ h\left[X\right]=\frac{\lambda e^{-\lambda}}{1-e^{-\lambda}}-\log\left(1-e^{-\lambda}\right)\]

Implementation: `scipy.stats.planck`
