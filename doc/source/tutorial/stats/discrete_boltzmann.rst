
.. _discrete-boltzmann:

Boltzmann (truncated Planck) Distribution
=========================================

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;N,\lambda\right) & = & \frac{1-e^{-\lambda}}{1-e^{-\lambda N}}\exp\left(-\lambda k\right)\quad k\in\left\{ 0,1,\ldots,N-1\right\} \\ F\left(x;N,\lambda\right) & = & \left\{ \begin{array}{cc} 0 & x<0\\ \frac{1-\exp\left[-\lambda\left(\left\lfloor x\right\rfloor +1\right)\right]}{1-\exp\left(-\lambda N\right)} & 0\leq x\leq N-1\\ 1 & x\geq N-1\end{array}\right.\\ G\left(q,\lambda\right) & = & \left\lceil -\frac{1}{\lambda}\log\left[1-q\left(1-e^{-\lambda N}\right)\right]-1\right\rceil \end{eqnarray*}

Define :math:`z=e^{-\lambda}`

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{z}{1-z}-\frac{Nz^{N}}{1-z^{N}}\\ \mu_{2} & = & \frac{z}{\left(1-z\right)^{2}}-\frac{N^{2}z^{N}}{\left(1-z^{N}\right)^{2}}\\ \gamma_{1} & = & \frac{z\left(1+z\right)\left(\frac{1-z^{N}}{1-z}\right)^{3}-N^{3}z^{N}\left(1+z^{N}\right)}{\left[z\left(\frac{1-z^{N}}{1-z}\right)^{2}-N^{2}z^{N}\right]^{3/2}}\\ \gamma_{2} & = & \frac{z\left(1+4z+z^{2}\right)\left(\frac{1-z^{N}}{1-z}\right)^{4}-N^{4}z^{N}\left(1+4z^{N}+z^{2N}\right)}{\left[z\left(\frac{1-z^{N}}{1-z}\right)^{2}-N^{2}z^{N}\right]^{2}}\end{eqnarray*}

.. math::

    M\left(t\right)=\frac{1-e^{N\left(t-\lambda\right)}}{1-e^{t-\lambda}}\frac{1-e^{-\lambda}}{1-e^{-\lambda N}}

Implementation: `scipy.stats.boltzmann`
