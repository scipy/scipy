
.. _continuous-dweibull:

Double Weibull Distribution
===========================

This is a signed form of the Weibull distribution.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{c}{2}\left|x\right|^{c-1}\exp\left(-\left|x\right|^{c}\right)\\ F\left(x;c\right) & = & \left\{ \begin{array}{ccc} \frac{1}{2}\exp\left(-\left|x\right|^{c}\right) &  & x\leq0\\ 1-\frac{1}{2}\exp\left(-\left|x\right|^{c}\right) &  & x>0\end{array}\right.\\ G\left(q;c\right) & = & \left\{ \begin{array}{ccc} -\log^{1/c}\left(\frac{1}{2q}\right) &  & q\leq\frac{1}{2}\\ \log^{1/c}\left(\frac{1}{2q-1}\right) &  & q>\frac{1}{2}\end{array}\right.\end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\mu_{n}=\begin{cases} \Gamma\left(1+\frac{n}{c}\right) & n\mathrm{ even}\\ 0 & n\mathrm{ odd}\end{cases}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=\mu & = & 0\\ \mu_{2} & = & \Gamma\left(\frac{c+2}{c}\right)\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & \frac{\Gamma\left(1+\frac{4}{c}\right)}{\Gamma^{2}\left(1+\frac{2}{c}\right)}\\ m_{d} & = & \mathrm{NA bimodal}\end{eqnarray*}

Implementation: `scipy.stats.dweibull`
