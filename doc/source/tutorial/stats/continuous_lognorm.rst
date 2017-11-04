
.. _continuous-lognorm:

Log Normal (Cobb-Douglass) Distribution
=======================================

Has one shape parameter :math:`\sigma` >0. (Notice that the "Regress ":math:`A=\log S` where :math:`S` is the scale parameter and :math:`A` is the mean of the underlying normal distribution). The standard form
is :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\sigma\right) & = & \frac{1}{\sigma x\sqrt{2\pi}}\exp\left[-\frac{1}{2}\left(\frac{\log x}{\sigma}\right)^{2}\right]\\ F\left(x;\sigma\right) & = & \Phi\left(\frac{\log x}{\sigma}\right)\\ G\left(q;\sigma\right) & = & \exp\left\{ \sigma\Phi^{-1}\left(q\right)\right\} \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \exp\left(\sigma^{2}/2\right)\\ \mu_{2} & = & \exp\left(\sigma^{2}\right)\left[\exp\left(\sigma^{2}\right)-1\right]\\ \gamma_{1} & = & \sqrt{p-1}\left(2+p\right)\\ \gamma_{2} & = & p^{4}+2p^{3}+3p^{2}-6\quad\quad p=e^{\sigma^{2}}\end{eqnarray*}

Notice that using JKB notation we have :math:`\theta=L,` :math:`\zeta=\log S` and we have given the so-called antilognormal form of the
distribution. This is more consistent with the location, scale
parameter description of general probability distributions.

.. math::

     h\left[X\right]=\frac{1}{2}\left[1+\log\left(2\pi\right)+2\log\left(\sigma\right)\right].

Also, note that if :math:`X` is a log-normally distributed random-variable with :math:`L=0` and :math:`S` and shape parameter :math:`\sigma.` Then, :math:`\log X` is normally distributed with variance :math:`\sigma^{2}` and mean :math:`\log S.`

Implementation: `scipy.stats.lognorm`
