
.. _continuous-foldnorm:

Folded Normal Distribution
==========================

If :math:`Z` is Normal with mean :math:`L` and :math:`\sigma=S` , then :math:`\left|Z\right|` is a folded normal with shape parameter :math:`c=\left|L\right|/S` , location parameter :math:`0` and scale parameter :math:`S` . This is a special case of the non-central chi distribution with one-
degree of freedom and non-centrality parameter :math:`c^{2}.` Note that :math:`c\geq0` . The standard form of the folded normal is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \sqrt{\frac{2}{\pi}}\cosh\left(cx\right)\exp\left(-\frac{x^{2}+c^{2}}{2}\right)\\ F\left(x;c\right) & = & \Phi\left(x-c\right)-\Phi\left(-x-c\right)=\Phi\left(x-c\right)+\Phi\left(x+c\right)-1\\ G\left(\alpha;c\right) & = & F^{-1}\left(x;c\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\exp\left[\frac{t}{2}\left(t-2c\right)\right]\left(1+e^{2ct}\right)

.. math::
   :nowrap:

    \begin{eqnarray*} k & = & \mathrm{erf}\left(\frac{c}{\sqrt{2}}\right)\\ p & = & \exp\left(-\frac{c^{2}}{2}\right)\\ \mu & = & \sqrt{\frac{2}{\pi}}p+ck\\ \mu_{2} & = & c^{2}+1-\mu^{2}\\ \gamma_{1} & = & \frac{\sqrt{\frac{2}{\pi}}p^{3}\left(4-\frac{\pi}{p^{2}}\left(2c^{2}+1\right)\right)+2ck\left(6p^{2}+3cpk\sqrt{2\pi}+\pi c\left(k^{2}-1\right)\right)}{\pi\mu_{2}^{3/2}}\\ \gamma_{2} & = & \frac{c^{4}+6c^{2}+3+6\left(c^{2}+1\right)\mu^{2}-3\mu^{4}-4p\mu\left(\sqrt{\frac{2}{\pi}}\left(c^{2}+2\right)+\frac{ck}{p}\left(c^{2}+3\right)\right)}{\mu_{2}^{2}}\end{eqnarray*}

Implementation: `scipy.stats.foldnorm`
