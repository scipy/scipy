
.. _discrete-zipf:

Zipf (Zeta) Distribution
========================

A random variable has the zeta distribution (also called the zipf
distribution) with parameter :math:`\alpha>1` if it's probability mass function is given by

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;\alpha\right) & = & \frac{1}{\zeta\left(\alpha\right)k^{\alpha}}\quad k\geq1\end{eqnarray*}

where

.. math::

    \zeta\left(\alpha\right)=\sum_{n=1}^{\infty}\frac{1}{n^{\alpha}}

is the Riemann zeta function. Other functions of this distribution are

.. math::
   :nowrap:

    \begin{eqnarray*} F\left(x;\alpha\right) & = & \frac{1}{\zeta\left(\alpha\right)}\sum_{k=1}^{\left\lfloor x\right\rfloor }\frac{1}{k^{\alpha}}\\ \mu & = & \frac{\zeta_{1}}{\zeta_{0}}\quad\alpha>2\\ \mu_{2} & = & \frac{\zeta_{2}\zeta_{0}-\zeta_{1}^{2}}{\zeta_{0}^{2}}\quad\alpha>3\\ \gamma_{1} & = & \frac{\zeta_{3}\zeta_{0}^{2}-3\zeta_{0}\zeta_{1}\zeta_{2}+2\zeta_{1}^{3}}{\left[\zeta_{2}\zeta_{0}-\zeta_{1}^{2}\right]^{3/2}}\quad\alpha>4\\ \gamma_{2} & = & \frac{\zeta_{4}\zeta_{0}^{3}-4\zeta_{3}\zeta_{1}\zeta_{0}^{2}+12\zeta_{2}\zeta_{1}^{2}\zeta_{0}-6\zeta_{1}^{4}-3\zeta_{2}^{2}\zeta_{0}^{2}}{\left(\zeta_{2}\zeta_{0}-\zeta_{1}^{2}\right)^{2}}.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & \frac{\textrm{Li}_{\alpha}\left(e^{t}\right)}{\zeta\left(\alpha\right)}\end{eqnarray*}

where :math:`\zeta_{i}=\zeta\left(\alpha-i\right)` and :math:`\textrm{Li}_{n}\left(z\right)` is the :math:`n^{\textrm{th}}` polylogarithm function of :math:`z` defined as

.. math::

    \textrm{Li}_{n}\left(z\right)\equiv\sum_{k=1}^{\infty}\frac{z^{k}}{k^{n}}

.. math::

    \mu_{n}^{\prime}=\left.M^{\left(n\right)}\left(t\right)\right|_{t=0}=\left.\frac{\textrm{Li}_{\alpha-n}\left(e^{t}\right)}{\zeta\left(a\right)}\right|_{t=0}=\frac{\zeta\left(\alpha-n\right)}{\zeta\left(\alpha\right)}

Implementation: `scipy.stats.zipf`
