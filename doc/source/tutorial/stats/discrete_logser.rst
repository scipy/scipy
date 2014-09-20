
.. _discrete-logser:

Logarithmic (Log-Series, Series) Distribution
=============================================

The logarimthic distribution with parameter :math:`p` has a probability mass function with terms proportional to the Taylor
series expansion of :math:`\log\left(1-p\right)`

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k;p\right) & = & -\frac{p^{k}}{k\log\left(1-p\right)}\quad k\geq1\\ F\left(x;p\right) & = & -\frac{1}{\log\left(1-p\right)}\sum_{k=1}^{\left\lfloor x\right\rfloor }\frac{p^{k}}{k}=1+\frac{p^{1+\left\lfloor x\right\rfloor }\Phi\left(p,1,1+\left\lfloor x\right\rfloor \right)}{\log\left(1-p\right)}\end{eqnarray*}

where

.. math::
   :nowrap:

    \[ \Phi\left(z,s,a\right)=\sum_{k=0}^{\infty}\frac{z^{k}}{\left(a+k\right)^{s}}\]

is the Lerch Transcendent. Also define :math:`r=\log\left(1-p\right)`

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & -\frac{p}{\left(1-p\right)r}\\ \mu_{2} & = & -\frac{p\left[p+r\right]}{\left(1-p\right)^{2}r^{2}}\\ \gamma_{1} & = & -\frac{2p^{2}+3pr+\left(1+p\right)r^{2}}{r\left(p+r\right)\sqrt{-p\left(p+r\right)}}r\\ \gamma_{2} & = & -\frac{6p^{3}+12p^{2}r+p\left(4p+7\right)r^{2}+\left(p^{2}+4p+1\right)r^{3}}{p\left(p+r\right)^{2}}.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} M\left(t\right) & = & -\frac{1}{\log\left(1-p\right)}\sum_{k=1}^{\infty}\frac{e^{tk}p^{k}}{k}\\  & = & \frac{\log\left(1-pe^{t}\right)}{\log\left(1-p\right)}\end{eqnarray*}

Thus,

.. math::
   :nowrap:

    \[ \mu_{n}^{\prime}=\left.M^{\left(n\right)}\left(t\right)\right|_{t=0}=\left.\frac{\textrm{Li}_{1-n}\left(pe^{t}\right)}{\log\left(1-p\right)}\right|_{t=0}=-\frac{\textrm{Li}_{1-n}\left(p\right)}{\log\left(1-p\right)}.\]

Implementation: `scipy.stats.logser`
