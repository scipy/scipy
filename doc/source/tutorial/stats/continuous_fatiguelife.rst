
.. _continuous-fatiguelife:

Fatigue Life (Birnbaum-Saunders) Distribution
=============================================

This distribution's pdf is the average of the inverse-Gaussian :math:`\left(\mu=1\right)` and reciprocal inverse-Gaussian pdf :math:`\left(\mu=1\right)` . We follow the notation of JKB here with :math:`\beta=S.` for :math:`x>0`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{x+1}{2c\sqrt{2\pi x^{3}}}\exp\left(-\frac{\left(x-1\right)^{2}}{2xc^{2}}\right)\\ F\left(x;c\right) & = & \Phi\left(\frac{1}{c}\left(\sqrt{x}-\frac{1}{\sqrt{x}}\right)\right)\\ G\left(q;c\right) & = & \frac{1}{4}\left[c\Phi^{-1}\left(q\right)+\sqrt{c^{2}\left(\Phi^{-1}\left(q\right)\right)^{2}+4}\right]^{2}\end{eqnarray*}

.. math::

     M\left(t\right)=c\sqrt{2\pi}\exp\left[\frac{1}{c^{2}}\left(1-\sqrt{1-2c^{2}t}\right)\right]\left(1+\frac{1}{\sqrt{1-2c^{2}t}}\right)

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{c^{2}}{2}+1\\ \mu_{2} & = & c^{2}\left(\frac{5}{4}c^{2}+1\right)\\ \gamma_{1} & = & \frac{4c\sqrt{11c^{2}+6}}{\left(5c^{2}+4\right)^{3/2}}\\ \gamma_{2} & = & \frac{6c^{2}\left(93c^{2}+41\right)}{\left(5c^{2}+4\right)^{2}}\end{eqnarray*}

Implementation: `scipy.stats.fatiguelife`
