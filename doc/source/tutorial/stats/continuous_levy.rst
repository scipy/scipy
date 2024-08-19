
.. _continuous-levy:

Lévy Distribution
==================

A special case of Lévy-stable distributions with :math:`\alpha=\frac{1}{2}`
and :math:`\beta=1` and support :math:`x\geq0`.  In standard form it is defined for :math:`x>0` as

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{x\sqrt{2\pi x}}\exp\left(-\frac{1}{2x}\right)\\ F\left(x\right) & = & 2\left[1-\Phi\left(\frac{1}{\sqrt{x}}\right)\right]\\ G\left(q\right) & = & \left[\Phi^{-1}\left(1-\frac{q}{2}\right)\right]^{-2}.\end{eqnarray*}

It has no finite moments.

Implementation: `scipy.stats.levy`
