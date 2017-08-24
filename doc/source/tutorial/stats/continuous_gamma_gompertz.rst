
.. _continuous-gamma_gompertz:

Gamma/Gompertz Distribution
========================================

For :math:`x\geq0`, :math:`c>0`, and :math:`\beta>0` . 

.. math::

    f\left(x;c,\beta \right) = \frac{c {e^x} {{\beta }^c} }{ {\left( \beta - 1 + {e ^ x} \right) }^ {c+1} }
    
If :math:`\beta=1,` the distribution reduces to the exponential distribution.

Implementation: `scipy.stats.gamma_gompertz`
