
.. _continuous-landau:

Landau distribution
===================

A special case of Lévy-stable distributions with :math:`\alpha=1`
and :math:`\beta=1` and support :math:`-\infty < x < \infty`. The probability
density function is given by

.. math::
   :nowrap:

    f(x) = \frac{1}{2\pi i}\int_{a-i\infty}^{a+i\infty}e^{s\log(s)+xs}\ ds

Implementation: `scipy.stats.landau`
