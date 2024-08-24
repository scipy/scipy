
.. _continuous-landau:

Landau distribution
===================

A special case of Lévy-stable distributions with :math:`\alpha=1`
and :math:`\beta=1` and support :math:`-\infty < x < \infty`. The probability
density function is given by

.. math::

    f(x) = \frac{1}{\pi}\int_0^\infty \exp(-t \log t - xt)\sin(\pi t) dt

The differential entropy is 2.37263644000448182, and the moments are undefined.

Implementation: `scipy.stats.landau`
