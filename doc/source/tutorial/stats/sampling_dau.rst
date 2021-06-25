
.. _sampling-dau:

Discrete Alias Urn (DAU)
========================

.. currentmodule:: scipy.stats

* Required: probability vector (PV) or the PMF along with a finite domain
* Speed:

    * Set-up: slow (linear with the vector-length)
    * Sampling: very fast 


DAU samples from distributions with arbitrary but finite probability vectors
(PV) of length N. The algorithm is based on an ingenious method by A.J.
Walker and requires a table of size (at least) N. It needs one random number
and only one comparison for each generated random variate. The setup time for
constructing the tables is O(N).

    >>> import numpy as np
    >>> from scipy.stats import DiscreteAliasUrn
    >>> 
    >>> pv = [0.18, 0.02, 0.8]
    >>> urng = np.random.default_rng(0x8383b295ecbf0874ac910b13adcc85a0)
    >>> rng = DiscreteAliasUrn(pv, seed=urng)
    >>> rng.rvs()
    2

By default, the probability vector is indexed starting at 0. However, this
can be changed by passing a ``domain`` parameter. When ``domain`` is given
in combination with the PV, it has the effect of relocating the
distribution from ``(0, N)`` to ``(domain[0]``, ``domain[0] + N)``.
``domain[1]`` is ignored in this case.

   >>> rng = DiscreteAliasUrn(pv, domain=(10, 13), seed=urng)
   >>> rng.rvs()
   12

The method also works when no probability vector but a PMF is given. However
then additionally a bounded (finite) domain must also be given.

    >>> class Distribution:
    ...     def pmf(self, x, c):
    ...         return x**c
    ... 
    >>> dist = Distribution()
    >>> rng = DiscreteAliasUrn(dist=dist, domain=(1, 10), params=(2, ),
    ...                        seed=urng)
    >>> rng.rvs()
    9

.. note:: As :class:`~DiscreteAliasUrn` expects PMF with signature
          ``def pmf(x: float, ...) -> float``, it first vectorizes the
          PMF using ``np.vectorize`` and then evaluates it over all the
          points in the domain. But if the PMF is already vectorized,
          it is much faster to just evaluate it at each point in the domain
          and pass the obtained PV instead along with the domain.
          For example, ``pmf`` methods of SciPy's discrete distributions
          are vectorized and a PV can be obtained by doing:

          >>> from scipy.stats import binom, DiscreteAliasUrn
          >>> dist = binom(10, 0.2)  # distribution object
          >>> domain = dist.support()  # the domain of your distribution
          >>> x = np.arange(domain[0], domain[1] + 1)
          >>> pv = dist.pmf(x)
          >>> rng = DiscreteAliasUrn(pv, domain=domain)

          Domain is required here to relocate the distribution.

The performance can be slightly influenced by setting the size of the used
table which can be changed by passing a ``urn_factor`` parameter.

    >>> # use a table twice the length of PV.
    >>> urn_factor = 2
    >>> rng = DiscreteAliasUrn(pv, urn_factor=urn_factor, seed=urng)
    >>> rng.rvs()
    2

.. note:: It is recommended to keep this parameter under 2.


References
----------
.. [1] UNU.RAN reference manual, Section 5.8.2,
       "DAU - (Discrete) Alias-Urn method",
       http://statmath.wu.ac.at/software/unuran/doc/unuran.html#DAU
.. [2] A.J. Walker (1977). An efficient method for generating discrete
       random variables with general distributions, ACM Trans. Math.
       Software 3, pp. 253-256.
