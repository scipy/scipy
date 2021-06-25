.. _non-uniform-random-number-sampling:

=====================================================
Universal Non-Uniform Random Number Sampling in SciPy
=====================================================

.. currentmodule:: scipy.stats

SciPy provides an interface to many universal non-uniform random number
generators to sample random variates from a wide variety of univariate
continuous and discrete distributions. Implementations of a fast C library
called `UNU.RAN <http://statmath.wu.ac.at/software/unuran/>`__ are used
for speed and performance. Please look at
`UNU.RAN's documentation <http://statmath.wu.ac.at/software/unuran/doc/unuran.html>`__
for an in-depth explanation of these methods. It is heavily referred
for writing this tutorial and the documentation of all the generators.


Introduction
------------

Random variate generation is the small field of research that deals with
algorithms to generate random variates from various distributions. It is
common to assume that a uniform random number generator is available.
This is a program that produces a sequence of independent and identically
distributed continuous U(0,1) random variates (i.e. uniform random variates
on the interval (0,1)). Of course, real-world computers can never generate
ideal random numbers and they cannot produce numbers of arbitrary precision
but state-of-the-art uniform random number generators come close to this
aim. Thus random variate generation deals with the problem of transforming
such a sequence of U(0,1) random numbers into non-uniform random variates.
These methods are universal and work in a black-box fashion.

Some methods to do that are:

* The Inversion method: When the inverse :math:`F^{-1}` of the cumulative
  distribution function is known, then random variate generation is easy.
  We just generate a uniformly U(0,1) distributed random number U and
  return :math:`X = F^{-1}(U)`. As closed form solutions for the inverse
  are rarely available, one usually needs to rely on approximations of
  the inverse. See :class:`~NumericalInverseHermite`.
* The Rejection Method: The rejection method, often called
  acceptance-rejection method, has been suggested by John von Neumann in
  1951 [1]_. It involves computing an upper bound to the PDF (also called the
  hat function) and using the inversion method to generate a random
  variate, say Y, from this bound. Then a uniform random number can be
  drawn between 0 to the value of the upper bound at Y. If this number
  is less than the PDF at Y, return the sample otherwise reject it. See
  :class:`~TransformedDensityRejection`.
* The Ratio-of-Uniforms Method: This is a type of acceptance-rejection
  method which is uses minimal bounding rectangles to construct the hat
  function. See :func:`~rvs_ratio_uniforms`.
* Inversion for Discrete Distributions: The difference compared to the
  continuous case is that :math:`F` is now a step-function. To realize
  this in a computer, a search algorithm is used, the simplest of which
  is *sequential search*. A uniform random number is generated from
  U(0, 1) and probabilities are summed until the cumulative probability
  exceeds the uniform random number. The index at which this happens is
  the required random variate and is returned.


More details on these algorithms can be found in the `appendix of the UNU.RAN
user manual <http://statmath.wu.ac.at/software/unuran/doc/unuran.html#RVG>`__.


Basic concepts of the Interface
-------------------------------

Every generator needs to be set up before one can start sampling from it.
This can be done by instantiating an object of that class. Most of the
generators take a distribution object as input which contains the implementation
of required methods like PDF, CDF, etc. In addition to the distribution
object, one can also pass parameters used to set up the generator. It is also
possible to truncate the distributions using a ``domain`` parameter.  All
generators need a stream of uniform random numbers that are transformed into
random variates of the given distribution. This is done by passing a ``seed``
parameter with a NumPy BitGenerator as the uniform random number generator.
``seed`` can either be a integer, `np.random.Generator`,
`np.random.BitGenerator`, `np.random.RandomState`, or
`np.random.SeedSequence`.

.. warning:: Use of NumPy < 1.19.0 is discouraged as it doesn't have a fast
             Cython API for generating uniform random numbers and might be
             too slow to practically rely on.

All the generators have a common ``rvs`` method that can be used to draw
samples from the given distribution.

An example of this interface is shown below:

    >>> from scipy.stats import TransformedDensityRejection
    >>> from math import exp
    >>> 
    >>> class StandardNormal:
    ...     def pdf(self, x: float) -> float:
    ...         # note that the normalization constant isn't required
    ...         return exp(-0.5 * x*x)
    ...     def dpdf(self, x: float) -> float:
    ...         return -x * exp(-0.5 * x*x)
    ... 
    >>> dist = StandardNormal()
    >>> 
    >>> urng = np.random.default_rng(0x6be00a9336fc82258fe73bfa672d3d52)
    >>> rng = TransformedDensityRejection(dist, seed=urng)

As shown in the example, we first initialize a distribution object that
contains an implementation of the methods required by the generator. In
our case, we use the :class:`~TransformedDensityRejection` method which
requires a PDF and its derivative w.r.t. x (i.e. the variate). Note that
the PDF doesn't need to be vectorized. It should accept and return floats.

.. note:: One can also pass the SciPy distributions as arguments but it can
          be slow due to validations and expensive NumPy operations.
          Moreover, it doesn't always have all the information required
          by some generators like derivative of PDF for TDR method. Also,
          most of the distributions in SciPy provide a ``rvs`` method which
          can be used instead.

In the above example, we have set up an object of the
:class:`~TransformedDensityRejection` method to sample from a
standard normal distribution. Now, we can start sampling from our
distribution by calling the ``rvs`` method:

    >>> rng.rvs()
    -0.6853213192313732
    >>> rng.rvs((5, 3))
    array([[ 0.24129542,  1.28709556,  0.43244247],
           [ 0.88977668,  1.22981663, -0.09820373],
           [ 1.38588147, -1.06629264,  0.54475257],
           [-0.55153436, -0.67281406,  0.01830143],
           [ 1.37709649, -1.33036774,  0.62131967]])

.. note:: Please note the difference between the ``rvs`` method of the
          distributions present in :mod:`scipy.stats` and the one provided
          by these generators: The general aim of these generators is to
          provide a fast generation method in :class:`~rv_continuous`. The
          tools presented here allow the user to experiment with different
          methods that can lead to performance improvements in certain
          situations, e.g., if a large number of samples are needed or if
          particular shape parameters are used. Also, even if the same URNG
          (``seed``) is used, the resulting rvs will be different in general.

We can pass a ``domain`` parameter to truncate the distribution:

    >>> rng = TransformedDensityRejection(dist, domain=(-1, 1), seed=urng)
    >>> rng.rvs((5, 3))
    array([[ 0.89612408,  0.03455796,  0.02369241],
           [ 0.60741875, -0.86184683,  0.49225344],
           [-0.18349351,  0.55651665, -0.65822471],
           [ 0.53502455, -0.91848951,  0.72540482],
           [-0.63404952,  0.28427289,  0.98542246]])

Invalid and bad arguments are handled either by SciPy or by UNU.RAN. The
latter throws a ``RuntimeError`` that follows a common format:

``RuntimeError: [objid: <object id>] <error code>: <reason> => <type of error>``

where:

* ``<object id>`` is the ID of the object given by UNU.RAN
* ``<error code>`` is an error code representing a type of error.
* ``<reason>`` is the reason why the error occurred.
* ``<type of error>`` is a short description of the type of error.

The ``<reason>`` shows what caused the error. This, by itself, should contain
enough information to help debug the error. In addition, ``<error id>`` and
``<type of error>`` can be used to investigate different classes of error in
UNU.RAN. A complete list of all the error codes and their descriptions can be
found in the `Section 8.4 of the UNU.RAN user manual
<http://statmath.wu.ac.at/software/unuran/doc/unuran.html#Errno>`__.

An example of an error generated by UNU.RAN is shown below:

``RuntimeError: [objid: TDR.003] 50 : PDF(x) < 0.! => (generator) (possible) invalid data``

This shows that UNU.RAN failed to initialize an object with ID ``TDR.003``
because the PDF was < 0. i.e. negative. This falls under the type
"possible invalid data for the generator" and has error code 50.

Warnings thrown by UNU.RAN also follow the same format.


Generators in :mod:`scipy.stats`
--------------------------------
.. toctree::
   :maxdepth: 1

   sampling_tdr
   sampling_dau


References
~~~~~~~~~~

.. [1] Von Neumann, John. "13. various techniques used in connection with
       random digits." Appl. Math Ser 12.36-38 (1951): 3.
