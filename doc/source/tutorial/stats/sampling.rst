.. _non-uniform-random-number-sampling:

=====================================================
Universal Non-Uniform Random Number Sampling in SciPy
=====================================================

.. currentmodule:: scipy.stats

SciPy provides an interface to many universal non-uniform random number
samplers to sample random variates from a wide variety of univariate
continuous and discrete distributions. Implementations of a fast C library
called `UNU.RAN <http://statmath.wu.ac.at/software/unuran/>`__ are used
for speed and performance. Please look at
`UNU.RAN's documentation <http://statmath.wu.ac.at/software/unuran/doc/unuran.html>`__
for an in-depth explanation of these methods. It is heavily referred
for writing this tutorial and the documentation of all the samplers.


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

Some methods to do that are:

* The Inversion method: When the inverse :math:`F^{-1}` of the cumulative
  distribution function is known, then random variate generation is easy.
  We just generate a uniformly U(0,1) distributed random number U and
  return :math:`X = F^{-1}(U)`. See
  :class:`~NumericalInverseHermite`.
* The Rejection Method: The rejection method, often called
  acceptance-rejection method, has been suggested by John von Neumann in
  1951. It involves computing an upper bound to the PDF (also called the
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

Every sampler needs to be set up before one can start sampling from it.
This can be done by instantiating an object of that class. Most of the
samplers take a distribution object as input which contains the implementation
of required methods like PDF, CDF, etc. In addition to the distribution
object, one can also pass parameters used to set up the sampler. It is also
possible to truncate the distributions using a ``domain`` parameter.  All
samplers need a stream of uniform random numbers that are transformed into
random variates of the given distribution. This is done by passing a ``seed``
parameter with a NumPy BitGenerator as the uniform random number generator.
``seed`` can either be a integer, `np.random.Generator`,
`np.random.BitGenerator`, `np.random.RandomState`, or
`np.random.SeedSequence`.

.. warning:: Use of NumPy < 1.19.0 is discouraged as it doesn't have a fast
             Cython API for generating uniform random numbers and might be
             too slow to practically rely on.

All the samplers have a common ``rvs`` method that can be used to draw
samples from the given distribution.

An example of this interface is shown below:

    >>> import numpy as np
    >>> from scipy.stats import TransformedDensityRejection
    >>> 
    >>> class StandardNormal:
    ...     def pdf(self, x):
    ...         # note that the normalization contant isn't required
    ...         return np.exp(-0.5 * x*x)
    ...     def dpdf(self, x):
    ...         return -x * np.exp(-0.5 * x*x)
    ... 
    >>> dist = StandardNormal()
    >>> 
    >>> urng = np.random.default_rng(0x6be00a9336fc82258fe73bfa672d3d52)
    >>> rng = TransformedDensityRejection(dist, seed=urng)

.. note:: One can also pass the SciPy distributions as arguments but it can
          be slow due to validations and expensive NumPy operations.
          Moreover, it doesn't always have all the information required
          by some samplers like derivative of PDF for the TDR sampler. Also,
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


Transformed Density Rejection (TDR)
-----------------------------------

* Required: T-concave PDF, dPDF
* Optional: mode, center
* Speed:

  * Set-up: slow
  * Sampling: fast 


TDR is an acceptance/rejection method that uses the concavity of a transformed
density to construct hat function and squeezes automatically. Such PDFs are
called T-concave. Currently the following transformations are implemented:

.. math::

    c = 0.: T(x) &= \log(x)\\
    c = -0.5: T(x) &= \frac{1}{\sqrt{x}} \text{ (Default)}

In addition to the PDF, it also requires the derivative of the PDF w.r.t x
(i.e. the variate). These functions must be present as methods of a python
object which can then be passed to the samplers to instantiate their object.

Three variants of this method are available:

* GW: squeezes between construction points.
* PS: squeezes proportional to hat function. (Default)
* IA: same as variant PS but uses a composition method with
  "immediate acceptance" in the region below the squeeze.

This can be changed by passing a ``variant`` parameter.

An example of using this method is shown below:

    >>> from scipy.stats import TransformedDensityRejection
    >>> 
    >>> class Distribution:
    ...     def pdf(self, x):
    ...         if abs(x) <= 1:
    ...             return 1 - x*x
    ...         return 0
    ...     def dpdf(self, x):
    ...         if abs(x) <= 1:
    ...             return -2*x
    ...         return 0
    ... 
    >>> dist = Distribution()
    >>> 
    >>> urng = np.random.default_rng(0xbb51a838c8a087854ad3643b2b268d43)
    >>> rng = TransformedDensityRejection(dist, seed=urng)
    >>> rng.rvs()
    0.39606414669189344

It is also possible to evaluate the inverse of the CDF of the hat distribution
directly using the ``ppf_hat`` method.

    >>> rng.ppf_hat(0.5)
    -5.7844272705533106e-05
    >>> u = np.linspace(0, 1, num=10)
    >>> rng.ppf_hat(u)
    array([       -inf, -0.58619335, -0.39060804, -0.22634331, -0.07436531,
            0.07424874,  0.22622123,  0.39047127,  0.58601716,         inf])

For densities with modes not close to 0 it is suggested to set either the
mode or the center of the distribution by passing ``mode`` or ``center``
parameters. The latter is the approximate location of the mode or the mean
of the distribution. This location provides some information about the main
part of the PDF and is used to avoid numerical problems.

    >>> # mode = 0 for our distribution
    >>> # if exact mode is not available, pass 'center' parameter instead
    >>> rng = TransformedDensityRejection(dist, mode=0.)

By default, the method uses 30 construction points to construct the hat.
This can be changed by passing a ``cpoints`` parameter which can either
be an array of construction points or an integer representing the number
of construction points to use.

    >>> rng = TransformedDensityRejection(dist, cpoints=[-0.5, 0., 0.5])

This method accepts many other set-up parameters. See the documentation for
an exclusive list. More information of the parameters and the method can be
found in `Section 5.3.16 of the UNU.RAN user manual
<http://statmath.wu.ac.at/software/unuran/doc/unuran.html#TDR>`__.


Discrete Alias Urn (DAU)
------------------------

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
