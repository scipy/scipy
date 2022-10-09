
.. _sampling-tdr:

Transformed Density Rejection (TDR)
===================================

.. currentmodule:: scipy.stats.sampling

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

In addition to the PDF, it also requires the derivative of the PDF w.r.t ``x``
(i.e. the variate). These functions must be present as methods of a python
object which can then be passed to the generators to instantiate their object.
The variant that is implemented uses squeezes proportional to hat function ([1]_).

An example of using this method is shown below:

    >>> from scipy.stats.sampling import TransformedDensityRejection
    >>> from scipy.stats import norm
    >>> 
    >>> class StandardNormal:
    ...     def pdf(self, x):
    ...         # note that the normalization constant is not required
    ...         return np.exp(-0.5 * x*x)
    ...     def dpdf(self, x):
    ...         return -x * np.exp(-0.5 * x*x)
    ... 
    >>> dist = StandardNormal()
    >>> 
    >>> urng = np.random.default_rng()
    >>> rng = TransformedDensityRejection(dist, random_state=urng)
    >>> rng.rvs()
    -1.526829048388144

In the above example, we have used the TDR method to sample from the standard
normal distribution. Note that we can drop the normalization constant while
computing the PDF. This usually helps speed up the sampling stage. Also, note
that the PDF doesn't need to be vectorized. It should accept and return a
scalar.

It is also possible to evaluate the inverse of the CDF of the hat distribution
directly using the ``ppf_hat`` method.

    >>> rng.ppf_hat(0.5)
    -0.00018050266342362759
    >>> norm.ppf(0.5)
    0.0
    >>> u = np.linspace(0, 1, num=10)
    >>> rng.ppf_hat(u)
    array([       -inf, -1.22227372, -0.7656556 , -0.43135731, -0.14002921,
            0.13966423,  0.43096141,  0.76517113,  1.22185606,         inf])
    >>> norm.ppf(u)
    array([       -inf, -1.22064035, -0.76470967, -0.4307273 , -0.1397103 ,
            0.1397103 ,  0.4307273 ,  0.76470967,  1.22064035,         inf])

Apart from the PPF method, other attributes can be accessed
to see how well the generator fits the given distribution. These are:

* 'squeeze_hat_ratio': (area below squeeze) / (area below hat) for the generator. It
  is a number between 0 and 1. Closer to 1 means that the hat and the squeeze
  functions tightly envelop the distribution and fewer PDF evaluations are
  required to generate samples. The expected number of evaluations of the
  density is bounded by ``(1/squeeze_hat_ratio) - 1`` per sample. By default, it is
  kept above 0.99 but that can be changed by passing a ``max_squeeze_hat_ratio``
  parameter.
* 'hat_area': area below the hat for the generator.
* 'squeeze_area': area below the squeeze for the generator.

    >>> rng.squeeze_hat_ratio
    0.9947024204884917
    >>> rng.hat_area
    2.510253139791547
    >>> rng.squeeze_area
    2.4969548741894876
    >>> rng.squeeze_hat_ratio == rng.squeeze_area / rng.hat_area
    True

The distribution can be truncated by passing a domain parameter:

    >>> urng = np.random.default_rng()
    >>> rng = TransformedDensityRejection(dist, domain=[0, 1], random_state=urng)
    >>> rng.rvs(10)
    array([0.05452512, 0.97251362, 0.49955877, 0.82789729, 0.33048885,
           0.55558548, 0.23168323, 0.13423275, 0.73176575, 0.35739799])

If the domain is not specified, the ``support`` method of the ``dist`` object
is used to determine the domain:

    >>> class StandardNormal:
    ...     def pdf(self, x):
    ...         return np.exp(-0.5 * x*x)
    ...     def dpdf(self, x):
    ...         return -x * np.exp(-0.5 * x*x)
    ...     def support(self):
    ...         return -np.inf, np.inf
    ... 
    >>> dist = StandardNormal()
    >>> 
    >>> urng = np.random.default_rng()
    >>> rng = TransformedDensityRejection(dist, random_state=urng)
    >>> rng.rvs(10)
    array([-1.52682905,  2.06206883,  0.15205036,  1.11587367, -0.30775562,
           0.29879802, -0.61858268, -1.01049115,  0.78853694, -0.23060766])

If the ``dist`` object does not provide a ``support`` method, the domain
is assumed to be ``(-np.inf, np.inf)``.

To increase ``squeeze_hat_ratio``, pass ``max_squeeze_hat_ratio``:

    >>> dist = StandardNormal()
    >>> rng = TransformedDensityRejection(dist, max_squeeze_hat_ratio=0.999,
    ...                                   random_state=urng)
    >>> rng.squeeze_hat_ratio
    0.999364900465214

Let's see how this affects the callbacks to the PDF method of the
distribution:

    >>> from copy import copy
    >>> class StandardNormal:
    ...     def __init__(self):
    ...         self.callbacks = 0
    ...     def pdf(self, x):
    ...         self.callbacks += 1
    ...         return np.exp(-0.5 * x*x)
    ...     def dpdf(self, x):
    ...         return -x * np.exp(-0.5 * x*x)
    ... 
    >>> dist1 = StandardNormal()
    >>> urng1 = np.random.default_rng()
    >>> urng2 = copy(urng1)
    >>> rng1 = TransformedDensityRejection(dist1, random_state=urng1)
    >>> dist1.callbacks  # evaluations during setup
    139
    >>> dist1.callbacks = 0  # don't consider evaluations during setup
    >>> rvs = rng1.rvs(100000)
    >>> dist1.callbacks  # evaluations during sampling
    527
    >>> dist2 = StandardNormal()
    >>> # use the same stream of uniform random numbers
    >>> rng2 = TransformedDensityRejection(dist2, max_squeeze_hat_ratio=0.999,
    ...                                    random_state=urng2)
    >>> dist2.callbacks  # evaluations during setup
    467
    >>> dist2.callbacks = 0  # don't consider evaluations during setup
    >>> rvs = rng2.rvs(100000)
    >>> dist2.callbacks  # evaluations during sampling
    84

As we can see, far fewer PDF evaluations are required during sampling when
we increase the ``squeeze_hat_ratio``. The PPF-hat function is also more accurate:

    >>> abs(norm.ppf(0.975) - rng1.ppf_hat(0.975))
    0.0027054565421578136
    >>> abs(norm.ppf(0.975) - rng2.ppf_hat(0.975))
    0.00047824084476300044

Though, notice that this comes at the cost of increased PDF evaluations
during setup.

For densities with modes not close to 0, it is suggested to set either the
mode or the center of the distribution by passing ``mode`` or ``center``
parameters. The latter is the approximate location of the mode or the mean
of the distribution. This location provides some information about the main
part of the PDF and is used to avoid numerical problems.

    >>> # mode = 0 for our distribution
    >>> # if exact mode is not available, pass 'center' parameter instead
    >>> rng = TransformedDensityRejection(dist, mode=0.)

By default, the method uses 30 construction points to construct the hat.
This can be changed by passing a ``construction_points`` parameter which
can either be an array of construction points or an integer representing
the number of construction points to use.

    >>> rng = TransformedDensityRejection(dist,
    ...                                   construction_points=[-5, 0, 5])

This method accepts many other set-up parameters. See the documentation for
an exclusive list. More information of the parameters and the method can be
found in `Section 5.3.16 of the UNU.RAN user manual
<http://statmath.wu.ac.at/software/unuran/doc/unuran.html#TDR>`__.


Please see [1]_ and [2]_ for more details on this method.


References
----------

.. [1] UNU.RAN reference manual, Section 5.3.16,
       "TDR - Transformed Density Rejection",
       http://statmath.wu.ac.at/software/unuran/doc/unuran.html#TDR
.. [2] Hörmann, Wolfgang. "A rejection technique for sampling from
       T-concave distributions." ACM Transactions on Mathematical
       Software (TOMS) 21.2 (1995): 182-193

