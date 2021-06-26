
.. _sampling-tdr:

Transformed Density Rejection (TDR)
===================================

.. currentmodule:: scipy.stats

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
object which can then be passed to the generators to instantiate their object.

Three variants of this method are available:

* GW: squeezes between construction points.
* PS: squeezes proportional to hat function. (Default)
* IA: same as variant PS but uses a composition method with
  "immediate acceptance" in the region below the squeeze.

This can be changed by passing a ``variant`` parameter.

An example of using this method is shown below:

    >>> from scipy.stats import TransformedDensityRejection
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
    >>> urng = np.random.default_rng(0xbb51a838c8a087854ad3643b2b268d43)
    >>> rng = TransformedDensityRejection(dist, seed=urng)
    >>> rng.rvs()
    0.7771420061983989

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

* 'sqhratio': (area below squeeze) / (area below hat) for the generator. It
  is a number between 0 and 1. Closer to 1 means that the hat and the squeeze
  functions tightly envelop the distribution and fewer PDF evaluations are
  required to generate samples. The expected number of evaluations of the
  density is bounded by ``(1/sqhratio) - 1`` per sample. By default, it is
  kept above 0.99 but that can be changed by passing a ``max_sqhratio``
  parameter.
* 'hat_area': area below the hat for the generator.
* 'squeeze_area': area below the squeeze for the generator.

    >>> rng.sqhratio
    0.9947024204884917
    >>> rng.hat_area
    2.510253139791547
    >>> rng.squeeze_area
    2.4969548741894876
    >>> rng.sqhratio == rng.squeeze_area / rng.hat_area
    True

To increase ``sqhratio``, pass ``max_sqhratio``:

    >>> rng = TransformedDensityRejection(dist, max_sqhratio=0.999,
    ...                                   max_intervals=1000, seed=urng)
    >>> rng.sqhratio
    0.999364900465214

Note that we need to increase the ``max_intervals`` parameter when we want
a higher ``sqhratio``. This is because more construction points are required
to fit the distribution more tightly.

Let's see how this affects the callbacks to the PDF method of the
distribution:

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
    >>> urng = np.random.default_rng(0x7344aca86f1dafd8820be998b639551c)
    >>> rng = TransformedDensityRejection(dist1, seed=urng)
    >>> dist1.callbacks  # evaluations during setup
    139
    >>> dist1.callbacks = 0  # don't consider evaluations during setup
    >>> rvs = rng.rvs(100000)
    >>> dist1.callbacks  # evaluations during sampling
    551
    >>> dist2 = StandardNormal()
    >>> # use the same stream of uniform random numbers
    >>> urng = np.random.default_rng(0x7344aca86f1dafd8820be998b639551c)
    >>> rng = TransformedDensityRejection(dist2, max_sqhratio=0.999,
    ...                                   max_intervals=1000, seed=urng)
    >>> dist2.callbacks  # evaluations during setup
    467
    >>> dist2.callbacks = 0  # don't consider evaluations during setup
    >>> rvs = rng.rvs(100000)
    >>> dist2.callbacks  # evaluations during sampling
    63

As we can see, far fewer PDF evaluations are required during sampling when
we increase the ``sqhratio``. Also, notice that this comes at the cost of
increased PDF evaluations during setup.

For densities with modes not close to 0, it is suggested to set either the
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

    >>> rng = TransformedDensityRejection(dist, cpoints=[-5, 0, 5])

It is also possible to change the variant of the method used by passing
a ``variant`` parameter:

    >>> rng = TransformedDensityRejection(dist, variant='ia')

This method accepts many other set-up parameters. See the documentation for
an exclusive list. More information of the parameters and the method can be
found in `Section 5.3.16 of the UNU.RAN user manual
<http://statmath.wu.ac.at/software/unuran/doc/unuran.html#TDR>`__.


References
----------

.. [1] UNU.RAN reference manual, Section 5.3.16,
       "TDR - Transformed Density Rejection",
       http://statmath.wu.ac.at/software/unuran/doc/unuran.html#TDR
.. [2] Hörmann, Wolfgang. "A rejection technique for sampling from
       T-concave distributions." ACM Transactions on Mathematical
       Software (TOMS) 21.2 (1995): 182-193
.. [3] W.R. Gilks and P. Wild (1992). Adaptive rejection sampling for
       Gibbs sampling, Applied Statistics 41, pp. 337-348.
