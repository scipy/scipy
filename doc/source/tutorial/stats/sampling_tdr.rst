
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
