
.. _sampling-nrou:

Naive Ratio-Of-Uniforms (NROU) Method
=====================================

.. currentmodule:: scipy.stats

* Required: PDF
* Optional: mode, center, bounding rectangle for acceptance region
* Speed:

  * Set-up: fast (if the bounding rectangle is known), otherwise slow
  * Sampling: moderate

NROU is an implementation of the (generalized) ratio-of-uniforms method which
uses (minimal) bounding rectangles ([3]_, [5]_). It uses a positive control
parameter r for adjusting the algorithm to the given distribution to improve
performance and/or to make this method applicable. Larger values of r increase
the class of distributions for which the method works at the expense of a higher
rejection constant. For computational reasons r=1 should be used if possible.
Moreover, this implementation uses the center of the distribution. 

If ``r==1``, the method works as follows: for a given pdf and a constant
``center``, define the set
``A = {(u, v) : 0 < v <= sqrt(pdf(u/v + center))}``.
If (U, V) is a random vector uniformly distributed over A,
then ``U/V + center`` follows a distribution according to the pdf.

Typical choices of ``center`` are zero or the mode of the pdf. The set A
is a subset of the rectangle ``R = [u_min, u_max] x [0, v_max]`` where

* ``v_max = sup_x sqrt(pdf(x))``
* ``u_min = inf_x (x - center) sqrt(pdf(x))``
* ``u_max = sup_x (x - center) sqrt(pdf(x))``

In particular, these values are finite if pdf is bounded and
``x**2 * pdf(x)`` is bounded (i.e. subquadratic tails).
One can generate (U, V) uniformly on R and return
``U/V + center`` if (U, V) are also in A which can be directly
verified.

The algorithm is not changed if one replaces the pdf by a function
that is proportional to it, e.g., by dropping unnecessary normalization
factors.

Intuitively, the method works well if A fills up most of the
enclosing rectangle such that the probability is high that (U, V)
lies in A whenever it lies in R as the number of required
iterations become too large otherwise. To be more precise, note that
the expected number of iterations to draw (U, V) uniformly
distributed on R such that (U, V) is also in A is given by
the ratio ``area(R) / area(A) = 2 * v_max * (u_max - u_min) / area(pdf)``,
where ``area(pdf)`` is the integral of the pdf (which is equal to one if the
probability density function is used but can take on other values if a
function proportional to the density is used). The equality holds since
the area of A is equal to ``0.5 * area(pdf)`` (Theorem 7.1 in [4]_).

In the general case, the acceptance region is defined as
``A = {(u, v) : 0 < v <= (pdf(u/v**r + center))**(1/(r+1))}`` and the
bounding rectangle is given by

* ``v_max = sup_x pdf(x)**(1/(r+1))``
* ``u_min = inf_x (x - center) (pdf(x))**(r/(r+1))``
* ``u_max = sup_x (x - center) (pdf(x))**(r/(r+1))``

Note that the orginal method is a special case with ``r=1``. The sampling
then works as follows:

* Generate ``(U,V)`` uniformly in ``[u_min, u_max] x [0, v_max]``
* Set ``X = U/V**r + center``
* If ``V**(r+1) <= pdf(X)``, return ``X``
* Else try again

An example of using this method is shown below:

    >>> import numpy as np
    >>> from scipy.stats import NaiveRatioUniforms
    >>>
    >>> class StandardNormal:
    ...     def pdf(self, x):
    ...         # note that the normalization constant is not required
    ...         return np.exp(-0.5 * x*x) 
    >>> dist = StandardNormal()
    >>> u_bound = np.sqrt(dist.pdf(np.sqrt(2))) * np.sqrt(2)
    >>> urng = np.random.default_rng()
    >>> rng = NaiveRatioUniforms(dist, v_max=1, u_min=-u_bound,
    ...                          u_max=u_bound, random_state=urng)
    >>> rng.rvs()
    -0.5677819961248616

In the above example, we have used the NROU method to sample from the standard
normal distribution. Note that we dropped the normalization constant while
computing the PDF. This usually helps speed up the sampling stage. Also, note
that the PDF doesn't need to be vectorized. It should accept and return a
scalar.

Note that ``v_max`` is determined by the mode of the distribution. If the mode is
far away from zero, note that ``u_max`` will become large if ``center == 0``: it
is at least ``mode * pdf(mode)**(r/(r+1)``. This can be avoided by shifting the
distribution, e.g., if ``center == mode``, see [5]_. We illustrate this with the
Gamma distribution. For a given shape parameter ``p > 0``, the density is proportional
to ``x**(p-1) * exp(-x)``:

    >>> import math
    >>>
    >>> class Gamma:
    ...     def __init__(self, p):
    ...         self.p = p
    ...
    ...     def pdf(self, x):
    ...         if x < 0:
    ...             return 0.0
    ...         return x**(self.p - 1) * math.exp(-x)
    ...     @staticmethod
    ...     def support():
    ...         return 0, np.inf

If ``p < 1``, the pdf is unbounded and we cannot use NROU (``v_max`` needs to
be finite). If ``p >= 1`, the mode is at ``p-1``. We want to apply NROU with
``r=1``. It is possible to compute the points that determine ``u_min`` and
``u_max`` explicitly by finding the extrema of the following function:

    >>> def u_bound(x, p, center):
    ...     if x < 0:
    ...         return 0
    ...     return (x - center) * x**((p-1)/2) * math.exp(-x/2)

Depending on ``p`` and ``center``, the bounding rectangle is given by

    >>> def rectangle(p, center):
    ...     h = (p+1+center)/2
    ...     k = np.sqrt(h**2 - center*(p-1))
    ...     u_min, u_max = u_bound(h-k, p, center), u_bound(h+k, p, center)
    ...     v_max = math.sqrt((p-1)**(p-1) * math.exp(-(p-1)))
    ...     return u_min, u_max, v_max

We can compare the resulting rejection constants if we set ``center=0`` and
``center=p-1`` (the mode) for different values of ``p``. We use the formula
for the rejection constant stated above (note that we need to take into
account the normalization constant of the Gamma density):

    >>> from scipy import special as sc
    >>> ps = np.arange(1.0, 2.5, 0.1)
    >>> reject_const, reject_const_shift = [], []
    >>> for p in ps:
    ...     # no shift (center=0)
    ...     u_min, u_max, v_max = rectangle(p, 0)
    ...     reject_const.append(2*v_max*(u_max - u_min) / sc.gamma(p))
    ...     # mode shift (center=p-1)
    ...     u_min, u_max, v_max = rectangle(p, p-1)
    ...     reject_const_shift.append(2*v_max*(u_max - u_min) / sc.gamma(p))

We can see that it becomes advantageous to apply the mode shift once as the
mode moves further away from zero:

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(ps, reject_const, 'o', label='center=0')
    >>> ax.plot(ps, reject_const_shift, 'o', label='center=mode')
    >>> plt.xlabel('Shape parameter p of the Gamma distribution')
    >>> plt.title('NROU rejection constants for the Gamma distribution')
    >>> plt.legend()
    >>> plt.show()

.. plot:: tutorial/stats/plots/nrou_plot1.py
   :align: center
   :include-source: 0

To conclude this example, we generate random variates for ``Gamma(2.2)`` by
applying NROU with mode shift and create a histogram:

    >>> from scipy import stats
    >>> p = 2.2
    >>> dist = Gamma(p)
    >>> u_min, u_max, v_max = rectangle(p, p-1)
    >>> rng = NaiveRatioUniforms(dist, center=p-1, v_max=v_max,
    ...                          u_min=u_min, u_max=u_max, random_state=urng)
    >>> rvs = rng.rvs(1000)
    >>> x = np.linspace(rvs.min()-0.1, rvs.max()+0.1, num=500)
    >>> fx = stats.gamma.pdf(x, p)
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, fx, "r-", label=" Gamma({}) pdf".format(p))
    >>> ax.hist(rvs, bins=30, density=True, alpha=0.8, label="rvs")
    >>> plt.xlabel("x")
    >>> plt.title("Samples drawn using NROU method with mode shift.")
    >>> plt.legend()
    >>> plt.show()

.. plot:: tutorial/stats/plots/nrou_plot2.py
   :align: center
   :include-source: 0

Finally, we remark that the parameter ``r`` can be used to transform the
density. In general, it is recommended to use ``r=1`` for computational
reasons. However, larger values of r increase the class of distributions
for which the method works. For example, if the right tail of the pdf
decreases like ``x**(-1.5)``, ``u_max`` is not finite if ``r=1``
since ``x*sqrt(x**(-1.5)) = x**(0.25)`` is unbounded as ``x`` increases.
However, if we use ``r=2``, note that ``x * x**(-1.5*r/(r+1)) = 1``,
so ``u_max`` is bounded.

Please see [1]_, [2]_, [3]_, [4]_ and [5]_ for more details on this method.


References
----------

.. [1] UNU.RAN reference manual, Section 5.3.11,
       "NROU - Naive Ratio-Of-Uniforms method",
       http://statmath.wu.ac.at/software/unuran/doc/unuran.html#NROU
.. [2] W. Hörmann, J. Leydold and G. Derflinger (2004). Automatic
       Nonuniform Random Variate Generation, Springer, Berlin.
.. [3] A.J. Kinderman and J.F. Monahan, "Computer Generation of Random
       Variables Using the Ratio of Uniform Deviates",
       ACM Transactions on Mathematical Software, 3(3), p. 257--260, 1977.
.. [4] L. Devroye, "Non-Uniform Random Variate Generation",
       Springer-Verlag, 1986.
.. [5] J. C. Wakefield, A. E. Gelfand, and A. F. M. Smith. "Efficient
       generation of random variates via the ratio-of-uniforms method."
       Statistics and Computing 1.2, p. 129--133, 1991.

