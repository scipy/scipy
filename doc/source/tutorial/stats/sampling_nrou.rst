
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
uses (minimal) bounding rectangles. It uses a positive control parameter r for
adjusting the algorithm to the given distribution to improve performance and/or
to make this method applicable. Larger values of r increase the class of
distributions for which the method works at the expense of a higher
rejection constant. For computational reasons r=1 should be used if possible.
Moreover, this implementation uses the center of the distribution. 

If ``r==1``, the method works as follows: for a given ``pdf`` and a constant
`center`, define the set
``A = {(u, v) : 0 < v <= sqrt(pdf(u/v + center))}``.
If ``(U, V)`` is a random vector uniformly distributed over ``A``,
then ``U/V + center`` follows a distribution according to ``pdf``.

Typical choices of `center` are zero or the mode of `pdf`. The set ``A``
is a subset of the rectangle ``R = [u_min, u_max] x [0, v_max]`` where

* ``v_max = sup_x sqrt(pdf(x))``
* ``u_min = inf_x (x - center) sqrt(pdf(x))``
* ``u_max = sup_x (x - center) sqrt(pdf(x))``

In particular, these values are finite if ``pdf`` is bounded and
``x**2 * pdf(x)`` is bounded (i.e. subquadratic tails).
One can generate ``(U, V)`` uniformly on ``R`` and return
``U/V + center`` if ``(U, V)`` are also in ``A`` which can be directly
verified.

The algorithm is not changed if one replaces ``pdf`` by ``k * pdf`` for any
constant k > 0. Thus, it is often convenient to work with a function
that is proportional to the probability density function by dropping
unneccessary normalization factors.

Intuitively, the method works well if ``A`` fills up most of the
enclosing rectangle such that the probability is high that ``(U, V)``
lies in ``A`` whenever it lies in ``R`` as the number of required
iterations becomes too large otherwise. To be more precise, note that
the expected number of iterations to draw ``(U, V)`` uniformly
distributed on ``R`` such that ``(U, V)`` is also in ``A`` is given by
the ratio ``area(R) / area(A) = 2 * v_max * (u_max - u_min) / area(pdf)``,
where ``area(pdf)`` is the integral of ``pdf`` (which is equal to one if the
probability density function is used but can take on other values if a
function proportional to the density is used). The equality holds since
the area of ``A`` is equal to ``0.5 * area(pdf)`` (Theorem 7.1 in [4]_).

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
                                 u_max=u_bound, random_state=urng)
    >>> rng.rvs()
    -0.5677819961248616

In the above example, we have used the NROU method to sample from the standard
normal distribution. Note that we dropped the normalization constant while
computing the PDF. This usually helps speed up the sampling stage. Also, note
that the PDF doesn't need to be vectorized. It should accept and return a
scalar.


TODO: show mode shift, r != 0

Please see [1]_, [2]_, [3]_ and [4]_ for more details on this method.


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

