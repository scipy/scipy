# The portion below involving multigammaln has the following copyright
# information:
# Last Change: Sat Mar 21 02:00 PM 2009 J

# Copyright (c) 2001, 2002 Enthought, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   a. Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#   b. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#   c. Neither the name of the Enthought nor the names of its contributors
#      may be used to endorse or promote products derived from this software
#      without specific prior written permission.
#
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

# The Poisson binomial distribution functions copyright SciPy developers 2026.

"""Some more special functions which may be useful for multivariate statistical
analysis."""

import numpy as np

import scipy.special._gufuncs as _gufuncs
from scipy.special import gammaln as loggam
from scipy.special._ufunc_tools import _with_cache_optimization


__all__ = ['multigammaln', 'poisson_binom_cdf']


def multigammaln(a, d):
    r"""Returns the log of multivariate gamma, also sometimes called the
    generalized gamma.

    Parameters
    ----------
    a : ndarray
        The multivariate gamma is computed for each item of `a`.
    d : int
        The dimension of the space of integration.

    Returns
    -------
    res : ndarray
        The values of the log multivariate gamma at the given points `a`.

    Notes
    -----
    The formal definition of the multivariate gamma of dimension d for a real
    `a` is

    .. math::

        \Gamma_d(a) = \int_{A>0} e^{-tr(A)} |A|^{a - (d+1)/2} dA

    with the condition :math:`a > (d-1)/2`, and :math:`A > 0` being the set of
    all the positive definite matrices of dimension `d`.  Note that `a` is a
    scalar: the integrand only is multivariate, the argument is not (the
    function is defined over a subset of the real set).

    This can be proven to be equal to the much friendlier equation

    .. math::

        \Gamma_d(a) = \pi^{d(d-1)/4} \prod_{i=1}^{d} \Gamma(a - (i-1)/2).

    References
    ----------
    R. J. Muirhead, Aspects of multivariate statistical theory (Wiley Series in
    probability and mathematical statistics).

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import multigammaln, gammaln
    >>> a = 23.5
    >>> d = 10
    >>> multigammaln(a, d)
    454.1488605074416

    Verify that the result agrees with the logarithm of the equation
    shown above:

    >>> d*(d-1)/4*np.log(np.pi) + gammaln(a - 0.5*np.arange(0, d)).sum()
    454.1488605074416
    """
    a = np.asarray(a)
    # Support for 0d arrays is needed for array_api_strict and dask.
    d = np.asarray(d)[()]
    if not np.isscalar(d) or (np.floor(d) != d):
        raise ValueError("d should be a positive integer (dimension)")
    if np.any(a <= 0.5 * (d - 1)):
        raise ValueError(f"condition a ({a}) > 0.5 * (d-1) ({0.5 * (d-1)}) not met")

    res = (d * (d-1) * 0.25) * np.log(np.pi)
    res += np.sum(loggam([(a - (j - 1.)/2) for j in range(1, d+1)]), axis=0)
    return res


_poisson_binom_pmf_doc = (
    r"""Returns pmf of Poisson Binomial distribution.

    Parameters
    ----------
    k : array
        Number of successes at which to evaluate pmf.

    p : array
        Success probabilities of independent Bernoulli trials.


    """
)


_poisson_binom_pmf = _with_cache_optimization(
    name="_poisson_binom_pmf",
    arg_names=["k", "p"],
    docstring=_poisson_binom_pmf_doc,
    ufunc=_gufuncs._poisson_binom_pmf,
    cache_arg_indices=[1],
)


_poisson_binom_cdf_doc = (
    r"""Poisson binomial cumulative distribution function.

    The Poisson binomial distribution is the discrete probability
    distribution of a sum of independent Bernoulli trials that are not
    necessarily identically distributed [1]_.

    .. math::

        \sum_{l=0}^{k}\sum_{A \in F_l}\prod_{i \in A} p_i \prod_{j \in A^c} (1 - p_j)

    where :math:`F_l` is the set of all subsets of size :math:`l` of
    :math:`\{1, 2, \ldots, n\}` and :math:`A^c` is the complement of :math:`A` in
    :math:`\{1, 2, \ldots, n\}`.

    Parameters
    ----------
    k : array_like
        Number of successes at which to evaluate cdf.
    p : array_like
        Success probabilities of independent Bernoulli trials.
    out : ndarray, optional
        Optional output array for the function results.
    **kwargs
        For other keyword-only arguments, see the
        `NumPy ufunc docs <https://numpy.org/doc/stable/reference/ufuncs.html#optional-keyword-arguments>`_.

    Returns
    -------
    scalar or ndarray
        Values of the Poisson binomial cumulative distribution function.

    References
    ----------
    .. [1] Poisson binomial distribution,
           https://en.wikipedia.org/wiki/Poisson_binomial_distribution


    """
)

poisson_binom_cdf = _with_cache_optimization(
    name="poisson_binom_cdf",
    arg_names=["k", "p"],
    docstring=_poisson_binom_cdf_doc,
    ufunc=_gufuncs._poisson_binom_cdf,
    cache_arg_indices=[1],
    module="scipy.special._spfun_stats",
)
