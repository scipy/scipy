# -*- coding: utf-8 -*-
r"""
Quasi-Monte Carlo methods (:mod:`scipy.stats.qmc`)
==================================================

.. currentmodule:: scipy.stats.qmc

This module provides Quasi-Monte Carlo generators and associated helper
functions.


Quasi-Monte Carlo
=================

Engines
-------

.. autosummary::
   :toctree: generated/

   QMCEngine
   Sobol
   Halton
   OrthogonalLatinHypercube
   LatinHypercube
   OptimalDesign
   NormalQMC
   MultivariateNormalQMC

Helpers
-------

.. autosummary::
   :toctree: generated/

   discrepancy
   scale
   multinomial_qmc


Introduction to Quasi-Monte Carlo
=================================

Quasi-Monte Carlo (QMC) methods [1]_, [2]_, [3]_ provide an
:math:`n \times dim` matrix of numbers in :math:`[0,1]`. They can be used in
place of n points from the :math:`U[0,1]^{dim}` distribution. Compared to
random points, QMC points are designed to have fewer gaps and clumps. This is
quantified by discrepancy measures [4]_. From the Koksma-Hlawka
inequality [5]_ we know that low discrepancy reduces a bound on
integration error. Averaging a function :math:`f` over :math:`n` QMC points
can achieve an integration error close to :math:`O(n^{-1})` for well
behaved functions [2]_.

Most QMC constructions are designed for special values of :math:`n`
such as powers of 2 or large primes. Changing the sample
size by even one can degrade their performance, even their
rate of convergence [6]_. For instance :math:`n=100` points may give less
accuracy than :math:`n=64` if the method was designed for :math:`n=2^m`.

Some QMC constructions are extensible in n: we can find
another special sample size :math:`n' > n` and often an infinite
sequence of increasing special sample sizes. Some QMC
constructions are extensible in :math:`dim`: we can increase the dimension,
possibly to some upper bound, and typically without requiring
special values of d. Some QMC methods are extensible in
both :math:`n` and :math:`dim`.

QMC points are deterministic. That makes it hard to estimate
the accuracy of averages. Randomized QMC (RQMC) [7]_
points are constructed so that each point is individually :math:`U[0,1]^{dim}`
while collectively the n points retain their low discrepancy.
One can make :math:`R` independent replications of RQMC points to
see how stable a computation is. From :math:`R` independent values
a t test (or bootstrap t test [8]_) then gives approximate confidence
intervals on the mean value.  Some RQMC methods produce a
root mean squared error that is actually :math:`o(1/n)` and smaller than
the rate seen in unrandomized QMC.  An intuitive explanation is
that the error is a sum of many small ones and random errors
cancel in a way that deterministic ones do not. RQMC also
has advantages on integrands that are singular or, for other
reasons, fail to be Riemann integrable.

(R)QMC cannot beat Bahkvalov's curse of dimension (see [9]_). For
any random or deterministic method, there are worst case functions
that will give it poor performance in high dimensions. A worst
case function for QMC might be 0 at all n points but very
large elsewhere. Worst case analyses get very pessimistic
in high dimensions. (R)QMC can bring a great improvement over
MC when the functions on which it is used are not worst case.
For instance (R)QMC can be especially effective on integrands
that are well approximated by sums of functions of one or two
or some small number of their input variables at a time [10]_, [11]_.
That property is often a surprising finding about those functions.

Scrambled nets are a kind of RQMC that have some valuable robustness
properties [12]_. If the integrand is square integrable they give variance
:math:`var_{SNET} = o(1/n)`. There is a finite upper bound on
:math:`var_{SNET} / var_{MC}` that holds simultaneously for every square
integrable integrand. Scrambled nets satisfy a strong law of large numbers
for :math:`f` in :math:`L^p` when :math:`p>1`. In some
special cases there is a central limit theorem [13]_. For smooth enough
integrands they can achieve RMSE nearly :math:`O(n^{-3})`. See [12]_
for references about these properties.

The main kinds of QMC methods are lattice rules [14]_ and digital
nets and sequences [2]_, [15]_. The theories meet up in polynomial
lattice rules [16]_ which can produce digital nets. Lattice rules
require some form of search for good constructions. For digital
nets there are widely used default constructions.

The most widely used QMC methods are Sobol' sequences [17]_.
These are digital nets. They are extensible in both n and d.
They can be scrambled.  The special sample sizes are powers
of 2. Another popular method are Halton sequences [18]_.
The constructions resemble those of digital nets. The earlier
dimensions have much better equidistribution properties than
later ones. There are essentially no special sample sizes.
They are not thought to be as accurate as Sobol' sequences.
They can be scrambled. The nets of Faure [19]_ are also widely
used. All dimensions are equally good but the special sample
sizes grow rapidly with dimension d. They can be scrambled.
The nets of Niederreiter and Xing [20]_ have the best asymptotic
properties but have not shown good empirical performance [21]_.

Higher order digital nets are formed by a digit interleaving process
in the digits of the constructed points. They can achieve higher
levels of asymptotic accuracy given higher smoothness conditions on :math:`f`
and they can be scrambled [22]_. There is little or no empirical work
showing the improved rate to be attained.

Using QMC is like using the entire period of a small random
number generator. The constructions are similar and so
therefore are the computational costs [23]_.

(R)QMC is sometimes improved by passing the points through
a baker's transformation (tent function) prior to using them.
That function has the form :math:`1-2|x-1/2|`.  As x goes from 0 to
1 this function goes from 0 to 1 and then back.  It is very
useful to produce a periodic function for lattice rules [14]_
and sometimes it improves the convergence rate [24]_.

It is not straightforward to apply QMC methods to Markov
chain Monte Carlo (MCMC).  We can think of MCMC as using
n=1 one point in :math:`[0,1]^{dim}` for very large d, with ergodic results
corresponding to :math:`dim \to \infty`.  One proposal is in [25]_
and under strong conditions an improved rate of convergence
has been shown [26]_.

Returning to Sobol' points: there are many versions depending
on what are called direction numbers. Those are the result of
searches and are tabulated. A very widely used set of direction
numbers come from [27]_. It is extensible in dimension up to
:math:`dim=21201`.

References
----------
.. [1] Owen, Art B. "Monte Carlo Book: the Quasi-Monte Carlo parts." (2019).
.. [2] Niederreiter, Harald. Random number generation and quasi-Monte Carlo
   methods. Society for Industrial and Applied Mathematics, 1992.
.. [3] Dick, Josef, Frances Y. Kuo, and Ian H. Sloan. "High-dimensional
   integration: the quasi-Monte Carlo way." Acta Numerica 22 (2013): 133.
.. [4] Aho, A. V., C. Aistleitner, T. Anderson, K. Appel, V. Arnol'd, N.
   Aronszajn, D. Asotsky et al. "W. Chen et al.(eds.), A Panorama of
   Discrepancy Theory (2014): 679. Sringer International Publishing,
   Switzerland.
.. [5] Hickernell, Fred J. "Koksma-Hlawka Inequality." Wiley StatsRef:
   Statistics Reference Online (2014).
.. [6] Owen, Art B. "On dropping the first Sobol' point." arXiv preprint
   arXiv:2008.08051 (2020).
.. [7] L'Ecuyer, Pierre, and Christiane Lemieux. "Recent advances in randomized
   quasi-Monte Carlo methods." In Modeling uncertainty, pp. 419-474. Springer,
   New York, NY, 2002.
.. [8] DiCiccio, Thomas J., and Bradley Efron. "Bootstrap confidence
   intervals." Statistical science (1996): 189-212.
.. [9] Dimov, Ivan T. Monte Carlo methods for applied scientists. World
   Scientific, 2008.
.. [10] Caflisch, Russel E., William J. Morokoff, and Art B. Owen. Valuation of
   mortgage backed securities using Brownian bridges to reduce effective
   dimension. Journal of Computational Finance, (1997): 1, no. 1 27-46.
.. [11] Sloan, Ian H., and Henryk Wozniakowski. "When are quasi-Monte Carlo
   algorithms efficient for high dimensional integrals?." Journal of Complexity
   14, no. 1 (1998): 1-33.
.. [12] Owen, Art B., and Daniel Rudolf "A strong law of large numbers for
   scrambled net integration." SIAM Review, to appear.
.. [13] Loh, Wei-Liem. "On the asymptotic distribution of scrambled net
   quadrature." The Annals of Statistics 31, no. 4 (2003): 1282-1324.
.. [14] Sloan, Ian H. and S. Joe. Lattice methods for multiple integration.
   Oxford University Press, 1994.
.. [15] Dick, Josef, and Friedrich Pillichshammer. Digital nets and sequences:
   discrepancy theory and quasi-Monte Carlo integration. Cambridge University
   Press, 2010.
.. [16] Dick, Josef, F. Kuo, Friedrich Pillichshammer, and I. Sloan.
   "Construction algorithms for polynomial lattice rules for multivariate
   integration." Mathematics of computation 74, no. 252 (2005): 1895-1921.
.. [17] Sobol', Il'ya Meerovich. "On the distribution of points in a cube and
   the approximate evaluation of integrals." Zhurnal Vychislitel'noi Matematiki
   i Matematicheskoi Fiziki 7, no. 4 (1967): 784-802.
.. [18] Halton, John H. "On the efficiency of certain quasi-random sequences of
   points in evaluating multi-dimensional integrals." Numerische Mathematik 2,
   no. 1 (1960): 84-90.
.. [19] Faure, Henri. "Discrepance de suites associees a un systeme de
   numeration (en dimension s)." Acta arithmetica 41, no. 4 (1982): 337-351.
.. [20] Niederreiter, Harold, and Chaoping Xing. "Low-discrepancy sequences and
   global function fields with many rational places." Finite Fields and their
   applications 2, no. 3 (1996): 241-273.
.. [21] Hong, Hee Sun, and Fred J. Hickernell. "Algorithm 823: Implementing
   scrambled digital sequences." ACM Transactions on Mathematical Software
   (TOMS) 29, no. 2 (2003): 95-109.
.. [22] Dick, Josef. "Higher order scrambled digital nets achieve the optimal
   rate of the root mean square error for smooth integrands." The Annals of
   Statistics 39, no. 3 (2011): 1372-1398.
.. [23] Niederreiter, Harald. "Multidimensional numerical integration using
   pseudorandom numbers." In Stochastic Programming 84 Part I, pp. 17-38.
   Springer, Berlin, Heidelberg, 1986.
.. [24] Hickernell, Fred J. "Obtaining O (N-2+e) Convergence for Lattice
   Quadrature Rules." In Monte Carlo and Quasi-Monte Carlo Methods 2000,
   pp. 274-289. Springer, Berlin, Heidelberg, 2002.
.. [25] Owen, Art B., and Seth D. Tribble. "A quasi-Monte Carlo Metropolis
   algorithm." Proceedings of the National Academy of Sciences 102, no. 25
   (2005): 8844-8849.
.. [26] Chen, Su. "Consistency and convergence rate of Markov chain quasi Monte
   Carlo with examples." PhD diss., Stanford University, 2011.
.. [27] Joe, Stephen, and Frances Y. Kuo. "Constructing Sobol sequences with
   better two-dimensional projections." SIAM Journal on Scientific Computing
   30, no. 5 (2008): 2635-2654.

"""
from abc import ABC, abstractmethod
import math
import warnings

import numpy as np

from scipy.optimize import brute
from scipy._lib._util import check_random_state
from scipy.optimize import basinhopping
from scipy.stats import norm
from scipy.stats._sobol import (
    initialize_v, _cscramble, _fill_p_cumulative, _draw, _fast_forward,
    _categorize, initialize_direction_numbers, _MAXDIM, _MAXBIT
)

__all__ = ['scale', 'discrepancy', 'QMCEngine', 'Sobol', 'Halton',
           'OrthogonalLatinHypercube', 'LatinHypercube', 'OptimalDesign',
           'multinomial_qmc', 'NormalQMC', 'MultivariateNormalQMC']


def scale(sample, bounds, reverse=False):
    """Sample scaling from unit hypercube to bounds range.

    To convert a sample from :math:`[0, 1)` to :math:`[a, b), b>a`, the
    following transformation is used:

    .. math::

        (b - a) * sample + a

    Parameters
    ----------
    sample : array_like (n_samples, dim)
        Sample to scale.
    bounds : tuple or array_like ([min, dim], [max, dim])
        Desired range of transformed data. The transformation apply the bounds
        on the sample and not the theoretical space, unit cube. Thus min and
        max values of the sample will coincide with the bounds.
    reverse : bool
        Reverse the transformation, from bounds range to unit hypercube.

    Returns
    -------
    sample : array_like (n_samples, dim)
        Scaled-sample.

    Examples
    --------
    >>> from scipy.stats import qmc
    >>> bounds = [[-2, 0],
    ...           [6, 5]]
    >>> sample = [[0.5 , 0.5 ],
    ...           [0.75, 0.25]]
    >>> qmc.scale(sample, bounds)
    array([[2.  , 2.5 ],
           [4.  , 1.25]])

    """
    bounds = np.asarray(bounds)
    min_ = np.min(bounds, axis=0)
    max_ = np.max(bounds, axis=0)
    if not reverse:
        return sample * (max_ - min_) + min_
    else:
        return (sample - min_) / (max_ - min_)


def discrepancy(sample, iterative=False, method='CD'):
    """Discrepancy on a given sample.

    Parameters
    ----------
    sample : array_like (n_samples, dim)
        The sample to compute the discrepancy from.
    iterative : bool
        Must be False if not using it for updating the discrepancy.
    method : str
        Type of discrepancy, can be ['CD', 'WD', 'MD', 'star'].

    Returns
    -------
    discrepancy : float
        Centered discrepancy.

    Notes
    -----
    The discrepancy is a uniformity criterion used to assess the space filling
    of a number of samples in a hypercube.
    The discrepancy measures how the spread of the points deviates from a
    uniform distribution.
    The lower the value is, the better the coverage of the parameter space is.

    Taking a subspace of the parameter space we count the number of points in
    the subpace and compare it to the total number of the points of the
    sample. This value is substracted by the volume of the subspace. This
    process is done over all subspaces and the highest value is kept.

    Four methods are available:

    * ``CD``: Centered Discrepancy - subspace involves a corner of the
      hypercube
    * ``WD``: Wrap-around Discrepancy - subspace can wrap around bounds
    * ``MD``: Mixture Discrepancy - mix between CD/WD covering more criteria
    * ``star``: Star L2-discrepancy - like CD BUT variant to rotation

    References
    ----------
    .. [1] Fang et al. Design and modeling for computer experiments,
       Computer Science and Data Analysis Series, 2006.
    .. [2] Zhou Y.-D. et al. Mixture discrepancy for quasi-random point sets
       Journal of Complexity, 29 (3-4) , pp. 283-301, 2013.
    .. [3] T. T. Warnock. Computational investigations of low discrepancy point
       sets, Applications of Number Theory to Numerical
       Analysis, Academic Press, pp. 319-343, 1972.

    Examples
    --------
    Calculate the quality of the sample using the discrepancy:

    >>> from scipy.stats import qmc
    >>> space = np.array([[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]])
    >>> bounds = np.array([[0.5, 0.5], [6.5, 6.5]])
    >>> space = qmc.scale(space, bounds, reverse=True)
    >>> space
    array([[0.08333333, 0.41666667],
           [0.25      , 0.91666667],
           [0.41666667, 0.25      ],
           [0.58333333, 0.75      ],
           [0.75      , 0.08333333],
           [0.91666667, 0.58333333]])
    >>> qmc.discrepancy(space)
    0.008142039609053464

    We can also compute iteratively the discrepancy by using
    ``iterative=True``.

    >>> disc_init = qmc.discrepancy(space[:-1], iterative=True)
    >>> disc_init
    0.04769081147119336

    """
    sample = np.asarray(sample)

    n_samples, dim = sample.shape

    if iterative:
        n_samples += 1

    if method == 'CD':
        abs_ = abs(sample - 0.5)
        disc1 = np.sum(np.prod(1 + 0.5 * abs_ - 0.5 * abs_ ** 2, axis=1))

        prod_arr = 1
        for i in range(dim):
            s0 = sample[:, i]
            prod_arr *= (1 +
                         0.5 * abs(s0[:, None] - 0.5) + 0.5 * abs(s0 - 0.5) -
                         0.5 * abs(s0[:, None] - s0))
        disc2 = prod_arr.sum()

        return ((13.0 / 12.0) ** dim - 2.0 / n_samples * disc1 +
                1.0 / (n_samples ** 2) * disc2)
    elif method == 'WD':
        prod_arr = 1
        for i in range(dim):
            s0 = sample[:, i]
            x_kikj = abs(s0[:, None] - s0)
            prod_arr *= 3.0 / 2.0 - x_kikj + x_kikj ** 2

        return - (4.0 / 3.0) ** dim + 1.0 / (n_samples ** 2) * prod_arr.sum()
    elif method == 'MD':
        abs_ = abs(sample - 0.5)
        disc1 = np.sum(np.prod(5.0 / 3.0 - 0.25 * abs_ - 0.25 * abs_ ** 2,
                               axis=1))

        prod_arr = 1
        for i in range(dim):
            s0 = sample[:, i]
            prod_arr *= (15.0 / 8.0 -
                         0.25 * abs(s0[:, None] - 0.5) - 0.25 * abs(s0 - 0.5) -
                         3.0 / 4.0 * abs(s0[:, None] - s0) +
                         0.5 * abs(s0[:, None] - s0) ** 2)
        disc2 = prod_arr.sum()

        disc = (19.0 / 12.0) ** dim
        disc1 = 2.0 / n_samples * disc1
        disc2 = 1.0 / (n_samples ** 2) * disc2

        return disc - disc1 + disc2
    elif method == 'star':
        return np.sqrt(
            3 ** (-dim) - 2 ** (1 - dim) / n_samples
            * np.sum(np.prod(1 - sample ** 2, axis=1))
            + np.sum([
                np.prod(1 - np.maximum(sample[k, :], sample[j, :]))
                for k in range(n_samples) for j in range(n_samples)
            ]) / n_samples ** 2
        )
    else:
        raise ValueError('{} is not a valid method. Options are '
                         'CD, WD, MD, star.'.format(method))


def _update_discrepancy(x_new, sample, initial_disc):
    """Update the centered discrepancy with a new sample.

    Parameters
    ----------
    x_new : array_like (1, dim)
        The new sample to add in `sample`.
    sample : array_like (n_samples, dim)
        The initial sample.
    initial_disc : float
        Centered discrepancy of the `sample`.

    Returns
    -------
    discrepancy : float
        Centered discrepancy of the sample composed of `x_new` and `sample`.

    """
    sample = np.asarray(sample)
    x_new = np.asarray(x_new)

    n_samples = len(sample) + 1
    abs_ = abs(x_new - 0.5)

    disc1 = - 2 / n_samples * np.prod(1 + 1 / 2 * abs_ - 1 / 2 * abs_ ** 2)
    disc2 = 2 / (n_samples ** 2) * np.sum(np.prod(1 + 1 / 2 * abs_ +
                                                  1 / 2 * abs(sample - 0.5) -
                                                  1 / 2 * abs(x_new - sample),
                                                  axis=1))
    disc3 = 1 / (n_samples ** 2) * np.prod(1 + abs_)

    return initial_disc + disc1 + disc2 + disc3


def _perturb_discrepancy(sample, i1, i2, k, disc):
    """Centered discrepancy after and elementary perturbation on a LHS.

    An elementary perturbation consists of an exchange of coordinates between
    two points: ``sample[i1, k] <-> sample[i2, k]``. By construction,
    this operation conserves the LHS properties.

    Parameters
    ----------
    sample : array_like (n_samples, dim)
        The sample (before permutation) to compute the discrepancy from.
    i1 : int
        The first line of the elementary permutation.
    i2 : int
        The second line of the elementary permutation.
    k : int
        The column of the elementary permutation.
    disc : float
        Centered discrepancy of the design before permutation.

    Returns
    -------
    discrepancy : float
        Centered discrepancy.

    References
    ----------
    .. [1] Jin et al. "An efficient algorithm for constructing optimal design
       of computer experiments", Journal of Statistical Planning and
       Inference, 2005.

    """
    sample = np.asarray(sample)
    n_samples = sample.shape[0]

    z_ij = sample - 0.5

    # Eq (19)
    c_i1j = 1. / n_samples ** 2. * np.prod(0.5 * (2. + abs(z_ij[i1, :]) +
                                                  abs(z_ij) -
                                                  abs(z_ij[i1, :] - z_ij)),
                                           axis=1)
    c_i2j = 1. / n_samples ** 2. * np.prod(0.5 * (2. + abs(z_ij[i2, :]) +
                                                  abs(z_ij) -
                                                  abs(z_ij[i2, :] - z_ij)),
                                           axis=1)

    # Eq (20)
    c_i1i1 = (1. / n_samples ** 2 * np.prod(1 + abs(z_ij[i1, :])) -
              2. / n_samples * np.prod(1. + 0.5 * abs(z_ij[i1, :]) -
                                       0.5 * z_ij[i1, :] ** 2))
    c_i2i2 = (1. / n_samples ** 2 * np.prod(1 + abs(z_ij[i2, :])) -
              2. / n_samples * np.prod(1. + 0.5 * abs(z_ij[i2, :]) -
                                       0.5 * z_ij[i2, :] ** 2))

    # Eq (22), typo in the article in the denominator i2 -> i1
    num = (2 + abs(z_ij[i2, k]) + abs(z_ij[:, k]) -
           abs(z_ij[i2, k] - z_ij[:, k]))
    denum = (2 + abs(z_ij[i1, k]) + abs(z_ij[:, k]) -
             abs(z_ij[i1, k] - z_ij[:, k]))
    gamma = num / denum

    # Eq (23)
    c_p_i1j = gamma * c_i1j
    # Eq (24)
    c_p_i2j = c_i2j / gamma

    alpha = (1 + abs(z_ij[i2, k])) / (1 + abs(z_ij[i1, k]))
    beta = (2 - abs(z_ij[i2, k])) / (2 - abs(z_ij[i1, k]))

    g_i1 = np.prod(1. + abs(z_ij[i1, :]))
    g_i2 = np.prod(1. + abs(z_ij[i2, :]))
    h_i1 = np.prod(1. + 0.5 * abs(z_ij[i1, :]) - 0.5 * (z_ij[i1, :] ** 2))
    h_i2 = np.prod(1. + 0.5 * abs(z_ij[i2, :]) - 0.5 * (z_ij[i2, :] ** 2))

    # Eq (25), typo in the article g is missing
    c_p_i1i1 = ((g_i1 * alpha) / (n_samples ** 2) -
                2. * alpha * beta * h_i1 / n_samples)
    # Eq (26), typo in the article n ** 2
    c_p_i2i2 = ((g_i2 / ((n_samples ** 2) * alpha)) -
                (2. * h_i2 / (n_samples * alpha * beta)))

    # Eq (26)
    sum_ = c_p_i1j - c_i1j + c_p_i2j - c_i2j

    mask = np.ones(n_samples, dtype=bool)
    mask[[i1, i2]] = False
    sum_ = sum(sum_[mask])

    disc_ep = (disc + c_p_i1i1 - c_i1i1 + c_p_i2i2 - c_i2i2 + 2 * sum_)

    return disc_ep


def primes_from_2_to(n):
    """Prime numbers from 2 to *n*.

    Parameters
    ----------
    n : int
        Sup bound with ``n >= 6``.

    Returns
    -------
    primes : list(int)
        Primes in ``2 <= p < n``.

    References
    ----------
    .. [1] `StackOverflow <https://stackoverflow.com/questions/2068372>`_.

    """
    sieve = np.ones(n // 3 + (n % 6 == 2), dtype=bool)
    for i in range(1, int(n ** 0.5) // 3 + 1):
        k = 3 * i + 1 | 1
        sieve[k * k // 3::2 * k] = False
        sieve[k * (k - 2 * (i & 1) + 4) // 3::2 * k] = False
    return np.r_[2, 3, ((3 * np.nonzero(sieve)[0][1:] + 1) | 1)]


def n_primes(n):
    """List of the n-first prime numbers.

    Parameters
    ----------
    n : int
        Number of prime numbers wanted.

    Returns
    -------
    primes : list(int)
        List of primes.

    """
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
              61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
              131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
              197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
              271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
              353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431,
              433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
              509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599,
              601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673,
              677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761,
              769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857,
              859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
              953, 967, 971, 977, 983, 991, 997][:n]

    if len(primes) < n:
        big_number = 2000
        while 'Not enough primes':
            primes = primes_from_2_to(big_number)[:n]
            if len(primes) == n:
                break
            big_number += 1000

    return primes


def van_der_corput(n_samples, base=2, start_index=0, scramble=False,
                   seed=None):
    """Van der Corput sequence.

    Pseudo-random number generator based on a b-adic expansion.

    Scrambling uses permutations of the remainders (see [1]_). Multiple
    permutations are applied to construct a point. The sequence of
    permutations has to be the same for all points of the sequence.

    Parameters
    ----------
    n_samples : int
        Number of element of the sequence.
    base : int
        Base of the sequence.
    start_index : int
        Index to start the sequence from.
    scramble: bool, optional
        If True, use Owen scrambling.
    seed : {int or `numpy.random.RandomState` instance}, optional
        If `seed` is not specified the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``RandomState`` instance, then that
        instance is used.

    Returns
    -------
    sequence : list (n_samples,)
        Sequence of Van der Corput.

    References
    ----------
    .. [1] A. B. Owen. "A randomized Halton algorithm in R",
       arXiv:1706.02808, 2017.

    """
    rng = check_random_state(seed)
    sequence = np.zeros(n_samples)

    quotient = np.arange(start_index, start_index + n_samples)
    b2r = 1 / base

    while (1 - b2r) < 1:
        remainder = quotient % base

        if scramble:
            # permutation must be the same for all points of the sequence
            perm = rng.permutation(base)
            remainder = perm[np.array(remainder).astype(int)]

        sequence += remainder * b2r
        b2r /= base
        quotient = (quotient - remainder) / base

    return sequence


class QMCEngine(ABC):
    """A generic Quasi-Monte Carlo sampler class meant for subclassing.

    QMCEngine is a base class to construct specific Quasi-Monte Carlo sampler.
    It cannot be used directly as a sampler.

    Parameters
    ----------
    dim : int
        Dimension of the parameter space.
    seed : {int or `numpy.random.RandomState` instance}, optional
        If `seed` is not specified the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``RandomState`` instance, then that
        instance is used.

    Notes
    -----
    By convention samples are distributed over the half-open interval
    ``[0, 1)``. Instances of the class can access the attributes: ``dim`` for
    the dimension; and ``rng`` for the random number generator (used for the
    ``seed``).

    **Subclassing**

    When subclassing `QMCEngine` to create a new sampler,  ``__init__`` and
    ``random`` has to be redefined.

    * ``__init__(dim, seed=None)``: at least fix the dimension. If the sampler
      does not take advantage of a ``seed`` (deterministic methods like
      Halton), this parameter can be omitted.
    * ``random(n_samples)``: draw ``n_samples`` from the engine.

    Optionally, 2 other methods can be overwritten by subclasses:

    * ``reset``: Reset the engine to it's original state.
    * ``fast_forward``: It should be used as a way to fast-forward a sequence
      to a further state. If the sequence is deterministic (like Halton
      sequence), then ``fast_forward(n)`` is skipping the ``n`` first draw.

    Examples
    --------
    To create a random sampler based on ``np.random.random``, we would do the
    following:

    >>> from scipy.stats import qmc
    >>> class RandomEngine(qmc.QMCEngine):
    ...     def __init__(self, dim, seed):
    ...         super().__init__(dim=dim, seed=seed)
    ...
    ...
    ...     def random(self, n_samples=1):
    ...         return self.rng.random((n_samples, self.dim))
    ...

    We subclass `QMCEngine` by defining the sampling strategy we want to use.
    And we can create an instance to sample from.

    >>> engine = RandomEngine(2, seed=12345)
    >>> engine.random(5)
    array([[0.92961609, 0.31637555],
           [0.18391881, 0.20456028],
           [0.56772503, 0.5955447 ],
           [0.96451452, 0.6531771 ],
           [0.74890664, 0.65356987]])

    """

    @abstractmethod
    def __init__(self, dim, seed=None):
        self.dim = dim
        self.rng = check_random_state(seed)
        self.num_generated = 0

    @abstractmethod
    def random(self, n_samples=1):
        """Draw `n_samples` in the half-open interval ``[0, 1)``.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, dim)
            QMC sample.
        """
        # self.num_generated += n_samples

    def reset(self):
        """Reset the engine to base state.

        Returns
        -------
        engine: QMCEngine
            Engine reset to its base state.

        """
        self.num_generated = 0
        return self

    def fast_forward(self, n):
        """Fast-forward the sequence by `n` positions.

        Parameters
        ----------
        n: int
            Number of points to skip in the sequence.
        """
        self.num_generated += n
        return self


class Halton(QMCEngine):
    """Halton sequence.

    Pseudo-random number generator that generalize the Van der Corput sequence
    for multiple dimensions. Halton sequence use base-two Van der Corput
    sequence for the first dimension, base-three for its second and base-n for
    its n-dimension.

    Parameters
    ----------
    dim : int
        Dimension of the parameter space.
    scramble: bool, optional
        If True, use Owen scrambling.
    seed : {int or `numpy.random.RandomState` instance}, optional
        If `seed` is not specified the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``RandomState`` instance, then that
        instance is used.

    Notes
    -----
    The Halton sequence has severe striping artifacts for even modestly
    large dimensions. These can be ameliorated by scrambling. Scrambling
    also supports replication-based error estimates and extends
    applicabiltiy to unbounded integrands.

    References
    ----------
    .. [1] Halton, "On the efficiency of certain quasi-random sequences of
       points in evaluating multi-dimensional integrals", Numerische
       Mathematik, 1960.
    .. [2] A. B. Owen. "A randomized Halton algorithm in R",
       arXiv:1706.02808, 2017.

    Examples
    --------
    Generate samples from a low discrepancy sequence of Halton.

    >>> from scipy.stats import qmc
    >>> sampler = qmc.Halton(dim=2, scramble=False)
    >>> sample = sampler.random(n_samples=5)
    >>> sample
    array([[0.        , 0.        ],
           [0.5       , 0.33333333],
           [0.25      , 0.66666667],
           [0.75      , 0.11111111],
           [0.125     , 0.44444444]])

    Compute the quality of the sample using the discrepancy criterion.

    >>> qmc.discrepancy(sample)
    0.088893711419753

    If some wants to continue an existing design, extra points can be obtained
    by calling again `random`. Alternatively, you can skip some points like:

    >>> sampler.fast_forward(5)
    >>> sample_continued = sampler.random(n_samples=5)
    >>> sample_continued
    array([[0.3125    , 0.37037037],
           [0.8125    , 0.7037037 ],
           [0.1875    , 0.14814815],
           [0.6875    , 0.48148148],
           [0.4375    , 0.81481481]])

    Finally, samples can be scaled to bounds.

    >>> bounds = [[0, 2], [10, 5]]
    >>> qmc.scale(sample_continued, bounds)
    array([[3.125     , 3.11111111],
           [8.125     , 4.11111111],
           [1.875     , 2.44444444],
           [6.875     , 3.44444444],
           [4.375     , 4.44444444]])

    """

    def __init__(self, dim, scramble=True, seed=None):
        super().__init__(dim=dim)
        self.seed = seed
        self.base = n_primes(dim)
        self.scramble = scramble

    def random(self, n_samples=1):
        """Draw `n_samples` in the half-open interval ``[0, 1)``.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, dim)
            QMC sample.

        """
        # Generate a sample using a Van der Corput sequence per dimension.
        # important to have ``type(bdim) == int`` for performance reason
        sample = [van_der_corput(n_samples, int(bdim), self.num_generated,
                                 scramble=self.scramble, seed=self.seed)
                  for bdim in self.base]

        self.num_generated += n_samples
        return np.array(sample).T


class OrthogonalLatinHypercube(QMCEngine):
    """Orthogonal array-based Latin hypercube sampling (OA-LHS).

    On top of the constraints from the Latin Hypercube, an orthogonal array of
    size n_samples is defined and only one point is allowed per subspace.

    Parameters
    ----------
    dim : int
        Dimension of the parameter space.
    seed : {int or `numpy.random.RandomState` instance}, optional
        If `seed` is not specified the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``RandomState`` instance, then that
        instance is used.

    References
    ----------
    .. [1] Art B. Owen, "Orthogonal arrays for computer experiments,
       integration and visualization", Statistica Sinica, 1992.

    Examples
    --------
    Generate samples from an orthogonal latin hypercube generator.

    >>> from scipy.stats import qmc
    >>> sampler = qmc.OrthogonalLatinHypercube(dim=2, seed=12345)
    >>> sample = sampler.random(n_samples=5)
    >>> sample
    array([[0.18592322, 0.77846875],
           [0.64091206, 0.18474763],
           [0.91354501, 0.80535794],
           [0.43678376, 0.27478424],
           [0.26327511, 0.43099467]])

    Compute the quality of the sample using the discrepancy criterion.

    >>> qmc.discrepancy(sample)
    0.02004864477993462

    Finally, samples can be scaled to bounds.

    >>> bounds = [[0, 2], [10, 5]]
    >>> qmc.scale(sample, bounds)
    array([[1.85923219, 4.33540624],
           [6.40912056, 2.55424289],
           [9.13545006, 4.41607383],
           [4.36783762, 2.82435271],
           [2.63275111, 3.29298401]])

    """

    def __init__(self, dim, seed=None):
        super().__init__(dim=dim, seed=seed)

    def random(self, n_samples=1):
        """Draw `n_samples` in the half-open interval ``[0, 1)``.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, dim)
            OLHS sample.

        """
        sample = []
        step = 1.0 / n_samples

        for _ in range(self.dim):
            # Enforce a unique point per grid
            j = np.arange(n_samples) * step
            temp = j + self.rng.uniform(low=0, high=step, size=n_samples)
            self.rng.shuffle(temp)

            sample.append(temp)

        return np.array(sample).T


class LatinHypercube(QMCEngine):
    """Latin hypercube sampling (LHS).

    A Latin hypercube sample [1]_ generates ``n`` points in
    :math:`[0,1)^{dim}`. Each univariate marginal distribution is stratified,
    placing exactly one point in :math:`[j/n, (j+1)/n)` for
    :math:`j=0,1,...,n-1`. They are still applicable when :math:`n << dim`.
    LHS is extremely effective on integrands that are nearly additive [2]_.
    LHS on n points never has more variance than plain MC on :math:`n-1`
    points [3]_. There is a central limit theorem for plain LHS [4]_, but not
    necessarily for optimized LHS.

    Parameters
    ----------
    dim : int
        Dimension of the parameter space.
    centered : bool
        Center the point within the multi-dimensional grid.
    seed : {int or `numpy.random.RandomState` instance}, optional
        If `seed` is not specified the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``RandomState`` instance, then that
        instance is used.

    References
    ----------
    .. [1] Mckay et al., "A Comparison of Three Methods for Selecting Values
       of Input Variables in the Analysis of Output from a Computer Code",
       Technometrics, 1979.
    .. [2] M. Stein, "Large sample properties of simulations using Latin
       hypercube sampling." Technometrics 29, no. 2: 143-151, 1987.
    .. [3] A. B. Owen, "Monte Carlo variance of scrambled net quadrature."
       SIAM Journal on Numerical Analysis 34, no. 5: 1884-1910, 1997
    .. [4]  Loh, W.-L. "On Latin hypercube sampling." The annals of statistics
       24, no. 5: 2058-2080, 1996.

    Examples
    --------
    Generate samples from a latin hypercube generator.

    >>> from scipy.stats import qmc
    >>> sampler = qmc.LatinHypercube(dim=2, seed=12345)
    >>> sample = sampler.random(n_samples=5)
    >>> sample
    array([[0.01407678, 0.53672489],
           [0.36321624, 0.75908794],
           [0.28645499, 0.48089106],
           [0.2070971 , 0.46936458],
           [0.45021867, 0.66928603]])

    Compute the quality of the sample using the discrepancy criterion.

    >>> qmc.discrepancy(sample)
    0.1271335273223828

    Finally, samples can be scaled to bounds.

    >>> bounds = [[0, 2], [10, 5]]
    >>> qmc.scale(sample, bounds)
    array([[0.14076781, 3.61017467],
           [3.63216238, 4.27726383],
           [2.86454994, 3.44267318],
           [2.07097096, 3.40809374],
           [4.50218672, 4.00785808]])

    """

    def __init__(self, dim, centered=False, seed=None):
        super().__init__(dim=dim, seed=seed)
        self.centered = centered

    def random(self, n_samples=1):
        """Draw `n_samples` in the half-open interval ``[0, 1)``.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, dim)
            LHS sample.

        """
        if self.centered:
            r = 0.5
        else:
            r = self.rng.random_sample((n_samples, self.dim))

        q = self.rng.randint(low=1, high=n_samples,
                             size=(n_samples, self.dim))

        return 1. / n_samples * (q - r)


class OptimalDesign(QMCEngine):
    """Optimal design.

    Optimize the design by doing random permutations to lower the centered
    discrepancy. If `optimization` is False, `niter` design are generated and
    the one with lowest centered discrepancy is return. This option is faster.

    Centered discrepancy based design show better space filling robustness
    toward 2D and 3D subprojections. Distance based design better space
    filling but less robust to subprojections.

    Parameters
    ----------
    dim : int
        Dimension of the parameter space.
    start_design : array_like (n_samples, dim)
        Initial design of experiment to optimize.
    niter : int
        Number of iteration to perform.
    force : bool
        If `optimization`, force *basinhopping* optimization. Otherwise
        grid search is used.
    optimization : bool
        Optimal design using global optimization or random generation of
        `niter` samples.
    seed : {int or `numpy.random.RandomState` instance}, optional
        If `seed` is not specified the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``RandomState`` instance, then that
        instance is used.

    References
    ----------
    .. [1] Fang et al. Design and modeling for computer experiments,
       Computer Science and Data Analysis Series, 2006.
    .. [2] Damblin et al., "Numerical studies of space filling designs:
       optimization of Latin Hypercube Samples and subprojection properties",
       Journal of Simulation, 2013.

    Examples
    --------
    Generate samples from an optimal design.

    >>> from scipy.stats import qmc
    >>> sampler = qmc.OptimalDesign(dim=2, seed=12345)
    >>> sample = sampler.random(n_samples=5)
    >>> sample
    array([[0.64091206, 0.77846875],
           [0.43678376, 0.18474763],
           [0.18592322, 0.80535794],
           [0.91354501, 0.27478424],
           [0.26327511, 0.43099467]])

    Compute the quality of the sample using the discrepancy criterion.

    >>> qmc.discrepancy(sample)
    0.019688311022535432

    You can possibly improve the quality of the sample by performing more
    optimization iterations by using `niter`:

    >>> sampler_2 = qmc.OptimalDesign(dim=2, niter=2, seed=12345)
    >>> sample_2 = sampler_2.random(n_samples=5)
    >>> qmc.discrepancy(sample_2)
    0.019607673478802434

    Finally, samples can be scaled to bounds.

    >>> bounds = [[0, 2], [10, 5]]
    >>> qmc.scale(sample, bounds)
    array([[6.40912056, 4.33540624],
           [4.36783762, 2.55424289],
           [1.85923219, 4.41607383],
           [9.13545006, 2.82435271],
           [2.63275111, 3.29298401]])

    """

    def __init__(self, dim, start_design=None, niter=1, force=False,
                 optimization=True, seed=None):
        super().__init__(dim=dim, seed=seed)
        self.start_design = start_design
        self.niter = niter
        self.force = force
        self.optimization = optimization

        self.best_doe = self.start_design
        if self.start_design is not None:
            self.best_disc = discrepancy(self.start_design)
        else:
            self.best_disc = np.inf

        self.olhs = OrthogonalLatinHypercube(self.dim, seed=self.rng)

    def random(self, n_samples=1):
        """Draw `n_samples` in the half-open interval ``[0, 1)``.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, dim)
            Optimal sample.

        """
        if self.optimization:
            if self.best_doe is None:
                self.best_doe = self.olhs.random(n_samples)
                self.best_disc = discrepancy(self.best_doe)

            def _perturb_best_doe(x):
                """Perturb the DoE and keep track of the best DoE.

                Parameters
                ----------
                x : list of int
                    It is a list of:
                        idx : int
                            Index value of the components to compute

                Returns
                -------
                discrepancy : float
                    Centered discrepancy.

                """
                # Perturb the DoE
                doe = self.best_doe.copy()
                col, row_1, row_2 = np.round(x).astype(int)
                doe[row_1, col], doe[row_2, col] = doe[row_2, col],\
                    doe[row_1, col]

                disc = _perturb_discrepancy(self.best_doe, row_1, row_2, col,
                                            self.best_disc)

                if disc < self.best_disc:
                    self.best_disc = disc
                    self.best_doe = doe

                return disc

            # Total number of possible design
            complexity = self.dim * n_samples ** 2

            if (complexity > 1e6) or self.force:
                bounds_optim = ([0, self.dim - 1],
                                [0, n_samples - 1],
                                [0, n_samples - 1])
            else:
                bounds_optim = (slice(0, self.dim - 1, 1),
                                slice(0, n_samples - 1, 1),
                                slice(0, n_samples - 1, 1))

            for _ in range(self.niter):
                if (complexity > 1e6) or self.force:
                    minimizer_kwargs = {"method": "L-BFGS-B",
                                        "bounds": bounds_optim}
                    _ = basinhopping(_perturb_best_doe, [0, 0, 0], niter=100,
                                     minimizer_kwargs=minimizer_kwargs)
                else:
                    _ = brute(_perturb_best_doe, ranges=bounds_optim,
                              finish=None)
        else:
            for _ in range(self.niter):
                doe = self.olhs.random(n_samples)
                disc = discrepancy(doe)
                if disc < self.best_disc:
                    self.best_disc = disc
                    self.best_doe = doe

        return self.best_doe


class Sobol(QMCEngine):
    """Engine for generating (scrambled) Sobol' sequences.

    Sobol' sequences are low-discrepancy, quasi-random numbers. Points
    can be drawn using two methods:

    * `random_base2`: safely draw :math:`n=2^m` points. This method
      guaranty the balance properties of the sequence.
    * `random`: draw an arbitrary number of points from the
      sequence.

    Parameters
    ----------
    dim: int
        Dimensionality of the sequence. Max dimensionality is 21201.
    scramble: bool, optional
        If True, use Owen scrambling.
    seed : {int or `numpy.random.RandomState` instance}, optional
        If `seed` is not specified the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``RandomState`` instance, then that
        instance is used.

    Notes
    -----
    Sobol' sequences [1]_ provide :math:`n=2^m` low discrepancy points in
    :math:`[0,1)^{dim}`. Scrambling them [2]_ makes them suitable for singular
    integrands, provides a means of error estimation, and can improve their
    rate of convergence.

    There are many versions of Sobol' sequences depending on their
    'direction numbers'. This code uses direction numbers from [3]_. Hence,
    maximum number of dimension is 21201. The direction numbers have been
    precomputed with search criterion 6 and can be retrieved at
    https://web.maths.unsw.edu.au/~fkuo/sobol/.

    .. warning::

       Sobol' sequences are a quadrature rule and they lose their balance
       properties if one uses a sample size that is not a power of 2, or skips
       the first point, or thins the sequence [4]_.

       If :math:`n=2^m` points are not enough then one should take :math:`2^M`
       points for :math:`M>m`. When scrambling, the number R of independent
       replicates does not have to be a power of 2.

       Sobol' sequences are generated to some number :math:`B` of bits. Then
       after :math:`2^B` points have been generated, the sequence will repeat.
       Currently :math:`B=30`.

    References
    ----------
    .. [1] I. M. Sobol. The distribution of points in a cube and the accurate
       evaluation of integrals. Zh. Vychisl. Mat. i Mat. Phys., 7:784-802,
       1967.

    .. [2] Art B. Owen. Scrambling Sobol and Niederreiter-Xing points.
       Journal of Complexity, 14(4):466-489, December 1998.

    .. [3] S. Joe and F. Y. Kuo. Constructing sobol sequences with better
       two-dimensional projections. SIAM Journal on Scientific Computing,
       30(5):2635-2654, 2008.

    .. [4] Art B. Owen. On dropping the first Sobol' point. arXiv 2008.08051,
       2020.

    Examples
    --------
    Generate samples from a low discrepancy sequence of Sobol'.

    >>> from scipy.stats import qmc
    >>> sampler = qmc.Sobol(dim=2, scramble=False)
    >>> sample = sampler.random_base2(m=3)
    >>> sample
    array([[0.   , 0.   ],
           [0.5  , 0.5  ],
           [0.75 , 0.25 ],
           [0.25 , 0.75 ],
           [0.375, 0.375],
           [0.875, 0.875],
           [0.625, 0.125],
           [0.125, 0.625]])

    Compute the quality of the sample using the discrepancy criterion.

    >>> qmc.discrepancy(sample)
    0.013882107204860938

    If some wants to continue an existing design, extra points can be obtained
    by calling again `random_base2`. Alternatively, you can skip some
    points like:

    >>> sampler.reset()
    >>> sampler.fast_forward(4)
    >>> sample_continued = sampler.random_base2(m=2)
    >>> sample_continued
    array([[0.375, 0.375],
           [0.875, 0.875],
           [0.625, 0.125],
           [0.125, 0.625]])

    Finally, samples can be scaled to bounds.

    >>> bounds = [[0, 2], [10, 5]]
    >>> qmc.scale(sample_continued, bounds)
    array([[3.75 , 3.125],
           [8.75 , 4.625],
           [6.25 , 2.375],
           [1.25 , 3.875]])

    """

    MAXDIM = _MAXDIM
    MAXBIT = _MAXBIT

    def __init__(self, dim, scramble=True, seed=None):
        if dim > self.MAXDIM:
            raise ValueError(
                "Maximum supported dimensionality is {}.".format(self.MAXDIM)
            )
        super().__init__(dim=dim, seed=seed)

        # initialize direction numbers
        initialize_direction_numbers()

        # v is dim x MAXBIT matrix
        self._sv = np.zeros((dim, self.MAXBIT), dtype=int)
        initialize_v(self._sv, dim)

        if not scramble:
            self._shift = np.zeros(dim, dtype=int)
        else:
            self._scramble()

        self._quasi = self._shift.copy()
        self._first_point = (self._quasi / 2 ** self.MAXBIT).reshape(1, -1)

    def _scramble(self):
        """Scramble the sequence."""
        # Generate shift vector
        self._shift = np.dot(
            self.rng.randint(2, size=(self.dim, self.MAXBIT)),
            2 ** np.arange(self.MAXBIT),
        )
        self._quasi = self._shift.copy()
        # Generate lower triangular matrices (stacked across dimensions)
        ltm = np.tril(self.rng.randint(2, size=(self.dim,
                                                self.MAXBIT,
                                                self.MAXBIT)))
        _cscramble(self.dim, ltm, self._sv)
        self.num_generated = 0

    def random(self, n_samples=1):
        """Draw next point(s) in the Sobol' sequence.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, dim)
            Sobol' sample.

        """
        sample = np.empty((n_samples, self.dim), dtype=float)

        if self.num_generated == 0:
            # verify n_samples is 2**n
            if not (n_samples & (n_samples - 1) == 0):
                warnings.warn("The balance properties of Sobol' points require"
                              " n_samples to be a power of 2.")

            if n_samples == 1:
                sample = self._first_point
            else:
                _draw(n_samples - 1, self.num_generated, self.dim, self._sv,
                      self._quasi, sample)
                sample = np.concatenate([self._first_point,
                                         sample])[:n_samples]
        else:
            _draw(n_samples, self.num_generated - 1, self.dim, self._sv,
                  self._quasi, sample)

        self.num_generated += n_samples
        return sample

    def random_base2(self, m=1):
        """Draw point(s) from the Sobol' sequence.

        This function draws :math:`n=2^m` points in the parameter space
        ensuring the balance properties of the sequence.

        Parameters
        ----------
        m : int
            Logarithm in base 2 of the number of samples; i.e., n = 2^m.

        Returns
        -------
        sample : array_like (n_samples, dim)
            Sobol' sample.

        """
        n_samples = 2 ** m

        total_n_samples = self.num_generated + n_samples
        if not (total_n_samples & (total_n_samples - 1) == 0):
            raise ValueError("The balance properties of Sobol' points require "
                             "n to be a power of 2. {0} points have been "
                             "previously generated, then: n={0}+2**{1}={2}. "
                             "If you still want to do this, the function "
                             "'Sobol.random()' can be used."
                             .format(self.num_generated, m, total_n_samples))

        return self.random(n_samples=n_samples)

    def reset(self):
        """Reset the engine to base state.

        Returns
        -------
        engine: Sobol
            Engine reset to its base state.

        """
        self._quasi = self._shift.copy()
        self.num_generated = 0
        return self

    def fast_forward(self, n):
        """Fast-forward the sequence by `n` positions.

        Parameters
        ----------
        n: int
            Number of points to skip in the sequence.

        Returns
        -------
        engine: Sobol
            The fast-forwarded engine.

        """
        if self.num_generated == 0:
            _fast_forward(n - 1, self.num_generated, self.dim,
                          self._sv, self._quasi)
        else:
            _fast_forward(n, self.num_generated - 1, self.dim,
                          self._sv, self._quasi)
        self.num_generated += n
        return self


def multinomial_qmc(n_samples, pvals, engine=None, seed=None):
    """Draw low-discreancy quasi-random samples from multinomial distribution.

    Parameters
    ----------
    n_samples : int
        Number of experiments.
    pvals: Iterable[float]
        float vector of probabilities of size ``p``. Elements must be
        non-negative and sum to 1.
    engine: QMCEngine
        Quasi-Monte Carlo engine sampler. If None, Sobol' is used.
    seed : {int or `numpy.random.RandomState` instance}, optional
        If `seed` is not specified the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``RandomState`` instance, then that
        instance is used.

    Returns
    -------
    samples: array_like (pvals,)
        int vector of size ``p`` summing to ``n_samples``.

    """
    if np.min(pvals) < 0:
        raise ValueError('Elements of pvals must be non-negative.')
    if not np.isclose(np.sum(pvals), 1):
        raise ValueError('Elements of pvals must sum to 1.')

    if engine is None:
        engine = Sobol(1, scramble=True, seed=seed)
    draws = engine.random(n_samples).ravel()
    p_cumulative = np.empty_like(pvals, dtype=float)
    _fill_p_cumulative(np.array(pvals, dtype=float), p_cumulative)
    sample = np.zeros_like(pvals, dtype=int)
    _categorize(draws, p_cumulative, sample)
    return sample


class NormalQMC(QMCEngine):
    """Engine for QMC sampling from a multivariate normal :math:`N(0, I_d)`.

    Parameters
    ----------
    dim: int
        The dimension of the samples.
    inv_transform: bool
        If True, use inverse transform instead of Box-Muller.
    engine: QMCEngine
        Quasi-Monte Carlo engine sampler. If None, Sobol' is used.
    seed : {int or `numpy.random.RandomState` instance}, optional
        If `seed` is not specified the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``RandomState`` instance, then that
        instance is used.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from scipy.stats import qmc
    >>> engine = qmc.NormalQMC(2)
    >>> sample = engine.random(512)
    >>> plt.scatter(sample[:, 0], sample[:, 1])
    >>> plt.show()

    """

    def __init__(self, dim, inv_transform=True, engine=None, seed=None):
        super().__init__(dim=dim, seed=seed)
        self._inv_transform = inv_transform
        if not inv_transform:
            # to apply Box-Muller, we need an even number of dimensions
            engine_dim = 2 * math.ceil(dim / 2)
        else:
            engine_dim = dim

        if engine is None:
            self.engine = Sobol(dim=engine_dim, scramble=True, seed=seed)
        else:
            self.engine = engine

    def random(self, n_samples=1):
        """Draw `n_samples` QMC samples from the standard Normal.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, dim)
            Sample.

        """
        # get base samples
        samples = self.engine.random(n_samples)
        if self._inv_transform:
            # apply inverse transform
            # (values to close to 0/1 result in inf values)
            return norm.ppf(0.5 + (1 - 1e-10) * (samples - 0.5))
        else:
            # apply Box-Muller transform (note: indexes starting from 1)
            even = np.arange(0, samples.shape[-1], 2)
            Rs = np.sqrt(-2 * np.log(samples[:, even]))
            thetas = 2 * math.pi * samples[:, 1 + even]
            cos = np.cos(thetas)
            sin = np.sin(thetas)
            transf_samples = np.stack([Rs * cos, Rs * sin],
                                      -1).reshape(n_samples, -1)
            # make sure we only return the number of dimension requested
            return transf_samples[:, : self.dim]


class MultivariateNormalQMC(QMCEngine):
    r"""QMC sampling from a multivariate Normal :math:`N(\mu, \Sigma)`.

    Parameters
    ----------
    mean: array_like (dim,)
        The mean vector.
    cov: array_like (dim, dim)
        The covariance matrix.
    inv_transform: bool
        If True, use inverse transform instead of Box-Muller.
    engine: QMCEngine
        Quasi-Monte Carlo engine sampler. If None, Sobol' is used.
    seed : {int or `numpy.random.RandomState` instance}, optional
        If `seed` is not specified the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``RandomState`` instance, then that
        instance is used.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from scipy.stats import qmc
    >>> engine = qmc.MultivariateNormalQMC(mean=[0, 5], cov=[[1, 0], [0, 1]])
    >>> sample = engine.random(512)
    >>> plt.scatter(sample[:, 0], sample[:, 1])
    >>> plt.show()

    """

    def __init__(self, mean, cov, inv_transform=True, engine=None, seed=None):
        # check for square/symmetric cov matrix and mean vector has the same d
        mean = np.array(mean, copy=False, ndmin=1)
        cov = np.array(cov, copy=False, ndmin=2)
        if not mean.shape[0] == cov.shape[0]:
            raise ValueError("Dimension mismatch between mean and covariance.")
        if not np.allclose(cov, cov.transpose()):
            raise ValueError("Covariance matrix is not symmetric.")

        super().__init__(dim=mean.shape[0])
        self._mean = mean
        self._normal_engine = NormalQMC(
            dim=self.dim, inv_transform=inv_transform,
            engine=engine, seed=seed
        )
        # compute Cholesky decomp; if it fails, do the eigen decomposition
        try:
            self._corr_matrix = np.linalg.cholesky(cov).transpose()
        except np.linalg.LinAlgError:
            eigval, eigvec = np.linalg.eigh(cov)
            if not np.all(eigval >= -1.0e-8):
                raise ValueError("Covariance matrix not PSD.")
            eigval = np.clip(eigval, 0.0, None)
            self._corr_matrix = (eigvec * np.sqrt(eigval)).transpose()

    def random(self, n_samples: int = 1) -> np.ndarray:
        """Draw `n_samples` QMC samples from the multivariate Normal.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, dim)
            Sample.

        """
        base_samples = self._normal_engine.random(n_samples)
        qmc_samples = base_samples @ self._corr_matrix + self._mean
        return qmc_samples
