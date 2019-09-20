"""Quasi-Monte Carlo methods.

Define function to generate sample of points in the unit hypercube.

"""

from __future__ import division

import copy
import numpy as np
from scipy.optimize import brute
from scipy._lib._util import check_random_state

try:
    from scipy.optimize import basinhopping

    have_basinhopping = True
except ImportError:
    have_basinhopping = False


def discrepancy(sample, bounds=None, iterative=False):
    """Discrepancy.

    Compute the centered discrepancy on a given sample.
    It is a measure of the uniformity of the points in the parameter space.
    The lower the value is, the better the coverage of the parameter space is.

    Parameters
    ----------
    sample : array_like (n_samples, k_vars)
        The sample to compute the discrepancy from.
    bounds : tuple or array_like ([min, k_vars], [max, k_vars])
        Desired range of transformed data. The transformation apply the bounds
        on the sample and not the theoretical space, unit cube. Thus min and
        max values of the sample will coincide with the bounds.
    iterative : bool
        Must be False if not using it for updating the discrepancy.

    Returns
    -------
    discrepancy : float
        Centered discrepancy.

    References
    ----------
    [1] Fang et al. "Design and modeling for computer experiments",
      Computer Science and Data Analysis Series Science and Data Analysis
      Series, 2006.

    """
    sample = np.asarray(sample)

    n_samples, dim = sample.shape

    if iterative:
        n_samples += 1

    # Sample scaling from bounds to unit hypercube
    if bounds is not None:
        min_ = bounds.min(axis=0)
        max_ = bounds.max(axis=0)
        sample = (sample - min_) / (max_ - min_)

    abs_ = abs(sample - 0.5)
    disc1 = np.sum(np.prod(1 + 0.5 * abs_ - 0.5 * abs_ ** 2, axis=1))

    prod_arr = 1
    for i in range(dim):
        s0 = sample[:, i]
        prod_arr *= (1 +
                     0.5 * abs(s0[:, None] - 0.5) + 0.5 * abs(s0 - 0.5) -
                     0.5 * abs(s0[:, None] - s0))
    disc2 = prod_arr.sum()

    c2 = ((13.0 / 12.0) ** dim - 2.0 / n_samples * disc1 +
          1.0 / (n_samples ** 2) * disc2)

    return c2


def update_discrepancy(x_new, sample, initial_disc, bounds=None):
    """Update the discrepancy with a new sample.

    Parameters
    ----------
    x_new : array_like (1, k_vars)
        The new sample to add in `sample`.
    sample : array_like (n_samples, k_vars)
        The initial sample.
    initial_disc : float
        Centered discrepancy of the `sample`.
    bounds : tuple or array_like ([min, k_vars], [max, k_vars])
        Desired range of transformed data. The transformation apply the bounds
        on the sample and not the theoretical space, unit cube. Thus min and
        max values of the sample will coincide with the bounds.

    Returns
    -------
    discrepancy : float
        Centered discrepancy of the sample composed of `x_new` and `sample`.

    """
    # Sample scaling from bounds to unit hypercube
    if bounds is not None:
        min_ = bounds.min(axis=0)
        max_ = bounds.max(axis=0)
        sample = (sample - min_) / (max_ - min_)
        x_new = (x_new - min_) / (max_ - min_)

    n_samples = len(sample) + 1
    abs_ = abs(x_new - 0.5)

    disc1 = - 2 / n_samples * np.prod(1 + 1 / 2 * abs_ - 1 / 2 * abs_ ** 2)
    disc2 = 2 / (n_samples ** 2) * np.sum(np.prod(1 + 1 / 2 * abs_ +
                                                  1 / 2 * abs(sample - 0.5) -
                                                  1 / 2 * abs(x_new - sample),
                                                  axis=1))
    disc3 = 1 / (n_samples ** 2) * np.prod(1 + abs_)

    return initial_disc + disc1 + disc2 + disc3


def perturb_discrepancy(sample, i1, i2, k, disc, bounds=None):
    """Centered discrepancy after and elementary perturbation on a LHS.

    An elementary perturbation consists of an exchange of coordinates between
    two points: ``sample[i1, k] <-> sample[i2, k]``. By construction,
    this operation conserves the LHS properties.

    Parameters
    ----------
    sample : array_like (n_samples, k_vars)
        The sample (before permutation) to compute the discrepancy from.
    i1 : int
        The first line of the elementary permutation.
    i2 : int
        The second line of the elementary permutation.
    k : int
        The column of the elementary permutation.
    disc : float
        Centered discrepancy of the design before permutation.
    bounds : tuple or array_like ([min, k_vars], [max, k_vars])
        Desired range of transformed data. The transformation apply the bounds
        on the sample and not the theoretical space, unit cube. Thus min and
        max values of the sample will coincide with the bounds.

    Returns
    -------
    discrepancy : float
        Centered discrepancy.

    References
    ----------
    [1] Jin et al. "An efficient algorithm for constructing optimal design
        of computer experiments", Journal of Statistical Planning and
        Inference, 2005.

    """
    sample = np.asarray(sample)
    n_samples = sample.shape[0]

    # Sample scaling from bounds to unit hypercube
    if bounds is not None:
        min_ = bounds.min(axis=0)
        max_ = bounds.max(axis=0)
        sample = (sample - min_) / (max_ - min_)

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
    [1] `StackOverflow <https://stackoverflow.com/questions/2068372>`_.

    """
    sieve = np.ones(n // 3 + (n % 6 == 2), dtype=np.bool)
    for i in range(1, int(n ** 0.5) // 3 + 1):
        if sieve[i]:
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
        big_number = 10
        while 'Not enought primes':
            primes = primes_from_2_to(big_number)[:n]
            if len(primes) == n:
                break
            big_number += 1000

    return primes


def van_der_corput(n_samples, base=2, start_index=0):
    """Van der Corput sequence.

    Pseudo-random number generator based on a b-adic expansion.

    Parameters
    ----------
    n_samples : int
        Number of element of the sequence.
    base : int
        Base of the sequence.
    start_index : int
        Index to start the sequence from.

    Returns
    -------
    sequence : list (n_samples,)
        Sequence of Van der Corput.

    """
    sequence = []
    for i in range(start_index, start_index + n_samples):
        n_th_number, denom = 0., 1.
        quotient = i
        while quotient > 0:
            quotient, remainder = divmod(quotient, base)
            denom *= base
            n_th_number += remainder / denom
        sequence.append(n_th_number)

    return sequence


def halton(dim, n_samples, bounds=None, start_index=0):
    """Halton sequence.

    Pseudo-random number generator that generalize the Van der Corput sequence
    for multiple dimensions. Halton sequence use base-two Van der Corput
    sequence for the first dimension, base-three for its second and base-n for
    its n-dimension.

    Parameters
    ----------
    dim : int
        Dimension of the parameter space.
    n_samples : int
        Number of samples to generate in the parameter space.
    bounds : tuple or array_like ([min, k_vars], [max, k_vars])
        Desired range of transformed data. The transformation apply the bounds
        on the sample and not the theoretical space, unit cube. Thus min and
        max values of the sample will coincide with the bounds.
    start_index : int
        Index to start the sequence from.

    Returns
    -------
    sequence : array_like (n_samples, k_vars)
        Sequence of Halton.

    References
    ----------
    [1] Halton, "On the efficiency of certain quasi-random sequences of points
      in evaluating multi-dimensional integrals", Numerische Mathematik, 1960.

    Examples
    --------
    Generate samples from a low discrepancy sequence of Halton.

    >>> from scipy.stats import qmc
    >>> sample = qmc.halton(dim=2, n_samples=5)

    Compute the quality of the sample using the discrepancy criterion.

    >>> uniformity = qmc.discrepancy(sample)

    If some wants to continue an existing design, extra points can be obtained.

    >>> sample_continued = qmc.halton(dim=2, n_samples=5, start_index=5)

    """
    base = n_primes(dim)

    # Generate a sample using a Van der Corput sequence per dimension.
    sample = [van_der_corput(n_samples + 1, bdim, start_index) for bdim in base]
    sample = np.array(sample).T[1:]

    # Sample scaling from unit hypercube to feature range
    if bounds is not None:
        min_ = bounds.min(axis=0)
        max_ = bounds.max(axis=0)
        sample = sample * (max_ - min_) + min_

    return sample


global best_doe, best_disc


def orthogonal_latin_hypercube(dim, n_samples, bounds=None, seed=None):
    """Orthogonal array-based Latin hypercube sampling (OA-LHS).

    Samples are uniformly distributed over the half-open interval [low, high)
    (includes low, but excludes high).

    On top of the constraints from the Latin Hypercube, an orthogonal array of
    size n_samples is defined and only one point is allowed per subspace.

    Parameters
    ----------
    dim : int
        Dimension of the parameter space.
    n_samples : int
        Number of samples to generate in the parameter space.
    bounds : tuple or array_like ([min, k_vars], [max, k_vars])
        Desired range of transformed data. The transformation apply the bounds
        on the sample and not the theoretical space, unit cube. Thus min and
        max values of the sample will coincide with the bounds.
    seed : int or `np.random.RandomState`, optional
        If `seed` is not specified the `np.RandomState` singleton is used.
        If `seed` is an int, a new `np.random.RandomState` instance is used,
        seeded with seed.
        If `seed` is already a `np.random.RandomState instance`, then that
        `np.random.RandomState` instance is used.
        Specify `seed` for repeatable sampling.

    Returns
    -------
    sample : ndarray (n_samples, k_vars)
        Latin hypercube Sampling.

    References
    ----------
    [1] Art B. Owen, "Orthogonal arrays for computer experiments, integration
    and visualization", Statistica Sinica, 1992.

    """
    sample = []
    step = 1.0 / n_samples

    rng = check_random_state(seed)

    for _ in range(dim):
        # Enforce a unique point per grid
        j = np.arange(n_samples) * step
        temp = j + rng.uniform(low=0, high=step, size=n_samples)
        rng.shuffle(temp)

        sample.append(temp)

    sample = np.array(sample).T

    # Sample scaling from unit hypercube to feature range
    if bounds is not None:
        min_ = bounds.min(axis=0)
        max_ = bounds.max(axis=0)
        sample = sample * (max_ - min_) + min_

    return sample


def latin_hypercube(dim, n_samples, bounds=None, centered=False, seed=None):
    """Latin hypercube sampling (LHS).

    Samples are uniformly distributed over the half-open interval [low, high)
    (includes low, but excludes high).

    The parameter space is subdivided into an orthogonal grid of n_samples per
    dimension. Within this multi-dimensional grid, n_samples are selected by
    ensuring there is only one sample per row and column.

    Parameters
    ----------
    dim : int
        Dimension of the parameter space.
    n_samples : int
        Number of samples to generate in the parametr space.
    bounds : tuple or array_like ([min, k_vars], [max, k_vars])
        Desired range of transformed data. The transformation apply the bounds
        on the sample and not the theoretical space, unit cube. Thus min and
        max values of the sample will coincide with the bounds.
    centered : bool
        Center the point within the multi-dimensional grid.
    seed : int or `np.random.RandomState`, optional
        If `seed` is not specified the `np.RandomState` singleton is used.
        If `seed` is an int, a new `np.random.RandomState` instance is used,
        seeded with seed.
        If `seed` is already a `np.random.RandomState instance`, then that
        `np.random.RandomState` instance is used.
        Specify `seed` for repeatable sampling.

    Returns
    -------
    sample : ndarray (n_samples, k_vars)
        Latin hypercube Sampling.

    References
    ----------
    [1] Mckay et al., "A Comparison of Three Methods for Selecting Values of
    Input Variables in the Analysis of Output from a Computer Code",
    Technometrics, 1979.

    """
    rng = check_random_state(seed)
    if centered:
        r = 0.5
    else:
        r = rng.random_sample((n_samples, dim))

    q = rng.randint(low=1, high=n_samples, size=(n_samples, dim))

    sample = 1. / n_samples * (q - r)

    # Sample scaling from unit hypercube to feature range
    if bounds is not None:
        min_ = bounds.min(axis=0)
        max_ = bounds.max(axis=0)
        sample = sample * (max_ - min_) + min_

    return sample


def optimal_design(dim, n_samples, bounds=None, start_design=None, niter=1,
                   force=False, optimization=True, seed=None):
    """Optimal design.

    Optimize the design by doing random permutations to lower the centered
    discrepancy. If `optimization` is False, `niter` design are generated and
    the one with lowest centered discrepancy is return. This option is faster.

    Centered discrepancy based design show better space filling robustness
    toward 2D and 3D subprojections. Distance based design better space filling
    but less robust to subprojections.

    Parameters
    ----------
    dim : int
        Dimension of the parameter space.
    n_samples : int
        Number of samples to generate in the parametr space.
    bounds : tuple or array_like ([min, k_vars], [max, k_vars])
        Desired range of transformed data. The transformation apply the bounds
        on the sample and not the theoretical space, unit cube. Thus min and
        max values of the sample will coincide with the bounds.
    start_design : array_like (n_samples, k_vars)
        Initial design of experiment to optimize.
    niter : int
        Number of iteration to perform.
    force : bool
        If `optimization`, force *basinhopping* optimization. Otherwise
        grid search is used.
    optimization : bool
        Optimal design using global optimization or random generation of
        `niter` samples.
    seed : int or `np.random.RandomState`, optional
        If `seed` is not specified the `np.RandomState` singleton is used.
        If `seed` is an int, a new `np.random.RandomState` instance is used,
        seeded with seed.
        If `seed` is already a `np.random.RandomState instance`, then that
        `np.random.RandomState` instance is used.
        Specify `seed` for repeatable sampling.

    Returns
    -------
    sample : array_like (n_samples, k_vars)
        Optimal Latin hypercube Sampling.

    References
    ----------
    [1] Damblin et al., "Numerical studies of space filling designs:
    optimization of Latin Hypercube Samples and subprojection properties",
    Journal of Simulation, 2013.

    """
    global best_doe, best_disc
    best_doe = start_design
    best_disc = np.inf
    rng = check_random_state(seed)

    if (bounds is None) and (best_doe is not None):
        bounds = np.array([best_doe.min(axis=0), best_doe.max(axis=0)])
    if optimization:
        if best_doe is None:
            best_doe = orthogonal_latin_hypercube(dim, n_samples, bounds,
                                                  seed=rng)

        best_disc = discrepancy(best_doe, bounds)

        def _perturb_best_doe(x, bounds):
            """Perturbe the DoE and keep track of the best DoE.

            Parameters
            ----------
            x : list of int
                It is a list of:
                    idx : int
                        Index value of the components to compute
            bounds : tuple or array_like ([min, k_vars], [max, k_vars])
                Desired range of transformed data. The transformation apply the
                bounds on the sample and not the theoretical space, unit cube.
                Thus min and max values of the sample will coincide with the
                bounds.

            Returns
            -------
            discrepancy : float
                Centered discrepancy.

            """
            global best_doe, best_disc

            # Perturbe the DoE
            doe = copy.deepcopy(best_doe)
            col, row_1, row_2 = np.round(x).astype(int)
            doe[row_1, col], doe[row_2, col] = doe[row_2, col], doe[row_1, col]

            disc = perturb_discrepancy(best_doe, row_1, row_2, col,
                                       best_disc, bounds)

            if disc < best_disc:
                best_disc = disc
                best_doe = doe

            return disc

        # Total number of possible design
        complexity = dim * n_samples ** 2

        if have_basinhopping and ((complexity > 1e6) or force):
            bounds_optim = ([0, dim - 1],
                            [0, n_samples - 1],
                            [0, n_samples - 1])
        else:
            bounds_optim = (slice(0, dim - 1, 1), slice(0, n_samples - 1, 1),
                            slice(0, n_samples - 1, 1))

        for _ in range(niter):
            if have_basinhopping and ((complexity > 1e6) or force):
                minimizer_kwargs = {"method": "L-BFGS-B",
                                    "bounds": bounds_optim,
                                    "args": (bounds,)}
                _ = basinhopping(_perturb_best_doe, [0, 0, 0], niter=100,
                                 minimizer_kwargs=minimizer_kwargs)
            else:
                _ = brute(_perturb_best_doe, ranges=bounds_optim,
                          finish=None, args=(bounds,))

    else:
        for _ in range(niter):
            doe = orthogonal_latin_hypercube(dim, n_samples, bounds, seed=rng)
            disc = discrepancy(doe, bounds)
            if disc < best_disc:
                best_disc = disc
                best_doe = doe

    return best_doe
