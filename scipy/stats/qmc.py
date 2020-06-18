"""Quasi-Monte Carlo methods.

Define function to generate sample of points in the unit hypercube.

"""
from abc import ABC, abstractmethod
import copy
import numpy as np
from scipy.optimize import brute
from scipy._lib._util import check_random_state
from scipy.optimize import basinhopping
from scipy.stats._sobol import (
    SobolEngine, NormalQMCEngine, MultivariateNormalQMCEngine
)

__all__ = ['discrepancy', 'Halton', 'OrthogonalLatinHypercube',
           'LatinHypercube', 'OptimalDesign', 'SobolEngine',
           'NormalQMCEngine', 'MultivariateNormalQMCEngine']


def scale(sample, bounds, reverse=False):
    """Sample scaling from unit hypercube to bounds range.

    Parameters
    ----------
    sample : array_like (n_samples, k_vars)
        Sample to scale.
    bounds : tuple or array_like ([min, k_vars], [max, k_vars])
        Desired range of transformed data. The transformation apply the bounds
        on the sample and not the theoretical space, unit cube. Thus min and
        max values of the sample will coincide with the bounds.
    reverse : bool
        Reverse the transformation, from bounds range to unit hypercube.

    Returns
    -------
    sample : array_like (n_samples, k_vars)
        Scaled-sample.

    """
    bounds = np.asarray(bounds)
    min_ = np.min(bounds, axis=0)
    max_ = np.max(bounds, axis=0)
    if not reverse:
        return sample * (max_ - min_) + min_
    else:
        return (sample - min_) / (max_ - min_)


def discrepancy(sample, bounds=None, iterative=False):
    """Centered discrepancy on a given sample.

    The discrepancy is a uniformity criterion used to assess the space filling
    of a number of samples in a hypercube.
    The discrepancy measures how the spread of the points deviates from a
    uniform distribution.
    The lower the value is, the better the coverage of the parameter space is.

    Parameters
    ----------
    sample : array_like (n_samples, k_vars)
        The sample to compute the discrepancy from.
    bounds : tuple or array_like ([min, k_vars], [max, k_vars])
        Desired range of transformed data. The transformation applies the bounds
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
      Computer Science and Data Analysis Series, 2006.

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


def _update_discrepancy(x_new, sample, initial_disc, bounds=None):
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
        Desired range of transformed data. The transformation applies the bounds
        on the sample and not the theoretical space, unit cube. Thus min and
        max values of the sample will coincide with the bounds.

    Returns
    -------
    discrepancy : float
        Centered discrepancy of the sample composed of `x_new` and `sample`.

    """
    sample = np.asarray(sample)
    x_new = np.asarray(x_new)

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


def _perturb_discrepancy(sample, i1, i2, k, disc, bounds=None):
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


class QMCEngine(ABC):
    """Quasi-Monte Carlo engine sampler.

    Samples are distributed over the half-open interval [0, 1).

    Examples
    --------
    Generate samples from a low discrepancy sequence of Halton.

    >>> from scipy.stats import qmc
    >>> sampler = qmc.Halton(k_vars=2)
    >>> sample = sampler.random(n_samples=5)

    Compute the quality of the sample using the discrepancy criterion.

    >>> uniformity = qmc.discrepancy(sample)

    If some wants to continue an existing design, extra points can be obtained.

    >>> sampler.fast_forward(5)
    >>> sample_continued = sampler.random(n_samples=5)

    Finally, samples can be scaled to bounds.

    >>> bounds = np.array([[0, 2], [10, 5]])
    >>> sample_continued_scaled = qmc.scale(sample_continued, bounds)

    """

    @abstractmethod
    def __init__(self, k_vars, seed=None):
        """Engine initialization.

        Parameters
        ----------
        k_vars : int
            Dimension of the parameter space.
        seed : int or `np.random.RandomState`, optional
            If `seed` is not specified the `np.RandomState` singleton is used.
            If `seed` is an int, a new `np.random.RandomState` instance is used,
            seeded with seed.
            If `seed` is already a `np.random.RandomState instance`, then that
            `np.random.RandomState` instance is used.
            Specify `seed` for repeatable sampling.
        """
        self.k_vars = k_vars
        self.rng = check_random_state(seed)

    @abstractmethod
    def random(self, n_samples):
        """Draw n_samples in the half-open interval [0, 1).

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, k_vars)
            QMC sample.
        """
        pass

    def fast_forward(self, n):
        """Fast-forward the sequence by n positions.

        Parameters
        ----------
        n: int
            Number of points to skip in the sequence.
        """
        pass


class Halton(QMCEngine):
    """Halton sequence.

    Pseudo-random number generator that generalize the Van der Corput sequence
    for multiple dimensions. Halton sequence use base-two Van der Corput
    sequence for the first dimension, base-three for its second and base-n for
    its n-dimension.

    References
    ----------
    [1] Halton, "On the efficiency of certain quasi-random sequences of points
      in evaluating multi-dimensional integrals", Numerische Mathematik, 1960.

    Examples
    --------
    Generate samples from a low discrepancy sequence of Halton.

    >>> from scipy.stats import qmc
    >>> sampler = qmc.Halton(k_vars=2)
    >>> sample = sampler.random(n_samples=5)

    Compute the quality of the sample using the discrepancy criterion.

    >>> uniformity = qmc.discrepancy(sample)

    If some wants to continue an existing design, extra points can be obtained.

    >>> sampler.fast_forward(5)
    >>> sample_continued = sampler.random(n_samples=5)

    Finally, samples can be scaled to bounds.

    >>> bounds = np.array([[0, 2], [10, 5]])
    >>> sample_continued_scaled = qmc.scale(sample_continued, bounds)

    """

    def __init__(self, k_vars):
        """Engine initialization.

        Parameters
        ----------
        k_vars : int
            Dimension of the parameter space.

        """
        super().__init__(k_vars=k_vars)
        self.base = n_primes(k_vars)
        self.n_fast_forward = 0

    def random(self, n_samples):
        """Draw n_samples in the half-open interval [0, 1).

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, k_vars)
            QMC sample.
        """
        # Generate a sample using a Van der Corput sequence per dimension.
        sample = [van_der_corput(n_samples + 1, int(bdim), self.n_fast_forward)
                  for bdim in self.base]

        self.n_fast_forward += n_samples
        return np.array(sample).T[1:]

    def fast_forward(self, n):
        self.n_fast_forward += n


class OrthogonalLatinHypercube(QMCEngine):
    """Orthogonal array-based Latin hypercube sampling (OA-LHS).

    Samples are uniformly distributed over the half-open interval [low, high)
    (includes low, but excludes high).

    On top of the constraints from the Latin Hypercube, an orthogonal array of
    size n_samples is defined and only one point is allowed per subspace.

    References
    ----------
    [1] Art B. Owen, "Orthogonal arrays for computer experiments, integration
    and visualization", Statistica Sinica, 1992.

    """

    def __init__(self, k_vars, seed=None):
        """Engine initialization.

        Parameters
        ----------
        k_vars : int
            Dimension of the parameter space.
        seed : int or `np.random.RandomState`, optional
            If `seed` is not specified the `np.RandomState` singleton is used.
            If `seed` is an int, a new `np.random.RandomState` instance is used,
            seeded with seed.
            If `seed` is already a `np.random.RandomState instance`, then that
            `np.random.RandomState` instance is used.
            Specify `seed` for repeatable sampling.

        """
        super().__init__(k_vars=k_vars, seed=seed)

    def random(self, n_samples):
        """Draw n_samples in the half-open interval [0, 1).

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, k_vars)
            QMC sample.
        """
        sample = []
        step = 1.0 / n_samples

        for _ in range(self.k_vars):
            # Enforce a unique point per grid
            j = np.arange(n_samples) * step
            temp = j + self.rng.uniform(low=0, high=step, size=n_samples)
            self.rng.shuffle(temp)

            sample.append(temp)

        return np.array(sample).T


class LatinHypercube(QMCEngine):
    """Latin hypercube sampling (LHS).

    Samples are uniformly distributed over the half-open interval [low, high)
    (includes low, but excludes high).

    The parameter space is subdivided into an orthogonal grid of n_samples per
    dimension. Within this multi-dimensional grid, n_samples are selected by
    ensuring there is only one sample per row and column.

    References
    ----------
    [1] Mckay et al., "A Comparison of Three Methods for Selecting Values of
    Input Variables in the Analysis of Output from a Computer Code",
    Technometrics, 1979.

    """

    def __init__(self, k_vars, centered=False, seed=None):
        """Engine initialization.

        Parameters
        ----------
        k_vars : int
            Dimension of the parameter space.
        centered : bool
            Center the point within the multi-dimensional grid.
        seed : int or `np.random.RandomState`, optional
            If `seed` is not specified the `np.RandomState` singleton is used.
            If `seed` is an int, a new `np.random.RandomState` instance is used,
            seeded with seed.
            If `seed` is already a `np.random.RandomState instance`, then that
            `np.random.RandomState` instance is used.
            Specify `seed` for repeatable sampling.

        """
        super().__init__(k_vars=k_vars, seed=seed)
        self.centered = centered

    def random(self, n_samples):
        if self.centered:
            r = 0.5
        else:
            r = self.rng.random_sample((n_samples, self.k_vars))

        q = self.rng.randint(low=1, high=n_samples,
                             size=(n_samples, self.k_vars))

        return 1. / n_samples * (q - r)


class OptimalDesign(QMCEngine):
    """Optimal design.

    Optimize the design by doing random permutations to lower the centered
    discrepancy. If `optimization` is False, `niter` design are generated and
    the one with lowest centered discrepancy is return. This option is faster.

    Centered discrepancy based design show better space filling robustness
    toward 2D and 3D subprojections. Distance based design better space filling
    but less robust to subprojections.

    References
    ----------
    [1] Damblin et al., "Numerical studies of space filling designs:
    optimization of Latin Hypercube Samples and subprojection properties",
    Journal of Simulation, 2013.

    """

    def __init__(self, k_vars, start_design=None, niter=1, force=False,
                 optimization=True, seed=None):
        """Engine initialization.

        Parameters
        ----------
        k_vars : int
            Dimension of the parameter space.
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

        """
        super().__init__(k_vars=k_vars, seed=seed)
        self.start_design = start_design
        self.niter = niter
        self.force = force
        self.optimization = optimization

        self.best_doe = self.start_design
        self.best_disc = np.inf

        self.olhs = OrthogonalLatinHypercube(self.k_vars, seed=self.rng)

    def random(self, n_samples):
        """Draw n_samples in the half-open interval [0, 1).

        Parameters
        ----------
        n_samples : int
            Number of samples to generate in the parameter space.

        Returns
        -------
        sample : array_like (n_samples, k_vars)
            QMC sample.
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
                doe = copy.deepcopy(self.best_doe)
                col, row_1, row_2 = np.round(x).astype(int)
                doe[row_1, col], doe[row_2, col] = doe[row_2, col], doe[row_1, col]

                disc = _perturb_discrepancy(self.best_doe, row_1, row_2, col,
                                            self.best_disc)

                if disc < self.best_disc:
                    self.best_disc = disc
                    self.best_doe = doe

                return disc

            # Total number of possible design
            complexity = self.k_vars * n_samples ** 2

            if (complexity > 1e6) or self.force:
                bounds_optim = ([0, self.k_vars - 1],
                                [0, n_samples - 1],
                                [0, n_samples - 1])
            else:
                bounds_optim = (slice(0, self.k_vars - 1, 1), slice(0, n_samples - 1, 1),
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
