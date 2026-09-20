"""
Multivariate Distance Correlation, Distance Covariance, and Energy Distance.

References
----------
.. [1] G. J. Székely, M. L. Rizzo, and N. K. Bakirov, "Measuring and testing
       dependence by correlation of distances", The Annals of Statistics,
       Vol. 35, No. 6, pp. 2769-2794, 2007.
.. [2] G. J. Székely and M. L. Rizzo, "Partial distance correlation with
       methods for dissimilarities", The Annals of Statistics, Vol. 42,
       No. 6, pp. 2382-2412, 2014.
.. [3] M. L. Rizzo and G. J. Székely, "Energy distance", Wiley
       Interdisciplinary Reviews: Computational Statistics, Vol. 8, No. 1,
       pp. 27-38, 2016.
"""

import math
import numpy as np
from typing import NamedTuple, Optional, Union, Literal

from scipy._lib._array_api import xp_capabilities
from scipy._lib._util import check_random_state, MapWrapper, rng_integers
from scipy._lib._bunch import _make_tuple_bunch
from scipy.spatial.distance import pdist, cdist, squareform


__all__ = [
    'distance_covariance',
    'distance_correlation',
    'distance_covariance_test',
    'energy_distance_nd',
    'DistanceCovarianceResult',
]


DistanceCovarianceResult = _make_tuple_bunch(
    'DistanceCovarianceResult',
    ['statistic', 'pvalue'],
    ['null_distribution']
)


def _validate_arrays(x, y):
    """Validate sample arrays x and y for distance covariance/correlation."""
    x = np.asarray(x)
    y = np.asarray(y)

    if not (np.issubdtype(x.dtype, np.number) and np.issubdtype(y.dtype, np.number)):
        raise ValueError("Inputs x and y must contain numeric data.")

    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)

    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("Inputs x and y must contain only finite values (no NaN or inf).")

    if x.ndim == 1:
        x = x[:, np.newaxis]
    elif x.ndim != 2:
        raise ValueError(f"x must be a 1D or 2D array, got ndim={x.ndim}.")

    if y.ndim == 1:
        y = y[:, np.newaxis]
    elif y.ndim != 2:
        raise ValueError(f"y must be a 1D or 2D array, got ndim={y.ndim}.")

    n = x.shape[0]
    if y.shape[0] != n:
        raise ValueError(
            f"x and y must have the same number of observations (rows), "
            f"got {n} and {y.shape[0]}."
        )

    if n < 2:
        raise ValueError(
            f"At least 2 observations are required to compute distance covariance, "
            f"got {n}."
        )

    return x, y


def _double_center(D):
    """
    Perform classical double-centering on a distance matrix D:
    A_{ij} = D_{ij} - \bar{D}_{i.} - \bar{D}_{.j} + \bar{D}_{..}
    """
    row_mean = np.mean(D, axis=1, keepdims=True)
    col_mean = np.mean(D, axis=0, keepdims=True)
    grand_mean = np.mean(D)
    return D - row_mean - col_mean + grand_mean


def _u_center(D):
    r"""
    Perform U-centering on a distance matrix D for unbiased distance covariance
    (Székely & Rizzo 2014):
    A^*_{ij} = D_{ij} - \frac{1}{n-2} \sum_{l=1}^n D_{il}
                      - \frac{1}{n-2} \sum_{k=1}^n D_{kj}
                      + \frac{1}{(n-1)(n-2)} \sum_{k,l=1}^n D_{kl},  for i != j
    A^*_{ii} = 0
    """
    n = D.shape[0]
    if n < 4:
        raise ValueError(
            "Unbiased distance covariance requires at least 4 observations, "
            f"got {n}."
        )
    row_sum = np.sum(D, axis=1, keepdims=True)
    col_sum = np.sum(D, axis=0, keepdims=True)
    grand_sum = np.sum(D)

    A_star = (
        D
        - (row_sum / (n - 2))
        - (col_sum / (n - 2))
        + (grand_sum / ((n - 1) * (n - 2)))
    )
    np.fill_diagonal(A_star, 0.0)
    return A_star


@xp_capabilities(np_only=True)
def distance_covariance(
    x,
    y,
    *,
    method: Literal["auto", "biased", "unbiased"] = "auto",
    metric: str = "euclidean"
) -> float:
    r"""Compute the sample distance covariance between two sets of observations.

    Distance covariance :math:`\text{dCov}(X, Y)` is a measure of dependence
    between two random vectors :math:`X` and :math:`Y` of arbitrary dimensions.
    Importantly, :math:`\text{dCov}(X, Y) \ge 0`, and
    :math:`\text{dCov}(X, Y) = 0` if and only if :math:`X` and :math:`Y`
    are statistically independent.

    Parameters
    ----------
    x : array_like
        Sample array of shape ``(n,)`` or ``(n, p)`` representing :math:`n`
        observations of a :math:`p`-dimensional random variable.
    y : array_like
        Sample array of shape ``(n,)`` or ``(n, q)`` representing :math:`n`
        observations of a :math:`q`-dimensional random variable.
    method : {'auto', 'biased', 'unbiased'}, optional
        Estimator type for sample distance covariance:
        - ``'auto'`` / ``'biased'``: Classical sample distance covariance
          introduced in [1]_. Non-negative for all samples.
        - ``'unbiased'``: Unbiased estimator of population squared distance
          covariance introduced in [2]_. Can yield slightly negative values
          for finite independent samples, with expected value 0 under
          independence. Requires :math:`n \ge 4`.
    metric : str, optional
        The distance metric to use when computing pairwise distances between
        points. Default is ``'euclidean'``. Supported metrics include any
        valid metric string for :func:`scipy.spatial.distance.pdist`.

    Returns
    -------
    dcov : float
        The sample distance covariance. For ``method='biased'``, this is
        :math:`\text{dCov}_n(X, Y) = \sqrt{\text{dCov}_n^2(X, Y)} \ge 0`.
        For ``method='unbiased'``, this is the unbiased estimate of
        :math:`\text{dCov}^2(X, Y)`.

    See Also
    --------
    distance_correlation : Sample distance correlation.
    distance_covariance_test : Permutation test of independence.
    energy_distance_nd : Multivariate energy distance between two samples.

    Notes
    -----
    Given pairwise distance matrices :math:`D^X = (a_{ij})` and
    :math:`D^Y = (b_{ij})`, the classical double-centered matrices
    :math:`A` and :math:`B` have elements:

    .. math::

        A_{ij} = a_{ij} - \bar{a}_{i\cdot} - \bar{a}_{\cdot j} + \bar{a}_{\cdot\cdot}

    The classical sample squared distance covariance is:

    .. math::

        \text{dCov}_n^2(X, Y) = \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n A_{ij} B_{ij}

    References
    ----------
    .. [1] G. J. Székely, M. L. Rizzo, and N. K. Bakirov, "Measuring and testing
           dependence by correlation of distances", The Annals of Statistics,
           Vol. 35, No. 6, pp. 2769-2794, 2007.
    .. [2] G. J. Székely and M. L. Rizzo, "Partial distance correlation with
           methods for dissimilarities", The Annals of Statistics, Vol. 42,
           No. 6, pp. 2382-2412, 2014.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> rng = np.random.default_rng(42)

    Compute distance covariance for independent standard normals:

    >>> x = rng.normal(size=100)
    >>> y = rng.normal(size=100)
    >>> stats.distance_covariance(x, y)
    0.1192...

    Compute distance covariance for a non-linear deterministic relationship
    (e.g., :math:`Y = X^2`):

    >>> x = rng.uniform(-1, 1, size=200)
    >>> y = x ** 2
    >>> stats.distance_covariance(x, y)
    0.1652...
    """
    x, y = _validate_arrays(x, y)
    n = x.shape[0]

    if method not in ("auto", "biased", "unbiased"):
        raise ValueError(
            f"method must be one of 'auto', 'biased', or 'unbiased', got {method!r}."
        )

    D_x = squareform(pdist(x, metric=metric))
    D_y = squareform(pdist(y, metric=metric))

    if method == "unbiased":
        A = _u_center(D_x)
        B = _u_center(D_y)
        dcov2 = np.sum(A * B) / (n * (n - 3))
        return float(dcov2)
    else:
        A = _double_center(D_x)
        B = _double_center(D_y)
        dcov2 = np.sum(A * B) / (n * n)
        return float(np.sqrt(max(0.0, dcov2)))


@xp_capabilities(np_only=True)
def distance_correlation(
    x,
    y,
    *,
    method: Literal["auto", "biased", "unbiased"] = "auto",
    metric: str = "euclidean"
) -> float:
    r"""Compute the sample distance correlation between two sets of observations.

    Distance correlation :math:`\text{dCor}(X, Y)` is a scale- and
    dimension-free coefficient of multivariate dependence satisfying:

    .. math::

        0 \le \text{dCor}(X, Y) \le 1

    Crucially, :math:`\text{dCor}(X, Y) = 0` if and only if :math:`X` and
    :math:`Y` are statistically independent [1]_. This contrasts sharply with
    Pearson's correlation coefficient, which only detects linear associations
    and can vanish even for completely deterministic non-linear relationships.

    Parameters
    ----------
    x : array_like
        Sample array of shape ``(n,)`` or ``(n, p)`` representing :math:`n`
        observations of a :math:`p`-dimensional random variable.
    y : array_like
        Sample array of shape ``(n,)`` or ``(n, q)`` representing :math:`n`
        observations of a :math:`q`-dimensional random variable.
    method : {'auto', 'biased', 'unbiased'}, optional
        Estimator type for distance correlation:
        - ``'auto'`` / ``'biased'``: Classical sample distance correlation,
          guaranteed in :math:`[0, 1]`.
        - ``'unbiased'``: Modified distance correlation using unbiased U-centered
          matrices [2]_. Can yield small negative values near 0 for independent
          finite samples.
    metric : str, optional
        The distance metric to use when computing pairwise distances. Default
        is ``'euclidean'``.

    Returns
    -------
    dcor : float
        The sample distance correlation between :math:`x` and :math:`y`.

    See Also
    --------
    distance_covariance : Sample distance covariance.
    distance_covariance_test : Permutation hypothesis test of independence.

    Notes
    -----
    Sample distance correlation is defined as:

    .. math::

        \text{dCor}_n(X, Y) =
        \begin{cases}
            \frac{\text{dCov}_n(X, Y)}{\sqrt{\text{dVar}_n(X) \text{dVar}_n(Y)}},
            & \text{if } \text{dVar}_n(X)\text{dVar}_n(Y) > 0 \\
            0, & \text{if } \text{dVar}_n(X)\text{dVar}_n(Y) = 0
        \end{cases}

    where distance variance is :math:`\text{dVar}_n(X) = \text{dCov}_n(X, X)`.

    References
    ----------
    .. [1] G. J. Székely, M. L. Rizzo, and N. K. Bakirov, "Measuring and testing
           dependence by correlation of distances", The Annals of Statistics,
           Vol. 35, No. 6, pp. 2769-2794, 2007.
    .. [2] G. J. Székely and M. L. Rizzo, "Partial distance correlation with
           methods for dissimilarities", The Annals of Statistics, Vol. 42,
           No. 6, pp. 2382-2412, 2014.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> rng = np.random.default_rng(42)

    1. **Linear dependence**:
    >>> x = rng.normal(size=100)
    >>> y = 2 * x + 1
    >>> stats.distance_correlation(x, y)
    1.0

    2. **Non-linear dependence (Circle: :math:`X^2 + Y^2 = 1`)**:
    Pearson correlation is near 0, but distance correlation detects strong
    dependence:

    >>> theta = rng.uniform(0, 2 * np.pi, size=300)
    >>> x = np.cos(theta)
    >>> y = np.sin(theta)
    >>> np.corrcoef(x, y)[0, 1]  # Pearson
    -0.003...
    >>> stats.distance_correlation(x, y)
    0.234...
    """
    x, y = _validate_arrays(x, y)
    n = x.shape[0]

    if method not in ("auto", "biased", "unbiased"):
        raise ValueError(
            f"method must be one of 'auto', 'biased', or 'unbiased', got {method!r}."
        )

    D_x = squareform(pdist(x, metric=metric))
    D_y = squareform(pdist(y, metric=metric))

    if method == "unbiased":
        A = _u_center(D_x)
        B = _u_center(D_y)
        v_x = np.sum(A * A) / (n * (n - 3))
        v_y = np.sum(B * B) / (n * (n - 3))
        cov_xy = np.sum(A * B) / (n * (n - 3))
        denom = v_x * v_y
        if denom <= 0:
            return 0.0
        return float(cov_xy / np.sqrt(denom))
    else:
        A = _double_center(D_x)
        B = _double_center(D_y)
        v_x = np.sum(A * A) / (n * n)
        v_y = np.sum(B * B) / (n * n)
        denom = v_x * v_y
        if denom <= 0:
            return 0.0
        cov_xy2 = max(0.0, np.sum(A * B) / (n * n))
        dcor = np.sqrt(cov_xy2) / np.sqrt(np.sqrt(denom))
        return float(np.clip(dcor, 0.0, 1.0))


class _ParallelPermutation:
    """Helper callable for parallelizing permutation replications."""

    def __init__(self, A, B, random_states):
        self.A = A
        self.B = B
        self.random_states = random_states
        self.n = A.shape[0]

    def __call__(self, idx):
        perm = self.random_states[idx].permutation(self.n)
        # Fast permutation via row and column indexing of pre-centered matrix
        B_perm = self.B[perm, :][:, perm]
        return float(np.sum(self.A * B_perm) / self.n)


@xp_capabilities(np_only=True)
def distance_covariance_test(
    x,
    y,
    *,
    permutations: int = 999,
    metric: str = "euclidean",
    workers: int = 1,
    random_state: Optional[Union[int, np.random.RandomState, np.random.Generator]] = None
) -> DistanceCovarianceResult:
    r"""Nonparametric permutation test of independence based on distance covariance.

    Tests the null hypothesis:

    .. math::

        H_0: X \text{ and } Y \text{ are independent}

    against the general alternative of arbitrary (linear or non-linear) dependence.

    Parameters
    ----------
    x : array_like
        Sample array of shape ``(n,)`` or ``(n, p)`` representing :math:`n`
        observations of a :math:`p`-dimensional random variable.
    y : array_like
        Sample array of shape ``(n,)`` or ``(n, q)`` representing :math:`n`
        observations of a :math:`q`-dimensional random variable.
    permutations : int, optional
        Number of random permutations used to approximate the null distribution.
        Default is 999.
    metric : str, optional
        The distance metric to use when computing pairwise distances. Default
        is ``'euclidean'``.
    workers : int or map-like callable, optional
        If `workers` is an int, the population of permutations is evaluated in
        parallel across `workers` processes. Default is 1 (single process).
    random_state : {None, int, `numpy.random.Generator`, `numpy.random.RandomState`}, optional
        Pseudorandom number generator state used for reproducible permutations.

    Returns
    -------
    res : DistanceCovarianceResult
        An object with attributes:

        statistic : float
            The observed sample distance covariance :math:`\text{dCov}_n(X, Y)`.
        pvalue : float
            The permutation test :math:`p`-value, computed as:
            :math:`p = \frac{1 + \sum_{b=1}^B \mathbb{I}(T_b \ge T_{\text{obs}})}{1 + B}`.
        null_distribution : ndarray
            1D array of shape ``(permutations,)`` containing the permuted test
            statistics :math:`n \cdot \text{dCov}_n^2(X, Y^{(b)})` under :math:`H_0`.

    See Also
    --------
    distance_covariance : Compute sample distance covariance.
    distance_correlation : Compute sample distance correlation.

    References
    ----------
    .. [1] G. J. Székely, M. L. Rizzo, and N. K. Bakirov, "Measuring and testing
           dependence by correlation of distances", The Annals of Statistics,
           Vol. 35, No. 6, pp. 2769-2794, 2007.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> rng = np.random.default_rng(42)

    Test independence between independent bivariate uniforms:

    >>> x = rng.uniform(size=(60, 2))
    >>> y = rng.uniform(size=(60, 2))
    >>> res = stats.distance_covariance_test(x, y, permutations=499, random_state=42)
    >>> res.pvalue > 0.05
    True

    Test independence for non-linearly dependent variables (:math:`Y = X_1^2 + X_2^2`):

    >>> y = np.sum(x ** 2, axis=1)
    >>> res = stats.distance_covariance_test(x, y, permutations=499, random_state=42)
    >>> res.pvalue < 0.01
    True
    """
    x, y = _validate_arrays(x, y)
    n = x.shape[0]

    if permutations < 1:
        raise ValueError(f"permutations must be at least 1, got {permutations}.")

    D_x = squareform(pdist(x, metric=metric))
    D_y = squareform(pdist(y, metric=metric))

    A = _double_center(D_x)
    B = _double_center(D_y)

    stat_obs = float(np.sum(A * B) / n)
    dcov_obs = float(np.sqrt(max(0.0, stat_obs / n)))

    # Set up random states for parallel or sequential execution
    rng = check_random_state(random_state)
    seeds = [
        np.random.RandomState(rng_integers(rng, 1 << 32, size=4, dtype=np.uint32))
        for _ in range(permutations)
    ]

    worker_task = _ParallelPermutation(A=A, B=B, random_states=seeds)
    with MapWrapper(workers) as mapwrapper:
        null_dist = np.fromiter(
            mapwrapper(worker_task, range(permutations)),
            dtype=np.float64,
            count=permutations
        )

    count_greater = int(np.sum(null_dist >= stat_obs))
    pvalue = float((1.0 + count_greater) / (1.0 + permutations))

    return DistanceCovarianceResult(
        statistic=dcov_obs,
        pvalue=pvalue,
        null_distribution=null_dist
    )


@xp_capabilities(np_only=True)
def energy_distance_nd(
    u,
    v,
    *,
    metric: str = "euclidean"
) -> float:
    r"""Compute the energy distance between two multivariate empirical distributions.

    Energy distance is a statistical distance between the distributions of
    random vectors in arbitrary dimensions :math:`\mathbb{R}^d`. For two
    independent random vectors :math:`U, U' \sim F` and :math:`V, V' \sim G`,
    the energy distance satisfies :math:`\mathcal{E}(F, G) \ge 0`, with
    :math:`\mathcal{E}(F, G) = 0` if and only if :math:`F = G` [1]_.

    Parameters
    ----------
    u : array_like
        Sample array of shape ``(n,)`` or ``(n, d)`` representing :math:`n`
        observations in :math:`d` dimensions from distribution :math:`F`.
    v : array_like
        Sample array of shape ``(m,)`` or ``(m, d)`` representing :math:`m`
        observations in :math:`d` dimensions from distribution :math:`G`.
    metric : str, optional
        The distance metric to use when computing pairwise distances. Default
        is ``'euclidean'``.

    Returns
    -------
    dist : float
        The sample energy distance :math:`\mathcal{E}_{n,m}(u, v)`.

    See Also
    --------
    energy_distance : Energy distance between two 1D distributions.
    distance_covariance : Distance covariance between two sets of variables.

    Notes
    -----
    The empirical energy distance between samples :math:`u = \{u_1, \ldots, u_n\}`
    and :math:`v = \{v_1, \ldots, v_m\}` is:

    .. math::

        \mathcal{E}_{n,m}(u, v) = \frac{2}{nm} \sum_{i=1}^n \sum_{j=1}^m \|u_i - v_j\|
                                - \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n \|u_i - u_j\|
                                - \frac{1}{m^2} \sum_{i=1}^m \sum_{j=1}^m \|v_i - v_j\|

    References
    ----------
    .. [1] G. J. Székely and M. L. Rizzo, "Energy distance", Wiley
           Interdisciplinary Reviews: Computational Statistics, Vol. 8, No. 1,
           pp. 27-38, 2016.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> rng = np.random.default_rng(42)

    Energy distance between two samples from the same 3D distribution is near 0:

    >>> u = rng.normal(size=(100, 3))
    >>> v = rng.normal(size=(100, 3))
    >>> stats.energy_distance_nd(u, v)
    0.041...

    Energy distance between different distributions is strictly positive:

    >>> w = rng.normal(loc=2.0, size=(100, 3))
    >>> stats.energy_distance_nd(u, w)
    3.297...
    """
    u = np.asarray(u)
    v = np.asarray(v)

    if not (np.issubdtype(u.dtype, np.number) and np.issubdtype(v.dtype, np.number)):
        raise ValueError("Inputs u and v must contain numeric data.")

    u = np.ascontiguousarray(u, dtype=np.float64)
    v = np.ascontiguousarray(v, dtype=np.float64)

    if not (np.all(np.isfinite(u)) and np.all(np.isfinite(v))):
        raise ValueError("Inputs u and v must contain only finite values.")

    if u.ndim == 1:
        u = u[:, np.newaxis]
    elif u.ndim != 2:
        raise ValueError(f"u must be 1D or 2D, got ndim={u.ndim}.")

    if v.ndim == 1:
        v = v[:, np.newaxis]
    elif v.ndim != 2:
        raise ValueError(f"v must be 1D or 2D, got ndim={v.ndim}.")

    if u.shape[1] != v.shape[1]:
        raise ValueError(
            f"u and v must have the same dimensionality (columns), "
            f"got {u.shape[1]} and {v.shape[1]}."
        )

    n = u.shape[0]
    m = v.shape[0]
    if n < 1 or m < 1:
        raise ValueError("Samples must have at least 1 observation.")

    d_uv = cdist(u, v, metric=metric)
    d_uu = squareform(pdist(u, metric=metric)) if n > 1 else np.zeros((1, 1))
    d_vv = squareform(pdist(v, metric=metric)) if m > 1 else np.zeros((1, 1))

    e_dist = (
        (2.0 / (n * m)) * np.sum(d_uv)
        - (1.0 / (n * n)) * np.sum(d_uu)
        - (1.0 / (m * m)) * np.sum(d_vv)
    )
    return float(max(0.0, e_dist))
