from __future__ import annotations

import inspect
from dataclasses import dataclass
from typing import (
    Callable, Dict, List, Literal, Protocol, TYPE_CHECKING, Tuple
)

import numpy as np

from scipy.stats._common import ConfidenceInterval
from scipy.stats._qmc import check_random_state
from scipy.stats._resampling import BootstrapResult
from scipy.stats import qmc, bootstrap


if TYPE_CHECKING:
    import numpy.typing as npt
    from scipy._lib._util import DecimalNumber, IntNumber, SeedType


__all__ = [
    'sobol_indices',
    'f_ishigami'
]


def f_ishigami(x: npt.ArrayLike) -> np.ndarray:
    r"""Ishigami function.

    .. math::

        Y(\mathbf{x}) = \sin x_1 + 7 \sin^2 x_2 + 0.1 x_3^4 \sin x_1

    with :math:`\mathbf{x} \in [-\pi, \pi]^3`.

    Parameters
    ----------
    x : array_like ([x1, x2, x3], n)

    Returns
    -------
    f : array_like (n,)
        Function evaluation.

    References
    ----------
    .. [1] Ishigami, T. and T. Homma. "An importance quantification technique
       in uncertainty analysis for computer models." IEEE,
       :doi:`10.1109/ISUMA.1990.151285`, 1990.
    """
    x = np.atleast_2d(x)
    f_eval = (
        np.sin(x[0])
        + 7 * np.sin(x[1])**2
        + 0.1 * (x[2]**4) * np.sin(x[0])
    )
    return f_eval


def sample_A_B(
    n: IntNumber,
    dists: List[PPFDist],
    random_state: SeedType = None
) -> Tuple[np.ndarray, np.ndarray]:
    """Sample two matrices A and B.

    Uses a Sobol' sequence with 2`d` columns to have 2 uncorrelated matrices.
    This is more efficient than using 2 random draw of Sobol'.
    See sec. 5 from [1]_.

    Output shape is (d, n).

    References
    ----------
    .. [1] Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and
       S. Tarantola. "Variance based sensitivity analysis of model
       output. Design and estimator for the total sensitivity index."
       Computer Physics Communications, 181(2):259-270,
       :doi:`10.1016/j.cpc.2009.09.018`, 2010.
    """
    d = len(dists)
    A_B = qmc.Sobol(d=2*d, seed=random_state, bits=64).random(n)

    A, B = A_B[:, :d], A_B[:, d:]

    for d_, dist in enumerate(dists):
        try:
            A[:, d_] = dist.ppf(A[:, d_])
            B[:, d_] = dist.ppf(B[:, d_])
        except AttributeError as exc:
            message = (
                "The method `ppf` must be specified for all distributions of "
                "'dists'."
            )
            raise ValueError(message) from exc

    return A.T, B.T


def sample_AB(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """AB matrix.

    AB: columns of B into A. Shape (d, n*d). e.g in 2d
    Take A and replace 1st column with 1st column of B. You have (d, n)
    Then A and replace 2nd column with 2nd column of B. You have (d, n)
    Concatenate to get AB. Which means AB is (d, 2*n).
    """
    d, n = A.shape
    AB = np.tile(A, (d, 1, 1))
    i = np.arange(d)
    AB[i, i] = B[i]
    AB = np.concatenate(AB, axis=1)
    return AB


def saltelli_2010(
    f_A: np.ndarray, f_B: np.ndarray, f_AB: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Saltelli2010 formulation.

    .. math::

        S_i = \frac{1}{N} \sum_{j=1}^N
        f(\mathbf{B})_j (f(\mathbf{AB}^{(i)})_j - f(\mathbf{A})_j)

    .. math::

        S_{T_i} = \frac{1}{N} \sum_{j=1}^N
        (f(\mathbf{A})_j - f(\mathbf{AB}^{(i)})_j)^2

    Parameters
    ----------
    f_A, f_B, f_AB : array_like (s, n)
        Function evaluations ``n`` being the number of evaluations and ``s``
        a vector output.

    Returns
    -------
    s, st : array_like (s, d)
        First order and total order Sobol' indices.

    References
    ----------
    .. [1] Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and
       S. Tarantola. "Variance based sensitivity analysis of model
       output. Design and estimator for the total sensitivity index."
       Computer Physics Communications, 181(2):259-270,
       :doi:`10.1016/j.cpc.2009.09.018`, 2010.
    """
    s_, n = f_A.shape
    f_AB = f_AB.reshape((-1, s_, n))

    # Empirical variance calculated using output from A and B which are
    # independent. Output of AB is not independent and cannot be used
    var = np.var(np.concatenate([f_A, f_B], axis=1), axis=1)

    # We divide by the variance to have a ratio of variance
    # this leads to eq. 2
    s = np.mean(f_B * (f_AB - f_A), axis=-1) / var  # Table 2 (b)
    st = 0.5 * np.mean((f_A - f_AB) ** 2, axis=-1) / var  # Table 2 (f)

    return s.reshape(s_, -1), st.reshape(s_, -1)


@dataclass
class BootstrapSobolResult:
    first_order: BootstrapResult
    total_order: BootstrapResult


@dataclass
class SobolResult:
    first_order: np.ndarray
    total_order: np.ndarray
    _indices_method: Callable
    _f_A: np.ndarray
    _f_B: np.ndarray
    _f_AB: np.ndarray
    _A: np.ndarray | None = None
    _B: np.ndarray | None = None
    _AB: np.ndarray | None = None
    _bootstrap_result: BootstrapResult = None  # type: ignore[valid-type]

    def bootstrap(
        self,
        confidence_level: DecimalNumber = 0.95,
        n_resamples: IntNumber = 99
    ) -> BootstrapResult:  # type: ignore[valid-type]
        """Bootstrap Sobol' indices to provide confidence intervals.

        Parameters
        ----------
        confidence_level : float, default: ``0.95``
            The confidence level of the confidence intervals.
        n_resamples : int, default: ``99``
            The number of resamples performed to form the bootstrap
            distribution of the indices.

        Returns
        -------
        res : BootstrapSobolResult
            Bootstrap result containing the confidence intervals and the
            bootstrap distribution of the indices.

            An object with attributes:

            first_order : BootstrapResult
                Bootstrap result of the first order indices.
            total_order : BootstrapResult
                Bootstrap result of the total order indices.
            See `BootstrapResult` for more details.

        """
        def statistic(idx):
            f_A_ = self._f_A[:, idx]
            f_B_ = self._f_B[:, idx]

            s, _ = self._f_A.shape
            f_AB_ = self._f_AB.reshape((s, -1, n))
            f_AB_ = f_AB_[:, :, idx].reshape(s, -1)

            return self._indices_method(f_A_, f_B_, f_AB_)

        n = self._f_A.shape[1]

        res = bootstrap(
            [np.arange(n)], statistic=statistic, method="BCa",
            n_resamples=n_resamples,
            confidence_level=confidence_level,
            bootstrap_result=self._bootstrap_result
        )
        self._bootstrap_result = res

        first_order = BootstrapResult(
            confidence_interval=ConfidenceInterval(
                res.confidence_interval.low[0], res.confidence_interval.high[0]
            ),
            bootstrap_distribution=res.bootstrap_distribution[0],
            standard_error=res.standard_error[0],
        )
        total_order = BootstrapResult(
            confidence_interval=ConfidenceInterval(
                res.confidence_interval.low[1], res.confidence_interval.high[1]
            ),
            bootstrap_distribution=res.bootstrap_distribution[1],
            standard_error=res.standard_error[1],
        )

        return BootstrapSobolResult(
            first_order=first_order, total_order=total_order
        )


class PPFDist(Protocol):
    @property
    def ppf(self) -> Callable[..., float]:
        ...


def sobol_indices(
    *,
    func: Callable[[np.ndarray], npt.ArrayLike] |
          Dict[Literal['f_A', 'f_B', 'f_AB'], np.ndarray],  # noqa
    n: IntNumber,
    dists: List[PPFDist] | None = None,
    method: Callable | Literal['saltelli_2010'] = 'saltelli_2010',
    random_state: SeedType = None
) -> SobolResult:
    r"""Global sensitivity indices of Sobol'.

    Parameters
    ----------
    func : callable or dict(str, array_like)
        If `func` is a callable, function to compute the Sobol' indices from.
        It's signature must be::

            func(x: ArrayLike) -> ArrayLike

        with ``x`` of shape ``(d, n)`` and the output should have a shape
        ``(s, n)`` with ``s`` the number of outputs and ``n`` the number of
        samples.

        If `func` is a dictionnary, contains the function evaluations from 3
        different arrays. Keys must be: ``f_A``, ``f_B`` and ``f_AB``.
        ``f_A`` and ``f_B`` should have a shape ``(s, n)`` and ``f_AB``
        should have a shape ``(s, d*n)``.
    n : int
        Number of samples.
        Must be a power of 2. The total number of function calls will be
        ``n(d+2)``.
    dists : list(distributions), optional
        List of each parameter's marginal distribution. Each parameter being
        independently distributed.

        Distributions must be an instance of a class with a ``ppf``
        method.

        Must be specified if `func` is a callable, and ignored otherwise.
    method : Callable or str, default: 'saltelli_2010'
        Method used to compute the first and total Sobol' indices.

        If a callable, it's signature must be::

            func(f_A: np.ndarray, f_B: np.ndarray, f_AB: np.ndarray)
            -> Tuple[np.ndarray, np.ndarray]

        with ``f_A, f_B, f_AB`` of shape (s, n) and the output being a tuple
        of the first and total indices with shape (s, d).
    random_state : {None, int, `numpy.random.Generator`}, optional
        If `random_state` is an int or None, a new `numpy.random.Generator` is
        created using ``np.random.default_rng(random_state)``.
        If `random_state` is already a ``Generator`` instance, then the
        provided instance is used.

    Returns
    -------
    res : SobolResult
        An object with attributes:

        first_order : ndarray of shape (s, d)
            First order Sobol' indices.
        total_order : ndarray of shape (s, d)
            Total order Sobol' indices.

        And methods:

        bootstrap(confidence_level: float, n_resamples: int)
        -> BootstrapSobolResult

            A method providing confidence intervals on the indices.
            See `scipy.stats.bootstrap` for more details.

            The bootstrapping is done on both first and total order indices
            and they are available in `BootstrapSobolResult` as attributes
            ``first_order`` and ``total_order``.

    Notes
    -----
    Sobol' method [1]_, [2]_ is a variance-based Sensitivity Analysis which
    allows to obtain the contribution of the parameters on the variance of the
    quantities of interest (QoIs; i.e., the outputs of `func`).
    Respective contributions can then be used to rank the parameters.

    .. note::

        Parameters are assumed to be independently distributed. Each
        parameter can still follow any distribution. In fact, the distribution
        is very important and should match the real distribution of the
        parameters.

    It uses a functional decomposition of the variance of the function to
    explore

    .. math::

        \mathbb{V}(Y) &= \sum_{i}^{d} \mathbb{V}_i (Y) + \sum_{i<j}^{d}
        \mathbb{V}_{ij}(Y) + ... + \mathbb{V}_{1,2,...,d}(Y),

    introducing conditional variances:

    .. math::

        \mathbb{V}_i(Y) = \mathbb{\mathbb{V}}[\mathbb{E}(Y|x_i)]
        \qquad
        \mathbb{V}_{ij}(Y) = \mathbb{\mathbb{V}}[\mathbb{E}(Y|x_i x_j)]
        - \mathbb{V}_i(Y) - \mathbb{V}_j(Y),

    Sobol' indices are expressed as

    .. math::

        S_i = \frac{\mathbb{V}_i(Y)}{\mathbb{V}[Y]}
        \qquad
        S_{ij} =\frac{\mathbb{V}_{ij}(Y)}{\mathbb{V}[Y]}.

    :math:`S_{i}` corresponds to the first-order term which apprises the
    contribution of the i-th parameter, while :math:`S_{ij}` corresponds to the
    second-order term which informs about the correlations between the
    i-th and the j-th parameters. These equations can be generalized to compute
    higher order terms; however, they are expensive to compute and their
    interpretation is complex.
    This is why only first order indices are provided.

    Total indices represent the global contribution of the parameters on the
    variance of the QoI and are defined as:

    .. math::

        S_{T_i} = S_i + \sum_j S_{ij} + \sum_{j,k} S_{ijk} + ...
        = 1 - \frac{\mathbb{V}[\mathbb{E}(Y|x_{\sim i})]}{\mathbb{V}[Y]}.

    First oder indices sum to 1, while the sum of total order indices will be
    greater than 1 if there are interactions.

    .. warning::

        Negative Sobol' values are due to numerical errors. Increasing the
        number of sample should help.

        If the parameters are not independent, the first-order indices would
        not sum to 1.
        Numerical noise can also contribute to this.

    References
    ----------
    .. [1] Sobol, I. M.. "Sensitivity analysis for nonlinear mathematical
       models." Mathematical Modeling and Computational Experiment, 1:407-414,
       1993.
    .. [2] Sobol, I. M. (2001). "Global sensitivity indices for nonlinear
       mathematical models and their Monte Carlo estimates." Mathematics
       and Computers in Simulation, 55(1-3):271-280,
       :doi:`10.1016/S0378-4754(00)00270-6`, 2001.
    .. [3] Saltelli, A. "Making best use of model evaluations to
       compute sensitivity indices."  Computer Physics Communications,
       145(2):280-297, :doi:`10.1016/S0010-4655(02)00280-1`, 2002.
    .. [4] Saltelli, A., M. Ratto, T. Andres, F. Campolongo, J. Cariboni,
       D. Gatelli, M. Saisana, and S. Tarantola. "Global Sensitivity Analysis.
       The Primer." 2007.
    .. [5] Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and
       S. Tarantola. "Variance based sensitivity analysis of model
       output. Design and estimator for the total sensitivity index."
       Computer Physics Communications, 181(2):259-270,
       :doi:`10.1016/j.cpc.2009.09.018`, 2010.
    .. [6] Ishigami, T. and T. Homma. "An importance quantification technique
       in uncertainty analysis for computer models." IEEE,
       :doi:`10.1109/ISUMA.1990.151285`, 1990.

    Examples
    --------
    The following is an example with the Ishigami function [6]_

    .. math::

        Y(\mathbf{x}) = \sin x_1 + 7 \sin^2 x_2 + 0.1 x_3^4 \sin x_1,

    with :math:`\mathbf{x} \in [-\pi, \pi]^3`. This function exhibits strong
    non-linearity and non-monotonicity.

    Remember, Sobol' indices assumes that samples are independently
    distributed. In this case we use a uniform distribution on each marginals.

    >>> import numpy as  np
    >>> from scipy.stats import f_ishigami, sobol_indices, uniform
    >>> rng = np.random.default_rng()
    >>> indices = sobol_indices(
    ...     func=f_ishigami, n=1024,
    ...     dists=[
    ...         uniform(loc=-np.pi, scale=2*np.pi),
    ...         uniform(loc=-np.pi, scale=2*np.pi),
    ...         uniform(loc=-np.pi, scale=2*np.pi)
    ...     ],
    ...     random_state=rng
    ... )
    >>> indices.first_order
    array([[0.31499073, 0.44011056, 0.00167054]])
    >>> indices.total_order
    array([[0.55508078, 0.43995732, 0.23803014]])

    .. note::

        By default, `scipy.stats.uniform` has support ``[0, 1]``.
        Using the parameters ``loc`` and ``scale``, one obtains the uniform
        distribution on ``[loc, loc + scale]``.

    This result is particularly interesting because the first order index
    :math:`S_{x_3} = 0` whereas its total order is :math:`S_{T_{x_3}} = 0.244`.
    This means that higher order interactions with :math:`x_3` are responsible
    for the difference. Almost 25% of the observed variance
    on the QoI is due to the correlations between :math:`x_3` and :math:`x_1`,
    although :math:`x_3` by itself has no impact on the QoI.

    The following gives a visual explanation of Sobol' indices on this
    function. Let's generate 1024 samples in :math:`[-\pi, \pi]^3` and
    calculate the value of the output.

    >>> from scipy.stats import qmc
    >>> n_dim = 3
    >>> p_labels = ['$x_1$', '$x_2$', '$x_3$']
    >>> sample = qmc.Sobol(d=n_dim, seed=rng).random(1024)
    >>> sample = qmc.scale(
    ...     sample=sample,
    ...     l_bounds=[-np.pi, -np.pi, -np.pi],
    ...     u_bounds=[np.pi, np.pi, np.pi]
    ... )
    >>> output = f_ishigami(sample)

    Now we can do scatter plots of the output with respect to each parameter.
    This gives a visual way to understand how each parameter impact the
    output of the function.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(1, n_dim)
    >>> for i in range(n_dim):
    ...     xi = sample[:, i]
    ...     ax[i].scatter(xi, output, marker='+')
    ...     ax[i].set_xlabel(p_labels[i])
    >>> ax[0].set_ylabel('Y')
    >>> plt.tight_layout()
    >>> plt.show()

    Now Sobol' goes a step further:
    by conditioning the output value by given values of the parameter
    (black lines), the conditional output mean is computed. It corresponds to
    the term :math:`\mathbb{E}(Y|x_i)`. Taking the variance of this term gives
    the numerator of the Sobol' indices.

    >>> mini = np.min(output)
    >>> maxi = np.max(output)
    >>> n_bins = 10
    >>> bins = np.linspace(-np.pi, np.pi, num=n_bins, endpoint=False)
    >>> dx = bins[1] - bins[0]
    >>> fig, ax = plt.subplots(1, n_dim)
    >>> for i in range(n_dim):
    ...     xi = sample[:, i]
    ...     ax[i].scatter(xi, output, marker='+')
    ...     ax[i].set_xlabel(p_labels[i])
    ...     for bin_ in bins:
    ...         idx = np.where((bin_ <= xi) & (xi <= bin_ + dx))
    ...         xi_ = xi[idx]
    ...         y_ = output[idx]
    ...         ave_y_ = np.mean(y_)
    ...         ax[i].plot([bin_ + dx/2] * 2, [mini, maxi], c='k')
    ...         ax[i].scatter(bin_ + dx/2, ave_y_, c='r')
    >>> ax[0].set_ylabel('Y')
    >>> plt.tight_layout()
    >>> plt.show()

    Looking at :math:`x_3`, the variance
    of the mean is zero leading to :math:`S_{x_3} = 0`. But we can further
    observe that the variance of the output is not constant along the parameter
    values of :math:`x_3`. This heteroscedasticity is explained by higher order
    interactions. Moreover, an heteroscedasticity is also noticeable on
    :math:`x_1` leading to an interaction between :math:`x_3` and :math:`x_1`.
    On :math:`x_2`, the variance seems to be constant and thus null interaction
    with this parameter can be supposed.

    This case is fairly simple to analyse visually---although it is only a
    qualitative analysis. Nevertheless, when the number of input parameters
    increases such analysis becomes unrealistic as it would be difficult to
    conclude on high-order terms. Hence the benefit of using Sobol' indices.

    """
    random_state = check_random_state(random_state)

    n_ = int(n)
    if not (n_ & (n_ - 1) == 0) or n != n_:
        raise ValueError(
            "The balance properties of Sobol' points require 'n' "
            "to be a power of 2."
        )
    n = n_

    if not callable(method):
        indices_methods: Dict[str, Callable] = {
            "saltelli_2010": saltelli_2010,
        }
        try:
            method = method.lower()  # type: ignore[assignment]
            indices_method = indices_methods[method]
        except KeyError as exc:
            message = (
                f"{method!r} is not a valid 'method'. It must be one of"
                f" {set(indices_methods)!r} or a callable."
            )
            raise ValueError(message) from exc
    else:
        kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
        parameters = (
            inspect.Parameter('f_A', kind=kind),
            inspect.Parameter('f_B', kind=kind),
            inspect.Parameter('f_AB', kind=kind)
        )
        sig = inspect.Signature(parameters=parameters)
        indices_method = method

        if inspect.signature(indices_method) != sig:
            message = (
                "If 'method' is a callable, it must have the following"
                f" signature: {inspect.signature(saltelli_2010)}"
            )
            raise ValueError(message)

    if callable(func):
        if dists is None:
            raise ValueError(
                "'dists' must be defined when 'func' is a callable."
            )

        A, B = sample_A_B(n=n, dists=dists, random_state=random_state)
        AB = sample_AB(A=A, B=B)

        f_A = np.atleast_2d(func(A))

        if f_A.shape[1] != n or f_A.shape[0] < 1:
            raise ValueError(
                "'func' output should have a shape ``(s, -1)`` with ``s`` "
                "the number of output."
            )

        f_B = np.atleast_2d(func(B))
        f_AB = np.atleast_2d(func(AB))
    else:
        message = (
            "When 'func' is a dictionary, it must contain the following "
            "keys: 'f_A', 'f_B' and 'f_AB'."
            "'f_A' and 'f_B' should have a shape ``(s, n)`` and 'f_AB' "
            "should have a shape ``(s, d*n)``."
        )
        try:
            f_A = np.atleast_2d(func['f_A'])
            f_B = np.atleast_2d(func['f_B'])
            f_AB = np.atleast_2d(func['f_AB'])
        except KeyError as exc:
            raise ValueError(message) from exc

        if f_A.shape[1] != n or f_A.shape != f_B.shape or \
                f_AB.shape == f_A.shape or f_AB.shape[1] % n != 0:
            raise ValueError(message)

    # Compute indices
    first_order, total_order = indices_method(f_A=f_A, f_B=f_B, f_AB=f_AB)

    res = dict(
        first_order=first_order,
        total_order=total_order,
        _indices_method=indices_method,
        _f_A=f_A,
        _f_B=f_B,
        _f_AB=f_AB
    )

    if callable(func):
        res.update(
            dict(
                _A=A,
                _B=B,
                _AB=AB,
            )
        )

    return SobolResult(**res)
