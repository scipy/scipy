from __future__ import annotations
from typing import Callable, List, Optional, TYPE_CHECKING, Tuple

import numpy as np


from scipy.stats._qmc import check_random_state
from scipy.stats import qmc, bootstrap
from scipy.stats.sampling import NumericalInversePolynomial

if TYPE_CHECKING:
    import numpy.typing as npt
    from scipy.stats._unuran.unuran_wrapper import PINVDist
    from scipy._lib._util import IntNumber, SeedType


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
    x : array_like (n, [x1, x2, x3])

    Returns
    -------
    f : array_like (n, 1)
        Function evaluation.

    References
    ----------
    .. [1] Ishigami, T. and T. Homma. "An importance quantification technique
       in uncertainty analysis for computer models." IEEE,
       :doi:`10.1109/ISUMA.1990.151285`, 1990.
    """
    x = np.atleast_2d(x)
    f_eval = (
        np.sin(x[:, 0])
        + 7 * np.sin(x[:, 1])**2
        + 0.1 * (x[:, 2]**4) * np.sin(x[:, 0])
    )
    return f_eval.reshape(-1, 1)


def sample_A_B(
    d: int, n: int,
    dists: Optional[List[PINVDist]] = None,
    random_state: SeedType = None
) -> Tuple[np.ndarray,  np.ndarray]:
    """Sample two matrices A and B.

    Uses a Sobol' sequence with 2`d` columns to have 2 uncorrelated matrices.
    This is more efficient than using 2 random draw of Sobol'.
    See sec. 5 from [1]_.

    References
    ----------
    .. [1] Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and
       S. Tarantola. "Variance based sensitivity analysis of model
       output. Design and estimator for the total sensitivity index."
       Computer Physics Communications, 181(2):259-270,
       :doi:`10.1016/j.cpc.2009.09.018`, 2010.
    """
    A_B = qmc.Sobol(d=2*d, seed=random_state, bits=64).random(n)

    A, B = A_B[:, :d], A_B[:, d:]

    if dists is not None:
        for d_, dist in enumerate(dists):
            dist_rng = NumericalInversePolynomial(
                dist, random_state=random_state
            )
            A[:, d_] = dist_rng.ppf(A[:, d_])
            B[:, d_] = dist_rng.ppf(B[:, d_])

    return A, B


def sample_AB(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """AB matrix.

    AB: columns of B into A. Shape (n*d, d). e.g in 2d
    Take A and replace 1st column with 1st column of B. You have (n, d)
    Then A and replace 2nd column with 2nd column of B. You have (n, d)
    Concatenate to get AB. Which means AB is (2*n, d).
    """
    n, d = A.shape
    AB = np.empty((int(d*n), d))
    for i in range(d):
        AB[i*n:(i+1)*n, :] = np.column_stack((A[:, :i], B[:, i], A[:, i+1:]))

    return AB


def saltelli_2010(
    f_A: np.ndarray, f_B: np.ndarray, f_AB: np.ndarray
) -> Tuple[np.ndarray,  np.ndarray]:
    r"""Saltelli2010 formulation.

    .. math::

        S_i = \frac{1}{N} \sum_{j=1}^N
        f(\mathbf{B})_j (f(\mathbf{AB}^{(i)})_j - f(\mathbf{A})_j)

    .. math::

        S_{T_i} = \frac{1}{N} \sum_{j=1}^N
        (f(\mathbf{A})_j - f(\mathbf{AB}^{(i)})_j)^2

    References
    ----------
    .. [1] Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and
       S. Tarantola. "Variance based sensitivity analysis of model
       output. Design and estimator for the total sensitivity index."
       Computer Physics Communications, 181(2):259-270,
       :doi:`10.1016/j.cpc.2009.09.018`, 2010.
    """
    f_AB = f_AB.reshape(-1, *f_A.shape)

    var = np.var(np.vstack([f_A, f_B]), axis=0)

    s = np.mean(f_B * (f_AB - f_A), axis=1) / var
    st = 0.5 * np.mean((f_A - f_AB) ** 2, axis=1) / var

    return s, st


def sobol_indices(
    *,
    func: Callable,
    n: IntNumber,
    d: IntNumber,
    dists: Optional[List[PINVDist]] = None,
    l_bounds: Optional[npt.ArrayLike] = None,
    u_bounds: Optional[npt.ArrayLike] = None,
    random_state: SeedType = None
):
    r"""Global sensitivity indices of Sobol'.

    Parameters
    ----------
    func : callable
        Function to compute the Sobol' indices from. It's signature must be
        ``func(x: ArrayLike) -> ArrayLike``, with ``x`` of shape ``(-1, d)``
        and the output should have a shape ``(-1, s)`` with ``s`` the number
        of output.
    n : int
        Must be a power of 2. The total number of function call will
        ``n(d+2)``.
    d : int
        Dimension of the parameter space.
    dists : list(distributions), optional
        List of distributions of the parameters. It must be compatible with
        `scipy.stats.sampling.NumericalInversePolynomial`.
    l_bounds, u_bounds : array_like (d,), optional
        Lower and upper bounds (resp. :math:`a`, :math:`b`) of transformed
        data in the parameter space.
    random_state : {None, int, `numpy.random.Generator`}, optional
        If `random_state` is an int or None, a new `numpy.random.Generator` is
        created using ``np.random.default_rng(random_state)``.
        If `random_state` is already a ``Generator`` instance, then the
        provided instance is used.

    Notes
    -----
    Variance-based Sensitivity Analysis allows obtaining the contribution of
    theparameters on the quantity of interest’s (QoI) variance. Sobol'
    method [1]_, [2]_, gives not only a ranking but also quantifies the
    importance factor using the variance.

    .. note:: Parameters are assumed to be independant.

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
    higher order terms. However, the computational effort to converge them is
    most often not at hand and their analysis and interpretations are not
    simple. This is why only first order indices are provided.

    Total indices represent the global contribution of the parameters on the
    variance of the QoI and express as:

    .. math::

        S_{T_i} = S_i + \sum_j S_{ij} + \sum_{j,k} S_{ijk} + ...
        = 1 - \frac{\mathbb{V}[\mathbb{E}(Y|x_{\sim i})]}{\mathbb{V}[Y]}.

    The total number of function call is ``n+n+d*n=n(d+2)``.
    Three matrices are required for the computation of
    the indices: A, B and a permutation matrix AB based
    on both A and B.

    .. warning::

        Negative Sobol' values are due to numerical errors. Increasing the
        number of sample should help.

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

    >>> from scipy.stats import f_ishigami
    >>> indices = sobol_indices(
    ...    func=f_ishigami, n=1024, d=3,
    ...    l_bounds=[-np.pi, -np.pi, -np.pi], u_bounds=[np.pi, np.pi, np.pi]
    ... )
    >>> indices[:2]
    (array([0.28717719, 0.44694404, 0.00404842]),  # may vary
     array([0.53461512, 0.42775372, 0.23343319]))

    It is particularly interesting because the first order indice of
    :math:`S_{x_3} = 0` whereas its total order is :math:`S_{T_{x_3}} = 0.244`.
    It means that higher order interactions with :math:`x_3` are responsible
    for the difference. Almost 25% of the observed variance
    on the QoI is due to the correlations between :math:`x_3` and :math:`x_1`,
    although :math:`x_3` by itself has no impact on the QoI.

    The following gives a visual explanation of Sobol' indices. It shows
    scatter plots of the output with respect to each parameter. By conditioning
    the output value by given values of the parameter
    (black lines), the conditional output mean is computed. It corresponds to
    the term :math:`\mathbb{E}(Y|x_i)`. Taking the variance of this term gives
    the numerator of the Sobol' indices. Looking at :math:`x_3`, the variance
    of the mean is zero leading to :math:`S_{x_3} = 0`. But we can further
    observe that the variance of the output is not constant along the parameter
    values of :math:`x_3`. This heteroscedasticity is explained by higher order
    interactions. Moreover, an heteroscedasticity is also noticeable on
    :math:`x_1` leading to an interaction between :math:`x_3` and :math:`x_1`.
    On :math:`x_2`, the variance seems to be constant and thus null interaction
    with this parameter can be supposed. This case is fairly simple to analyse
    visually---although it is only a qualitative analysis. Nevertheless, when
    the number of input parameters increases such analysis becomes unrealistic
    as it would be difficult to conclude on high-order terms.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.stats import f_ishigami, qmc
    >>> n_dim = 3
    >>> p_labels = ['$x_1$', '$x_2$', '$x_3$']
    >>> sample = qmc.Sobol(d=n_dim).random(1024)
    >>> sample = qmc.scale(
    ...     sample=sample,
    ...     l_bounds=[-np.pi, -np.pi, -np.pi],
    ...     u_bounds=[np.pi, np.pi, np.pi]
    ... )
    >>> output = f_ishigami(sample)
    >>> mini = np.min(output)
    >>> maxi = np.max(output)
    >>> n_bins = 10
    >>> bins = np.linspace(-np.pi, np.pi, num=n_bins, endpoint=False)
    >>> dx = bins[1] - bins[0]
    >>> fig, ax = plt.subplots(1, n_dim)
    >>> for i in range(n_dim):
    >>>     xi = sample[:, i]
    >>>     ax[i].scatter(xi, output, marker='+')
    >>>     ax[i].set_xlabel(p_labels[i])
    >>>     for bin_ in bins:
    >>>         idx = np.where((bin_ <= xi) & (xi <= bin_ + dx))
    >>>         xi_ = xi[idx]
    >>>         y_ = output[idx]
    >>>         ave_y_ = np.mean(y_)
    >>>         ax[i].plot([bin_ + dx / 2] * 2, [mini, maxi], c='k')
    >>>         ax[i].scatter(bin_ + dx / 2, ave_y_, c='r')
    >>> ax[0].set_ylabel('Y')
    >>> plt.tight_layout()
    >>> plt.show()

    """
    random_state = check_random_state(random_state)
    A, B = sample_A_B(d=d, n=n, dists=dists, random_state=random_state)
    AB = sample_AB(A=A, B=B)

    if (l_bounds is not None) or (u_bounds is not None):
        A = qmc.scale(A, l_bounds=l_bounds, u_bounds=u_bounds)
        B = qmc.scale(B, l_bounds=l_bounds, u_bounds=u_bounds)
        AB = qmc.scale(AB, l_bounds=l_bounds, u_bounds=u_bounds)

    f_A = func(A)
    f_B = func(B)
    f_AB = func(AB)

    # Y = (Y - Y.mean()) / Y.std()

    s, st = saltelli_2010(f_A, f_B, f_AB)

    def saltelli_2010_(idx):
        f_A_ = f_A[idx]
        f_B_ = f_B[idx]
        f_AB_ = f_AB[idx]
        return saltelli_2010(f_A_, f_B_, f_AB_)

    ci = bootstrap(
        [np.arange(n)], saltelli_2010_, method="BCa",
        batch=int(n*0.7), n_resamples=99
    )

    return s, st, ci
