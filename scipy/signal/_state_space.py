"""Functions for working with state-space systems. """
import numbers
from typing import Literal

import numpy as np
import scipy._external.array_api_extra as xpx

from scipy import linalg
from scipy._lib._array_api import Array, array_namespace, xp_capabilities, \
    xp_promote, xp_result_device, xp_size
from ._filter_design import findfreqs
from ._lti_conversion import abcd_normalize

__all__ = ['freqs_ss']


def _diagsvd(s, m, n, xp=None):
    """Array API compatible version of `scipy.linalg.diagsvd`. """
    xp = xp or array_namespace(s)
    DD = linalg.diagsvd(np.asarray(s), m, n)
    return xp.asarray(DD, device=xp_result_device(s))


def _generalized_eigvals(a, b, overwrite_a=False, overwrite_b=False):
    """Array API compatible version of `scipy.linalg.eigvals` for calculating
    generalized homogeneous eigenvalues. """
    xp, device = array_namespace(a, b), xp_result_device(a, b)
    a_np, b_np = np.asarray(a), np.asarray(b)
    # pyrefly: ignore[deprecated]
    w = linalg.eigvals(a_np, b_np, overwrite_a=overwrite_a, overwrite_b=overwrite_b,
                       homogeneous_eigvals=True)
    return xp.asarray(w[0], device=device), xp.asarray(w[1], device=device)


def _geomspace(start, stop, num, xp, device):
    """Array API compatible version of `numpy.geomspace`. """
    a = np.geomspace(start, stop, num=num, dtype=np.float64)
    return xp.asarray(a, device=device)


def _schur_complex(A: Array, overwrite_a: bool = False) -> tuple[Array, Array]:
    """Array API compatible version of `scipy.linalg.schur` with complex output and
    without sorting. """
    xp, device = array_namespace(A), xp_result_device(A)
    T_np, Z_np = linalg.schur(np.asarray(A), output='complex', overwrite_a=overwrite_a,
                              sort=None)
    T, Z = xp.asarray(T_np, device=device), xp.asarray(Z_np, device=device)
    dtyp = xp.result_type(T, Z, 1j)  # We want allways the same return dtype
    return xp.astype(T, dtyp, copy=False), xp.astype(Z, dtyp, copy=False)


def _solve_triangular(a: Array, b: Array, lower: bool = False,
                      overwrite_b: bool = False) -> Array:
    """Array API compatible version of `scipy.linalg.solve_triangular`. """
    xp, device = array_namespace(a, b), xp_result_device(a, b)
    a_np, b_np = np.asarray(a), np.asarray(b)
    x_np = linalg.solve_triangular(a_np, b_np, trans=0, lower=lower,
                                   unit_diagonal=False, overwrite_b=overwrite_b,
                                   check_finite=True)
    return xp.asarray(x_np, device=device)


@xp_capabilities(cpu_only=True, skip_backends=[('dask.array', "Missing property `mT`")])
def conjT(M: Array) -> Array:
    """Conjugate transpose or Hermitian transpose of the matrix `M`.

    If `M` has more than two dimensions, the conjugate transpose is applied to the last
    two axes.
    """
    xp = array_namespace(M)
    return xp.conj(M.mT)


@xp_capabilities(cpu_only=True, skip_backends=[
    ('dask.array', "svd() does not support parameter `full_matrices=True`"),
    ('jax.numpy', "Using `s = xpx.at(s)[r:].set(0)` for `s[r:] = 0` fails.")])
def _svd_rank(A: Array, m_size: Literal['full', 'shape', 'rank'] = 'full', *,
              rtol: float = 1e-12, atol: float = 0.) -> tuple[Array, Array, Array, int]:
    """Calculate the SVD ``A = U @ diag(s) @ Vh``and the rank `r` of a matrix `A`.

    Parameters
    ----------
    A : array
        Input matrix of shape (m, n).
    m_size : {'full', 'shape', 'rank'}
        Determines the shape of the returned matrices:
        If 'full',  `U` and `Vh` are of shape (m, m) and (n, n).
        If 'shape', `U` and `Vh` are of shape (m, k) and (k, n), with k = min(m, n).
        If 'rank',  `U` and `Vh` are of shape (m, r) and (r, n), with `r` being the
        rank of `A`.
    rtol, atol : float
        Tolerance for considering the singular values of `A` to be zero.
        That is, ``s = 0 if abs(s) <= rtol * s_max + atol else s``, where `s_max` is the
        largest singular value of `A`.

    Returns
    -------
    U, s, Vh : array
       The singular value composition as three arrays of compatible shapes.
    r : int
        The rank of `A`, i.e., the number of non-zero singular values based on the given
        tolerance.
    """
    xp = array_namespace(A)
    # Note that svd() returns the singular values sorted in decreasing order.
    U, s, Vh = xp.linalg.svd(A, full_matrices=(m_size == 'full')) # A = U @ diag(s) @ Vh

    if xp_size(s):  # Search s for nearest value to tolerance:
        r = s.shape[0] - xp.searchsorted(xp.flip(s), rtol * s[0] + atol, side='right')
    else:
        r = 0

    if m_size == 'full' or m_size == 'shape':
        s[r:] = 0  # set all singular values below the threshold to zero
        return U, s, Vh, r
    elif m_size == 'rank':
        return U[:, :r], s[:r], Vh[:r, :], r
    else:
        raise ValueError(f"Parameter {m_size=} not in ('full', 'shape', 'rank')")

@xp_capabilities(cpu_only=True, skip_backends=[
    ('dask.array', "function _svd_rank() not supported."),
    ('jax.numpy',  "function _svd_rank() not supported.")])
def _tran_zeros_reduce(A: Array, B: Array, C: Array, D: Array, rtol: float = 1e-12,
                       atol: float = 0.) -> tuple[Array, Array, Array, Array]:
    """Reduce a system so that its feedthrough matrix `D` has full row rank.

    Parameters
    ----------
    A, B, C, D : array
        `A`, `B`, `C`, `D` are complex or real-valued arrays of
        shape (n, n), (n, p), (q, n), (q, p).
    rtol, atol : float
        Tolerance for considering a (positive) singular value to be zero; that is,
        ``s = 0 if s <= rtol * s_max + atol else s``, where `s_max` is the largest
        singular value of the matrix.

    Returns
    -------
    A_r, B_r, C_r, D_r : array
        Reduced system matrices of shape (n_r, n_r), (n_r, p_r), (q_r, n_r), and
        (q_r, p_r), respectively.

    Notes
    -----
    This function implements the system-matrix reduction algorithm presented in [1]_ as
    "REDUCE". It also draws on the internal function `_tzeros_reduce` in the `harold`_
    package, which implements the same algorithm.

    Note that this algorithm is numerically sensitive towards `atol`: E.g., for the unit
    test `TestTransmissionZeros.test_rectangular_D`, the entries of `C_tilde` have,
    depending on backend/architecture, a magnitude of < 1e-30 or are zero. Hence, the
    rank calculation is unstable for ``atol=0``.

    Currently, this helper is solely used by `_transmission_zeros()`.

    References
    ----------
    .. [1] Abbas Emami-Naeini, Paul M. van Dooren, "Computation of zeros of linear
           multivariable systems" Automatica 18.4 (1982): 0005-1098.
           :doi:`10.1016/0005-1098(82)90070-X`
    .. _harold: https://harold.readthedocs.io/
    """
    A, B, C, D = abcd_normalize(A, B, C, D)
    xp = array_namespace(A)  # abcd_normalize() already ensures a unified API backend.
    device = xp_result_device(A, B, C, D)
    nu_, (mu_, m) = A.shape[0], D.shape

    for _ in range(A.shape[0]):  # at most n iterations are needed.
        # Compress rows of D:
        if xp.any(D):  # D != 0:
            U_D, s_D, Vh_D, sigma_ = _svd_rank(D, 'full', rtol=rtol, atol=atol)
            if sigma_ == mu_:  # D has full row rank, so the reduction is complete.
                break
            D_bar = (_diagsvd(s_D, *D.shape) @ Vh_D)[:sigma_, :]
            C1 = conjT(U_D) @ C
        else:  # D == 0:
            C1, D_bar, sigma_ = C, xp.zeros((0, m), dtype=D.dtype, device=device), 0
        C_bar, C_tilde = C1[:sigma_, :], C1[sigma_:, :]   # partition C

        # Compress columns of C_tilde to [0, diag(s_C) @ Vh_C]:
        U_C, _, _, rho_ = _svd_rank(conjT(C_tilde), 'full', rtol=rtol, atol=atol)
        U_C = xp.flip(U_C, axis=1)
        nu_, mu_ = nu_ - rho_, rho_ + sigma_
        if nu_ == 0:  # C_tilde has full column rank → system has no zeros.
            return (xp.zeros((0, 0), dtype=A.dtype, device=device),) * 4
        if rho_ == 0:  # C_tilde is zero, so the reduction is complete.
            break
        Uh_C = conjT(U_C)
        ABCD = xp.concat([xp.concat([Uh_C @ A @ U_C, Uh_C @ B], axis=1),
                          xp.concat([   C_bar @ U_C,    D_bar], axis=1)], axis=0)
        A, B = ABCD[:nu_, :nu_], ABCD[:nu_, -m:]
        C, D = ABCD[nu_:, :nu_], ABCD[nu_:, -m:]
    return A, B, C, D

@xp_capabilities(cpu_only=True, skip_backends=[
    ('dask.array', "function _tran_zeros_reduce() not supported."),
    ('jax.numpy',  "function _tran_zeros_reduce() not supported.")])
def _transmission_zeros(A: Array, B: Array, C: Array, D: Array, inf_zeros: bool = True,
                        rtol: float = 1e-12, atol: float = 0.) -> Array:
    r"""Helper for calculating the transmission zeros of a state space system.

    Parameters
    ----------
    A, B, C, D : array
        `A`, `B`, `C`, `D` are complex or real-valued arrays of
        shape (n, n), (n, p), (q, n), (q, p).
    inf_zeros : bool
        If ``True``, also zeros at infinity, i.e., ``np.inf``, are returned.
    rtol, atol : float
        Tolerance for considering a (positive) singular value to be zero; that is,
        ``s = 0 if s <= rtol * s_max + atol else s``, where `s_max` is the largest
        singular value of the matrix. It is also used to determine whether the
        quotients of the homogeneous eigenvalues (the zeros) are finite or infinite.

    Returns
    -------
    t_zeros : array
        Complex-valued 1d array with transmission zeros sorted in no particular order.
        An empty array is returned if the system has no transmission zeros.

    Notes
    -----
    The transmission zeros are the values of :math:`s` at which the rank of the transfer
    function :math:`G(s) = C (sI - A)^{-1} B + D` drops below its normal rank.

    This function implements the system-matrix reduction algorithm presented in [1]_ as
    "ZEROS" and then determines the zeros by using generalized eigenvalues. This
    implementation is very similar to the harold_ function `transmission_zeros()` but
    not identical.

    References
    ----------
    .. [1] Abbas Emami-Naeini, Paul M. van Dooren, "Computation of zeros of linear
           multivariable systems" Automatica 18.4 (1982): 0005-1098.
           :doi:`10.1016/0005-1098(82)90070-X`
    .. _harold: https://harold.readthedocs.io/
    """
    A, B, C, D = abcd_normalize(A, B, C, D)
    xp = array_namespace(A)  # abcd_normalize() already ensures a unified API backend.
    device = xp_result_device(A, B, C, D)
    dtyp = xp.result_type(A, 1j)  # we want always the same return dtype

    # Zeros only exist if B != 0 and C != 0:
    if not (xp.any(B) and xp.any(C)):  # Using any() ensures B and C being non-empty
        return xp.asarray([], dtype=dtyp, device=device)

    # Reduce system to full row rank feedthrough matrix D:
    A, B, C, D = _tran_zeros_reduce(A, B, C, D, rtol=rtol, atol=atol)
    if A.shape[0] == 0:  # no zeros
        return xp.asarray([], dtype=dtyp, device=device)

    if D.shape[0] != D.shape[1]:  # If D not square, reduce to full column rank D:
        ACBD_conjT = _tran_zeros_reduce(conjT(A), conjT(C), conjT(B), conjT(D),
                                        rtol=rtol, atol=atol)
        A, C, B, D = (conjT(M_) for M_ in ACBD_conjT)
    if A.shape[0] == 0:  # no zeros
        return xp.asarray([], dtype=dtyp, device=device)

    # Find unitary W for which [C D] = [0 D_f] @ W holds:
    _, _, Wh = xp.linalg.svd(xp.concat((C, D), axis=1), full_matrices=True)
    W, n = conjT(xp.flip(Wh, axis=0)), A.shape[0],
    # The zeros are the generalized eigenvalues λ of A_f v = λ B_f v:
    A_f, B_f = xp.concat((A, B), axis=1) @ W[:, :n], W[:n, :n]
    alpha, beta = _generalized_eigvals(A_f, B_f, overwrite_a=True, overwrite_b=True)
    # The homogenous eigenvalues are used to be able to determine zeros at ∞:
    ii = abs(beta) > rtol * xp.max(xp.abs(beta)) + atol  # indices of non-∞ zeros
    if not inf_zeros:  # only return finite zeros:
        return alpha[ii] / beta[ii]
    t_zeros = xp.full_like(alpha, xp.inf)
    t_zeros[ii] = alpha[ii] / beta[ii]
    return t_zeros


@xp_capabilities(cpu_only=True, skip_backends=[
    ('jax.numpy', "Expression `xp.all(abs_lamsM.T > tol, axis=0)` fails."),
    ('dask.array', "Expression `xp.all(abs_lamsM.T > tol, axis=0)` fails."),])
def _calcSSTF_helper(T: Array, B: Array, C: Array, s_values: Array, rtol: float,
                     atol: float, overwrite_b: bool, xp, device) -> Array:
    """Helper for calculating G = C (s I - T)^{-1} B assuming T is upper triangular.

    If ``s I - T`` is invertible, `solve_triangular` is used to solve the corresponding
    linear system. Otherwise, a singular value decomposition is used. Zero-valued
    singular values are handled so that they propagate correctly to infinite entries of
    `G`.

    This helper is used only by `freqs_ss` and is factored out to ease unit testing.

    Parameters
    ----------
    T, B, C : array
        `T`, `B`, `C` are complex or real-valued arrays of shape (n, n), (n, p), (q, n).
        `T` must be an upper triangular matrix.
    s_values : array
        Array of shape (N,) providing the values for s.
    rtol, atol : float
        Tolerance for considering the eigenvalues and singular values of ``s I - T`` to
        be zero. For each value of `s`, values at or below ``rtol * abs(λ_max) + atol``
        are treated as zero, where ``λ_max`` is the eigenvalue of ``s I - T`` with the
        largest absolute value.
    overwrite_b : bool
        Toggle whether the parameter `B` is potentially overwritten while calculating
        the result.

    Returns
    -------
    G : array
        An array of shape (N, q, p) containing the results
        ``G[i] = C @ inv(s_values[i] * eye(n) - T) @ B``.

    Notes
    -----
    The parameters are first promoted to float64 / complex128 dtypes to ensure a
    consistent accuracy.

    Adding a finite value to an infinite value, or multiplying a nonzero finite value by
    an infinite value, yields an infinite result. However, ``inf - inf`` and
    ``0 * inf`` yield ``NaN``. Therefore, matrix multiplication needs special treatment.

    Working with (and returning) (N, q, p) arrays instead of (q, p, N) arrays,
    simplifies indexing in this function considerably (e.g., `xp.diagonal` always works
    on the last axes).
    """
    n, q, p, N = T.shape[0], C.shape[0], B.shape[1], s_values.shape[0]

    T, B, C, s_values = (xp.astype(X_, xp.result_type(X_, xp.float64), copy=False)
                         for X_ in (T, B, C, s_values))  # ensure float64/complex128

    # M = sI - T of shape (N, n, n):
    M = xp.reshape(s_values, (N, 1, 1)) * xp.eye(n, device=device) - T
    G = xp.zeros((N, q, p), dtype=M.dtype, device=device)  # return value

    abs_lamsM = xp.abs(xp.linalg.diagonal(M))  # (N, n) array of abs. eigenvalues
    tol = rtol * xp.max(abs_lamsM, axis=-1) + atol  # (N,) array of tolerances
    ii = xp.all(abs_lamsM.T > tol, axis=0)  # invertible matrices
    if xp.any(ii):  # Use triangular solver for invertible matrices (no Array API)
        M_ii = _solve_triangular(M[ii], B, lower=False, overwrite_b=overwrite_b)
        C, M_ii = (xp.astype(X_, G.dtype, copy=False) for X_ in (C, M_ii)) # for PyTorch
        G[ii] = C @ M_ii
    if not xp.all(ii):  # Use SVD for non-invertible matrices (no Array API):
        U, s, Vh = xp.linalg.svd(M[~ii])  # M[~ii] = U @ diag(s) @ Vh
        # for PyTorch:
        U, Vh, B, C = (xp.astype(X_, G.dtype, copy=False) for X_ in (U, Vh, B, C))
        V, Uh_B = xp.conj(Vh).mT, xp.conj(U).mT @ B

        # Split S^{-1} into a finite and an infinite part:
        inv_s, ii0_s = xp.zeros_like(s), (s.T > tol[~ii]).T
        inv_s[ii0_s] = 1 / s[ii0_s]  # finite
        V_iS = V * xp.reshape(inv_s, (-1, 1, n))  # V_iS = V @ S^{-1}
        G_ = C @ V_iS @ Uh_B  # the result of the finite part

        # Infinite part: Mark all ∞ values in S^{-1} as 1 and all others as 0.
        # Then propagate those marks to G:
        ii1_s = inv_s == 0
        inv_s[ii1_s], inv_s[~ii1_s] = 1, 0
        V_iS = V * xp.reshape(inv_s, (-1, 1, n))  # V_iS = V @ S^{-1}
        G_inf = C @ V_iS @ Uh_B  # values for G ≠ 0 are ∞

        G_[G_inf != 0] = xp.inf * (1+1j) # set both real and imaginary pat to ∞
        G[~ii] = G_
    return G


def freqs_ss(A, B, C, D, *, om_spec, f_scale='log', rtol=1e-12, atol=0.):
    r"""Frequency response of a continuous-time state space system.

    Parameters
    ----------
    A : array_like
        State matrix as a two-dimensional array of shape (n, n).
    B : array_like
        Input matrix as a two-dimensional array of shape (n, p).
    C : array_like
        Output matrix as a two-dimensional array of shape (q, n).
    D : array_like
        Feedthrough matrix as a two-dimensional array of shape (q, p).
    om_spec : int | (float, float, int) | array_like
        The angular frequency range specification `om_spec` may have one of three forms:

        If `om_spec` is an array, it is used directly as the angular frequency values.
        The array must be a real-valued one-dimensional array of at least length one.
        The parameter `f_scale` is ignored in this case.

        If `om_spec` is a tuple of the form ``(om0, om1, N)``, i.e., of type
        ``tuple[float, float, int]``, the angular frequencies range from ``om0`` to
        ``om1`` in ``N`` samples. The parameter `f_scale` selects linear or logarithmic
        spacing.

        If `om_spec` is a positive integer ``N``, then the angular frequencies are
        determined by passing the system poles and the system's transmission zeros to
        :func:`~scipy.signal.findfreqs`.
    f_scale : 'log' | 'linear'
        The frequency scaling can be either logarithmic (default) or linear.
        This parameter is ignored if `om_spec` is an array.
    rtol, atol : float
        Tolerance for considering an angular frequency equal to a system's pole (i.e.,
        `A`'s eigenvalue). For automatic range selection (i.e., ``om_spec = N``), it is
        also used for checking if a singular value equals is zero. Two values ``x1``,
        ``x2`` are considered equal if ``abs(x1 - x2) <= rtol * max(abs(x1), abs(x2)) +
        atol`` holds. The default values are ``rtol=1e-12`` and ``atol=0``.

    Returns
    -------
    om : array
        A float64 array of angular frequencies of shape (N,) for which the response
        is calculated.
        If `om_spec` is an array, ``om = np.astype(om_spec, np.float64, copy=False)``.
    G : array
        A complex128 array of shape (q, p, N) representing the frequency response.

    See Also
    --------
    freqs : Frequency response of a continuous-time transfer function.
    freqs_zpk : Frequency response of a continuous-time transfer function in
                zero-pole-gain form.
    freqz : Frequency response of a discrete-time transfer function.
    freqz_zpk : Frequency response of a discrete-time transfer function in
                zero-pole-gain form.
    freqz_sos : Frequency response of a discrete-time transfer function in
                second-order-section format.

    Notes
    -----
    :func:`~scipy.signal.abcd_normalize` is used to verify that `A`, `B`, `C`, and `D`
    are compatible. Furthermore, it is assumed that all entries of `A`, `B`, `C`, and
    `D` are finite. The numeric calculations are performed with `~numpy.float64` /
    `~numpy.complex128` dtypes for a consistent accuracy.

    The frequency response can be derived by applying the Laplace transform to the
    state space system, i.e.,

    .. math::
        % LaTeX Macros are valid for the remainder of the HTML page:
        \newcommand{\IC}{{\mathbb{C}}}  % set of complex numbers
        \newcommand{\conjT}[1]{\overline{#1^T}} % transposed complex conjugate
        \newcommand{\inv}[1]{\left(#1\right)^{\!-1}} % inverse, i.e. ()^{-1}
        \newcommand{\jj}{{\mathbb{j}}}  % imaginary unit
        %
        \newcommand{\vb}[1]{\mathbf{#1}} % vectors and matrices are bold:
        \newcommand{\mA}{\vb{A}}\newcommand{\mB}{\vb{B}}\newcommand{\mC}{\vb{C}}
        \newcommand{\mD}{\vb{D}}\newcommand{\mG}{\vb{G}}\newcommand{\mI}{\vb{I}}
        \newcommand{\mM}{\vb{M}}\newcommand{\mS}{\vb{S}}\newcommand{\mT}{\vb{T}}
        \newcommand{\vU}{\vb{U}}\newcommand{\mU}{\vb{U}}\newcommand{\mV}{\vb{V}}
        \newcommand{\vY}{\vb{Y}} \newcommand{\mZ}{\vb{Z}}
        %
        \vY(s) = \left( \mC \inv{s\mI - \mA} \mB  + \mD \right) \vU(s)
               =: \mG(s)\, \vU(s) \,.

    The entry :math:`G_{i,j}(s)` of the matrix function :math:`\mG(s): \IC \mapsto
    \IC^{q\times p}` represents the transfer function of the :math:`i`-th output and
    :math:`j`-th input (while the other inputs are zero). Hence, the frequency response
    for a given angular frequency :math:`\omega` can be expressed as
    :math:`\mG(\jj\omega)`.

    The transmission zeros (which are calculated for ``om_spec = N``) are defined as the
    values of :math:`s` for which :math:`\mG(s)` drops below its normal matrix rank.

    .. dropdown:: Implementation Notes
        :color: primary

        The complex Schur decomposition (i.e., :func:`~scipy.linalg.schur`) :math:`\mA
        = \mZ\mT\conjT{\mZ}`, with :math:`\mT` being triangular and :math:`\mZ` being
        unitary (i.e., :math:`\mZ\conjT{\mZ}= \mI`), is used to transform the transfer
        function into

        .. math::

            \mG(s) = \mC \mZ \inv{s\mI - \mT} \conjT{\mZ} \mB + \mD \,.

        For each given :math:`\omega`, the triangular matrix :math:`\mM_\omega :=
        \jj\omega \mI - \mT` is computed. This allows using the efficient
        :func:`~scipy.linalg.solve_triangular` function to calculate
        :math:`\mM_\omega^{-1} \conjT{\mZ} \mB` if :math:`\mM_\omega` is invertible.
        Otherwise, a singular value decomposition (i.e., :func:`~scipy.linalg.svd`) is
        used. Inverting the zero-valued singular values produces infinite values in
        :math:`\mG(\jj\omega)`. Remarks:

        * Since :math:`\mM_\omega` is triangular, the eigenvalues of :math:`\mM_\omega`
          are the diagonal values. Hence, if all diagonal values are non-zero,
          :math:`\mM_\omega` is invertible.
        * Inverting the SVD :math:`\mM_\omega = \mU \mS \conjT{V}`, i.e.
          :math:`\mM_\omega^{-1} = \mV \mS^{-1} \conjT{U}`, is achieved by splitting
          :math:`\mS^{-1}` into a finite and an infinite part to separately calculate a
          finite and infinite part of the transfer function.
        * There are at most :math:`n` distinct values of :math:`\jj\omega`, where the
          singular value decomposition needs to be used, since :math:`\mM_\omega` can
          have no more than :math:`n` zero-valued eigenvalues.

    Examples
    --------
    Consider the transfer functions of damped oscillators characterized by their
    undamped angular frequency :math:`\omega_i` and their damping ratio
    :math:`d_i \in [0, 1)`, i.e.,

    .. math::

        H_i(s)  = \frac{ 1 }{\frac{1}{\omega_i^2} s^2 + 2\frac{d_i}{\omega_i} s + 1}
                = \frac{ \omega_i^2 }{ (s-p^{(i)}_0) (s-p^{(i)}_1) } \,, \quad
        p^{(i)}_{0,1} = \omega_i \left( -d_i \pm \sqrt{d_i^2 - 1} \right)  \,.

    For :math:`d_i \leq 1/\sqrt{2}`, the magnitude response :math:`|H_i(\jj\omega)|` has
    a maximum

    .. math::

        |H_i(\jj\omega^{(i)}_\max)| = \frac{1}{2d_i \sqrt{1 - d_i^2}}
        \quad\text{at}\quad
        \omega^{(i)}_\max = \pm\omega_i \sqrt{1 - 2 d_i^2} \,.

    For an undamped system, i.e., :math:`d_i=0`, the maximum magnitude is unbounded,
    i.e., :math:`|H_i(\jj\omega^{(i)}_\max)| = \infty`. For this example two oscillators
    linked in series by the coupling factor :math:`\gamma = 1.2` are investigated. The
    state space system is given by

    .. math::

        \begin{bmatrix} \dot{x}_0(t)\\ \dot{x}_1(t)\\ \dot{x}_2(t)\\ \dot{x}_3(t)
        \end{bmatrix} &=
        \underbrace{\begin{bmatrix}
                          0 &             1 &            0  &             0 \\
              -\omega_0^{2} & -2d_0\omega_0 &            0  &             0 \\
                          0 &             0 &            0  &             1 \\
         \gamma\omega_1^{2} &             0 & -\omega_1^{2} & -2d_1\omega_1
        \end{bmatrix}}_{=:\mA}
        \begin{bmatrix} x_0(t)\\ x_1(t)\\ x_2(t)\\ x_3(t) \end{bmatrix} +
        \underbrace{\begin{bmatrix} 0 & 0 \\ \omega_0^2 & 0\\ 0 & 0 \\ 0 & \omega_1^2
        \end{bmatrix}}_{=:\mB}
        \begin{bmatrix} u_0(t)\\ u_1(t)\\  \end{bmatrix} ,\\
        \begin{bmatrix} y_0(t)\\ y_1(t)\\  \end{bmatrix} &=
        \underbrace{\begin{bmatrix}  1 & 0 & 0 & 0\\  0 & 0 & 1 & 0
        \end{bmatrix}}_{=:\mC}
        \begin{bmatrix} x_0(t)\\ x_1(t)\\ x_2(t)\\ x_3(t) \end{bmatrix} ,

    with the undamped angular frequencies :math:`\omega_0 = 5\,`\ rad/s,
    :math:`\omega_1 = 15\,`\ rad/s and the damping ratios :math:`d_0 = 0.1\,,`
    :math:`d_1 = 0\,.` The following plot depicts the state space's frequency response,
    which can be formulated as

    .. math::

        \mG(\jj\omega) =
        \begin{bmatrix} G_{0,0}(\jj\omega) & G_{0,1}(\jj\omega) \\
                        G_{1,0}(\jj\omega) & G_{1,1}(\jj\omega) \end{bmatrix} =
        \begin{bmatrix} H_0(\jj\omega) & 0 \\
               H_1(\jj\omega)\, \gamma\, H_0(\jj\omega) & H_1(\jj\omega)\end{bmatrix}\,.

    I.e., the transfer function of the first/second input and first/second output is
    characterized by the first/second oscillator. The transfer function of the first
    input and second output resembles the linked oscillators.

    .. plot::

        from itertools import product
        import matplotlib.pyplot as plt
        import numpy as np
        from scipy.signal import freqs_ss

        om0, om1, d0, d1, gamma  = 5, 15, 0.1, 0., 1.2
        AA = np.array([[           0,         1,       0,         0],
                       [     -om0**2, -2*d0*om0,       0,         0],
                       [           0,         0,       0,         1],
                       [gamma*om1**2,         0, -om1**2, -2*d1*om1]])
        BB = np.array([[0, om0**2, 0, 0], [0, 0, 0, om1**2]]).T
        CC, DD = np.array([[1, 0, 0, 0], [0, 0, 1, 0]]), np.zeros((2, 2))

        om, GG = freqs_ss(AA, BB, CC, DD, om_spec=(0, 20, 200), f_scale='linear')

        om0_max = om0 * np.sqrt(1 - 2*d0**2)
        mag_H0_max = 1 / (2*d0 * np.sqrt(1 - d0**2))  # |H_0(j 𝜔^0_max)|

        fig, ax = plt.subplots(constrained_layout=True)
        ax.set(title="Frequency Response", xlim=(om[0], om[-1]), ylim=(-0.1, 7),
               xlabel=r"Angular frequency $\omega$ in rad/s", xticks=(0, 5, 10, 15, 20),
               ylabel=r"Magnitude $|G_{i,j}(j\omega)|$")
        for i, j in product(range(2), range(2)):
            ax.plot(om, abs(GG[i, j]), f'C{i*2+j}-', alpha=0.7,
                    label=rf"$|G_{{{i},{j}}}(j\omega)|$")
        ax.plot(om0_max, mag_H0_max, 'C4x', alpha=.7,
                label=r'$|H_0(j\omega^{(0)}_\max)|$')
        ax.grid(True)
        ax.legend()
        plt.show()

    The plot shows the finite resonance of the damped oscillator at
    :math:`\omega^{(0)}_\max = 4.950\,`\ rad/s as well as the infinite resonance of the
    undamped oscillator at :math:`\omega^{(1)}_\max = 15\,`\ rad/s. The linked
    oscillators exhibit both resonances. Furthermore, the maximum
    :math:`|H_0(\jj\omega^{(0)}_\max)| = 5.025` is marked to illustrate that it
    coincides with the closed form expression given above.

    Note that `freqs_ss` may return infinite responses if :math:`\omega` is an
    eigenvalue of :math:`\mA`. The following trivial system with the eigenvalues
    :math:`\lambda_{0,1} = \jj/4, \jj/2` demonstrates this:

    >>> import numpy as np
    >>> from scipy.signal import freqs_ss
    ...
    >>> ABCD = np.diag([0.25j, 0.5j]), np.eye(2), np.diag([1, 0]), np.zeros((0, 0))
    >>> om = np.array([0, .25, 0.5])  # angular frequencies
    >>> _, G = freqs_ss(*ABCD, om_spec=om)  # the response
    ...
    >>> print(G[..., 0])  # G(j*om[0]) is finite
    [[0.+4.j 0.+0.j]
     [0.+0.j 0.+0.j]]
    >>> print(G[..., 1])  # infinite entry, since eigenvalue 1j*om[1] passes through C.
    [[inf+infj  0. +0.j]
     [ 0. +0.j  0. +0.j]]
    >>> print(G[..., 2])  # finite, since eigenvalue 1j*om[2] is blocked by C.
    [[ 0.-4.j  0.+0.j]
     [ 0.+0.j  0.+0.j]]
    """
    if f_scale not in ('linear', 'log'):
        raise ValueError(f"Parameter {f_scale=} not in ('linear', 'log')!")

    A, B, C, D = abcd_normalize(A, B, C, D)
    xp = array_namespace(A, B, C, D)
    device = xp_result_device(A, B, C, D)

    dtyp = xp.result_type(A, B, C, D, xp.float64)  # ensure float64 / complex128 dtype
    A, B, C, D = (xp.astype(X_, dtyp, copy=False) for X_ in (A, B, C, D))

    T, Z = _schur_complex(A)  # A = Z @ T @ Z^H; T, Z are complex-valued (n, n) arrays

    match om_spec:
        case N if isinstance(N, numbers.Integral) and N > 0:
            pp = xp.linalg.diagonal(T)  # system poles are the eigenvalues of A
            zz = _transmission_zeros(A, B, C, D, inf_zeros=False, rtol=rtol, atol=atol)
            om = findfreqs(zz, pp, N=N, kind='zp', scale=f_scale)
        case (om0, om1, N) if isinstance(om_spec, tuple):
            if not (isinstance(N, numbers.Integral) and N > 0 and
                    isinstance(om0, numbers.Real) and isinstance(om1, numbers.Real)
                    and om0 < om1):
                raise ValueError(f"Parameter {om_spec=} is not a triple (om0, om1, N) "
                                 f"with om0 < om1 and N being a positive integer.")
            if f_scale == 'linear':
                om = xp.linspace(om0, om1, N, dtype=xp.float64, device=device)
            else:  # f_scale == 'log'
                if om0 <= 0:
                    raise ValueError(f"Parameter {om_spec=} may contain only " +
                                     f"positive frequencies for parameter {f_scale=}!")
                om = _geomspace(om0, om1, N, xp=xp, device=device)
        case _:
            om = xpx.atleast_nd(xp.asarray(om_spec), ndim=1, xp=xp)
            if not (om.ndim == 1 and xp_size(om) and
                    xp.isdtype(om.dtype, ('real floating', 'integral'))):
                raise ValueError("Parameter `om_spec` must either be a positive int, " +
                                 "a valid tuple (float, float, int), " +
                                 "or a real-valued 1d array with at least one value.")
            om = xp.astype(om, xp.float64, copy=False)
    if xp.all(B == 0) or xp.all(C == 0):  # constant transfer function
        N = om.shape[0]
        return om, xp.moveaxis(xp.tile(D, (N, 1, 1)), 0, -1)  # return (q, p, N) array

    # Calculate G = C Z (j𝜔 I - T)^{-1} Z_h B + D:
    Z, C, B = xp_promote(Z, C, B, xp=xp) # PyTorch doesn't promote real → complex in `@`
    C_Z, B_Zh = C @ Z, conjT(Z) @ B
    G = _calcSSTF_helper(T, B_Zh, C_Z, 1j*om, rtol, atol, True, xp, device) + D
    return om, xp.moveaxis(G, 0, -1)
