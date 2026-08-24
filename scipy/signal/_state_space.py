"""Functions for working with state-space systems. """
import numbers

import numpy as np

import scipy._external.array_api_extra as xpx
from scipy._lib._array_api import array_namespace, xp_capabilities, xp_result_device
from scipy.linalg import schur, solve_triangular, svd

from ._lti_conversion import abcd_normalize

__all__ = ['freqs_ss']


def _geomspace(start, stop, num, xp):
    """Portable version of `numpy.geomspace` that works with all array API backends. """
    a = np.geomspace(start, stop, num=num)
    return xp.asarray(a)


@xp_capabilities(skip_backends=[(b_, "Not tested successfully") for b_ in
                                ("dask.array", "jax.numpy", "torch")])
def _calcSSTF_helper(T, B, C, s_values, rtol, atol, overwrite_b, xp, device):
    """Helper for calculating G = C (s I - T)^{-1} B assuming T is upper triangular.

    If `T` is invertible, `solve_triangular` is used for inversion. Else, a singular
    value decomposition is used. Care is taken that the zero-valued singular values are
    propagated correctly to the ∞-valued entries of G.

    Currently, this helper is only used by `freqs_ss`. This function was factored out
    to allow easier unit testing.

    Parameters
    ----------
    T, B, C : array
        `T`, `B`, `C` are complex or real-valued arrays of
        shape (n, n), (n, p), (q, n). `T` must be an upper triangular matrix.
    s_values : array
        Array of shape (N,) providing the values for s.
    rtol, atol : float
        Tolerance for considering the eigenvalues and singular values of `T` being zero.
        I.e., ``s = 0 if abs(s) < rtol * abs(lam_max) + atol else s`` with ``lam_max`
        being the eigenvalue of `T` with the largest absolute value.
    overwrite_b : bool
        Toggle whether the parameter `B` is potentially overwritten while calculating
        the result.

    Returns
    -------
    G : array
        An array of the shape (N, q, p) of the N matrices with the
        results of ``G[i] = C @ inv(s_values[i] * eye(n) - T) @ B``.

    Notes
    -----
    Adding / multiplying one finite and one infinite value always works, resulting
    in ±∞, while ∞ - ∞ = NaN and 0 * ∞ = NaN. Hence, matrix multiplication needs special
    treatment.

    Working with (and returning) (N, q, p) arrays instead of (q, p, N) arrays,
    simplifies indexing in this function considerably (e.g., `xp.diagonal` always works
    on the last axes).
    """
    n, q, p, N = T.shape[0], C.shape[0], B.shape[1], s_values.shape[0]

    # M = sI - T of shape (N, n, n):
    M = xp.reshape(s_values, (N, 1, 1)) * xp.eye(n, device=device) - T
    G = xp.zeros((N, q, p), dtype=M.dtype, device=device)  # return value

    abs_lamsM = xp.abs(xp.linalg.diagonal(M))  # absolute values of eigenvalues
    tol = rtol * xp.max(abs_lamsM) + atol
    ii = xp.all(abs_lamsM > tol, axis=-1)  # invertible matrices
    if xp.any(ii):  # Use triangular solver for invertible matrices (no Array API)
        M_np, B_np, C_np = np.asarray(M[ii]), np.asarray(B), np.asarray(C)
        G_ii = C_np @ solve_triangular(M_np, B_np, lower=False, overwrite_b=overwrite_b)
        G[ii] = xp.asarray(G_ii)
    if not xp.all(ii):  # Use SVD for non-invertible matrices (no Array API):
        U, s, Vh = svd(np.asarray(M[~ii]), overwrite_a=True)  # M = U @ S @ Vh
        U, s, Vh = xp.asarray(U), xp.asarray(s), xp.asarray(Vh)
        V, Uh_B = xp.conj(Vh).mT, xp.conj(U).mT @ B

        # Split S^{-1} into a finite and an infinite part:
        inv_s, ii0_s = xp.zeros_like(s), s > tol
        inv_s[ii0_s] = 1 / s[ii0_s]  # finite
        V_iS = V * xp.reshape(inv_s, (-1, 1, n))  # V_iS = V @ S^{-1}
        G_finite = C @ V_iS @ Uh_B

        # Infinite part: Mark all ∞ values in S^{-1} as 1 and all others as 0.
        # Then propagate those marks to G:
        ii1_s = inv_s == 0
        inv_s[ii1_s], inv_s[~ii1_s] = 1, 0
        V_iS = V * xp.reshape(inv_s, (-1, 1, n))  # V_iS = V @ S^{-1}
        G_inf = C @ V_iS @ Uh_B
        G_inf[G_inf != 0] = xp.inf  # non-zero values become ∞

        G[~ii] += G_finite + G_inf
    return G


@xp_capabilities()
def _parse_om_spec(om_spec, f_scale, om0, om1, xp, device):
    """Helper for generating frequencies.

    Currently, this function is only used by `freq_ss`. Consult it for a detailed
    parameter description. Factoring it out allows easier unit testing.

    Parameters
    ----------
    om_spec : int | (float, float, int) | array_like
        Specification of the frequencies. If an integer, it specifies the number of
        frequencies to use. If a tuple, it specifies the start and end frequencies and
        the number of frequencies to use. If an array, it specifies the frequencies
        directly.
    f_scale : 'linear' | 'log'
        The frequency scale can be either linear or logarithmic. This parameter is
        ignored if `om_spec` is an array.
    om0, om1 : float
        The minimum and maximum frequencies to use if `om_spec` is an integer.
    """
    if f_scale not in ('linear', 'log'):
        raise ValueError(f"Parameter {f_scale=} not in ('linear', 'log')!")
    match om_spec:
        case N if isinstance(N, numbers.Integral) and N > 0:
            if f_scale == 'linear':
                return xp.linspace(0, 1.5 * om1, N, device=device)
            return _geomspace(om0 / 5, 5 * om1, N, xp=xp)
        case (om0, om1, N) if isinstance(om_spec, tuple):
            if not (isinstance(N, numbers.Integral) and N > 0 and
                    isinstance(om0, numbers.Real) and isinstance(om1, numbers.Real)
                    and om0 < om1):
                raise ValueError(f"Parameter {om_spec=} is not a triple (om0, om1, N) "
                                 f"with om0 < om1 and N being a positive integer.")
            if f_scale == 'linear':
                return xp.linspace(om0, om1, N, device=device)
            return _geomspace(om0, om1, N, xp=xp)
        case _:
            err_msg = ("Parameter `om_spec` must either be a positive int, a valid " +
                       "tuple (float, float, int), or a  real-valued 1d array with " +
                       "at least one value.")
            try:
                om = xpx.atleast_nd(xp.asarray(om_spec), ndim=1, xp=xp)
            except TypeError as e:
                raise ValueError(err_msg) from e
            if not (om.ndim == 1 and om.shape[0] > 0 and
                    xp.isdtype(om.dtype, ('real floating', 'integral'))):
                raise ValueError(err_msg)
            return om


def freqs_ss(A, B, C, D, *, om_spec, f_scale='linear', rtol=1e-12, atol=0.):
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
        The frequency range specification `om_spec` may have the three different forms:

        If `om_spec` is an array, it is used directly as the angular frequency values.
        The array must be a real-valued one-dimensional array of at least length one.
        The parameter `f_scale` is ignored in this case.

        If `om_spec` is a tuple of the form ``(w0, w1, N)``, i.e., of type
        ``tuple[float, float, int]``, the angular frequency ranges from ``w0`` to ``w1``
        with ``N`` steps. The parameter `f_scale` toggles between "linear" and "log"
        scaling.

        If `om_spec` is a positive integer ``N``, ``w0`` and ``w1`` are heuristically
        selected: For linear scaling ``w0, w1 = 0, 1.5*ilam1`` and for logarithmic
        scaling ``om0, om1 = ilam0/5, ilam1*5``. Here, ``ilam0``, ``ilam1`` are the
        minimum and the maximum of the nonzero absolute values of `A`\ 's eigenvalues'
        imaginary parts. Missing values of ``ilam0``, ``ilam1`` are set to ``ilam1/2``
        and ``1`` as needed.
    f_scale : 'linear' | 'log', optional
        The frequency scale can be either linear (default) or logarithmic.
        This parameter is ignored if `om_spec` is an array.
    rtol, atol : float, optional
        Eigenvalues or singular values of `A` whose absolute values are less than ``rtol
        * abs(lam_max) + atol`` are considered to be zero. Here, ``lam_max`` is the
        eigenvalue of `A` with the largest absolute value.
        The default is ``rtol, atol = 1e-12, 0``.

    Returns
    -------
    om : array
        A real-valued array of angular frequencies of shape (N,) for which the response
        is calculated.
    G : array
        A complex-valued array of shape (q, p, N) representing the frequency response.

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
    `D` are finite.

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
                = \frac{ \omega_i^2 }{ (s-p^{(i)}_0) (s-p^{(i)}_1) } \,, \qquad
        p^{(i)}_{0,1} = \omega_i \left( -d_i ± \sqrt{d_i^2 - 1} \right)  \,.

    The maximum magnitude can be expressed as

    .. math::

        |H_i(\jj\omega^{(i)}_\max)| = \frac{1}{2d_i \sqrt{1 - d_i^2}}
        \quad\text{with}\quad
        \omega^{(i)}_\max = \omega_i \sqrt{1 - 2 d_i^2} \,.

    For an undamped system, i.e., :math:`d_i=0`, the maximum magnitude is unbounded,
    i.e., :math:`|H_i(\jj\omega_\max)| = \infty`. As an example, two oscillators,
    linked in series by the coupling factor :math:`\gamma = 1`, are investigated. The
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
        \begin{bmatrix} u_0(t)\\ u_1(t)\\  \end{bmatrix} \,, \\
        \begin{bmatrix} y_0(t)\\ y_1(t)\\  \end{bmatrix} &=
        \underbrace{\begin{bmatrix}  1 & 0 & 0 & 0\\  0 & 0 & 1 & 0
        \end{bmatrix}}_{=:\mC}
        \begin{bmatrix} x_0(t)\\ x_1(t)\\ x_2(t)\\ x_3(t) \end{bmatrix} \,,

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

        om0, om1, d0, d1, gamma  = 5, 15, 0.1, 0., 1.
        AA = np.array([[           0,         1,       0,         0],
                       [     -om0**2, -2*d0*om0,       0,         0],
                       [           0,         0,       0,         1],
                       [gamma*om1**2,         0, -om1**2, -2*d1*om1]])
        BB = np.array([[0, om0**2, 0, 0], [0, 0, 0, om1**2]]).T
        CC, DD = np.array([[1, 0, 0, 0], [0, 0, 1, 0]]), np.zeros((2, 2))

        om, GG = freqs_ss(AA, BB, CC, DD, om_spec=(0, 20, 200))  # the response

        om0_max = om0 * np.sqrt(1 - 2*d0**2)
        mH0_max = 1 / (2*d0 * np.sqrt(1 - d0**2))  # |H_0(j 𝜔^0_max)|

        fig, ax = plt.subplots(constrained_layout=True)
        ax.set(title="Frequency Response", xlim=(om[0], om[-1]), ylim=(-0.1, 7),
               xlabel=r"Angular frequency $\omega$ in rad/s", xticks=(0, 5, 10, 15, 20),
               ylabel=r"Magnitude $|G_{i,j}(j\omega)|$")
        for i, j in product(range(2), range(2)):
            ax.plot(om, abs(GG[i, j]), f'C{i*2+j}-', alpha=0.7,
                    label=rf"$|G_{{{i},{j}}}(j\omega)|$")
        ax.plot(om0_max, mH0_max, 'C4x', alpha=.7,
                label=r'$|H_0(j\omega^{(0)}_\max)|$')
        ax.grid(True)
        ax.legend()
        plt.show()

    The plot shows the finite resonance of the damped oscillator at
    :math:`\omega^{(0)}_\max = 4.950\,`\ rad/s as well as the infinite resonance of the
    undamped oscillator at :math:`\omega^{(1)}_\max = 15\,`\ rad/s. The product of both
    oscillators exhibits both resonances. Furthermore, the maximum
    :math:`|H_0(\jj\omega^{(0)}_\max)| = 5.025` is marked to illustrate that it
    coincides with the closed form expression given above.

    Note that `freq_ss` may return infinite responses if :math:`\omega` is an
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
    [[inf+0.j  0.+0.j]
     [ 0.+0.j  0.+0.j]]
    >>> print(G[..., 2])  # finite, since eigenvalue 1j*om[2] is blocked by C.
    [[ 0.-4.j  0.+0.j]
     [ 0.+0.j  0.+0.j]]
    """
    A, B, C, D = abcd_normalize(A, B, C, D)
    xp = array_namespace(A, B, C, D)
    device = xp_result_device(A, B, C, D)


    A_np = np.asarray(A)  # No Array API support yet for schur()
    T, Z = schur(A_np, output='complex')  # A = Z @ T @ Z^h — T, Z have shape (n, n)
    T, Z = xp.asarray(T), xp.asarray(Z)

    # The maximum and minimum of the imaginary parts' absolute values of A's non-zero
    # eigenvalues are used for determining automatic frequency range:
    im_abs_lam = xp.abs(xp.imag(xp.linalg.diagonal(T)))
    om_max = xp.max(xp.abs(im_abs_lam)) if im_abs_lam.shape[0] > 0 else 0.
    im_abs_lam = im_abs_lam[im_abs_lam >= rtol*om_max + atol]  # only the positive parts
    om_min = xp.min(xp.abs(im_abs_lam)) if im_abs_lam.shape[0] > 0 else 0.

    if om_min == om_max == 0.:  # all eigenvalues are real
        om_min, om_max = 1., 1.
    elif om_min == 0.:  # some eigenvalues are real
        om_min = om_max / 2.
    om = _parse_om_spec(om_spec, f_scale, om_min, om_max, xp, device)  # frequencies

    if xp.all(B == 0) or xp.all(C == 0):  # constant transfer function
        N = om.shape[0]
        return om, xp.moveaxis(xp.tile(D, (N, 1, 1)), 0, -1)  # return (N, q, p) array

    # Calculate G = C Z (j𝜔 I - T)^{-1} Z_h B + D:
    C_Z, B_Zh = C @ Z, xp.conj(Z).mT @ B
    G = _calcSSTF_helper(T, B_Zh, C_Z, 1j * om, rtol, atol, True, xp, device) + D
    return om, xp.moveaxis(G, 0, -1)
