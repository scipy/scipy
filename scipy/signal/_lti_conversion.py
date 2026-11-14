"""
ltisys -- a collection of functions to convert linear time invariant systems
from one representation to another.
"""

import numpy as np
from scipy import linalg

from scipy._lib._array_api import (array_namespace, xp_size, xp_promote,
                                   xp_result_type, xp_device)
from scipy._external.array_api_compat import numpy as np_compat
import scipy._external.array_api_extra as xpx
from scipy.signal import _polyutils as _pu
from ._filter_design import tf2zpk, zpk2tf, normalize


__all__ = ['tf2ss', 'abcd_normalize', 'ss2tf', 'zpk2ss', 'ss2zpk',
           'cont2discrete']


def _char_poly(M, *, xp):
    """Characteristic polynomial of a square matrix, like 2-D `np.poly`."""
    return _pu.poly(_pu.eigvals(M, xp=xp), xp=xp)


def tf2ss(num, den):
    r"""Transfer function to state-space representation.

    Parameters
    ----------
    num, den : array_like
        Sequences representing the coefficients of the numerator and
        denominator polynomials, in order of descending degree. The
        denominator needs to be at least as long as the numerator.

    Returns
    -------
    A, B, C, D : ndarray
        State space representation of the system, in controller canonical
        form.

    Examples
    --------
    Convert the transfer function:

    .. math:: H(s) = \frac{s^2 + 3s + 3}{s^2 + 2s + 1}

    >>> num = [1, 3, 3]
    >>> den = [1, 2, 1]

    to the state-space representation:

    .. math::

        \dot{\textbf{x}}(t) =
        \begin{bmatrix} -2 & -1 \\ 1 & 0 \end{bmatrix} \textbf{x}(t) +
        \begin{bmatrix} 1 \\ 0 \end{bmatrix} \textbf{u}(t) \\

        \textbf{y}(t) = \begin{bmatrix} 1 & 2 \end{bmatrix} \textbf{x}(t) +
        \begin{bmatrix} 1 \end{bmatrix} \textbf{u}(t)

    >>> from scipy.signal import tf2ss
    >>> A, B, C, D = tf2ss(num, den)
    >>> A
    array([[-2., -1.],
           [ 1.,  0.]])
    >>> B
    array([[ 1.],
           [ 0.]])
    >>> C
    array([[ 1.,  2.]])
    >>> D
    array([[ 1.]])
    """
    # Controller canonical state-space representation.
    #  if M+1 = len(num) and K+1 = len(den) then we must have M <= K
    #  states are found by asserting that X(s) = U(s) / D(s)
    #  then Y(s) = N(s) * X(s)
    #
    #   A, B, C, and D follow quite naturally.
    #
    num, den = normalize(num, den)   # Strips zeros, checks arrays
    try:
        xp = array_namespace(num, den)
    except TypeError:
        # object arrays, test_ltisys.py::TestSS2TF::test_simo_round_trip
        xp = np_compat
    device = xp_device(den)

    if num.ndim == 1:
        num = num[None, :]
    M = num.shape[1]
    K = den.shape[0]
    if M > K:
        msg = "Improper transfer function. `num` is longer than `den`."
        raise ValueError(msg)
    if M == 0 or K == 0:  # Null system
        return tuple(xp.asarray([], dtype=xp.float64, device=device)
                     for _ in range(4))

    # pad numerator to have same number of columns has denominator
    num = xp.concat(
        (xp.zeros((num.shape[0], K - M), dtype=num.dtype, device=device), num),
        axis=1)

    if num.shape[-1] > 0:
        D = xpx.atleast_nd(num[:, 0], ndim=2, xp=xp)
    else:
        # We don't assign it an empty array because this system
        # is not 'null'. It just doesn't have a non-zero D
        # matrix. Thus, it should have a non-zero shape so that
        # it can be operated on by functions like 'ss2tf'
        D = xp.asarray([[0.0]], device=device)

    if K == 1:
        D = xp.reshape(D, num.shape)

        return (xp.zeros((1, 1), dtype=D.dtype, device=device),
                xp.zeros((1, D.shape[1]), dtype=D.dtype, device=device),
                xp.zeros((D.shape[0], 1), dtype=D.dtype, device=device),
                D)

    frow = -den[1:][None, :]
    # a sliced square eye, dask chunks eye(K - 2, K - 1) per column
    eye = xp.eye(K - 1, dtype=frow.dtype, device=device)
    A = xp.concat((frow, eye[:K - 2, :]), axis=0)
    B = eye[:, :1]
    C = num[:, 1:] - num[:, 0][:, None] * den[1:][None, :]
    D = xp.reshape(D, (C.shape[0], B.shape[1]))

    return A, B, C, D


def abcd_normalize(A=None, B=None, C=None, D=None):
    r"""Check state-space matrices compatibility and ensure they are 2d arrays.

    First, the input matrices are converted into two-dimensional arrays with
    appropriate dtype as needed. Then the dimensions n, q, p are determined by
    investigating the array shapes. If an input is ``None``, or has size zero, it is
    set to an array of zeros of compatible shape. Finally, it is verified that all
    parameter shapes are compatible with each other. If that fails, a ``ValueError`` is
    raised. Note that the dimensions n, q, p are allowed to be zero.

    Parameters
    ----------
    A : array_like, optional
        Two-dimensional array of shape (n, n).
    B : array_like, optional
        Two-dimensional array of shape (n, p).
    C : array_like, optional
        Two-dimensional array of shape (q, n).
    D : array_like, optional
        Two-dimensional array of shape (q, p).

    Returns
    -------
    A, B, C, D : array
        State-space matrices as two-dimensional arrays with identical dtype.
        The result dtype is determined based on the standard
        `dtype promotion rules <https://numpy.org/doc/2.3/reference/arrays.promotion.html>`_
        except for when the inputs are all of integer dtype, in which case the returned
        arrays will have the default floating point dtype of ``float64``.

    Raises
    ------
    ValueError
        If the dimensions n, q, or p could not be determined or if the shapes are
        incompatible with each other.

    Notes
    -----
    If a matrix is not modified, the original matrix (not a copy) is returned.

    The :ref:`tutorial_signal_state_space_representation` section of the
    :ref:`user_guide` presents the corresponding definitions of continuous-time and
    disrcete time state space systems.

    See Also
    --------
    StateSpace: Linear Time Invariant system in state-space form.
    dlti: Discrete-time linear time invariant system base class.
    tf2ss: Transfer function to state-space representation.
    ss2tf: State-space to transfer function.
    ss2zpk: State-space representation to zero-pole-gain representation.
    cont2discrete: Transform a continuous to a discrete state-space system.

    Examples
    --------
    The following example demonstrates that the passed lists are converted into
    two-dimensional arrays:

    >>> from scipy.signal import abcd_normalize
    >>> AA, BB, CC, DD = abcd_normalize(A=[[1, 2], [3, 4]], B=[[-1], [5]],
    ...                                 C=[[4, 5]], D=2.5)
    >>> AA.shape, BB.shape, CC.shape, DD.shape
    ((2, 2), (2, 1), (1, 2), (1, 1))

    In the following, the missing parameter C is assumed to be an array of zeros
    with shape (1, 2):

    >>> from scipy.signal import abcd_normalize
    >>> AA, BB, CC, DD = abcd_normalize(A=[[1, 2], [3, 4]], B=[[-1], [5]], D=2.5)
    >>> AA.shape, BB.shape, CC.shape, DD.shape
    ((2, 2), (2, 1), (1, 2), (1, 1))
    >>> CC
    array([[0., 0.]])

    """
    if A is None and B is None and C is None:
        raise ValueError("Dimension n is undefined for parameters A = B = C = None!")
    if B is None and D is None:
        raise ValueError("Dimension p is undefined for parameters B = D = None!")
    if C is None and D is None:
        raise ValueError("Dimension q is undefined for parameters C = D = None!")

    xp = array_namespace(A, B, C, D)
    # after `xp_promote` all non-None inputs share one device; at least one
    # of A, B, C is not None per the checks above
    A, B, C, D = xp_promote(A, B, C, D, xp=xp, force_floating=True)
    dtype = xp_result_type(A, B, C, D, xp=xp)
    device = xp_device(next(M_ for M_ in (A, B, C, D) if M_ is not None))

    # convert inputs into 2d arrays (zero-size 2d array if None):
    A, B, C, D = (
        xpx.atleast_nd(xp.asarray(M_), ndim=2, xp=xp)
        if M_ is not None else xp.zeros((0, 0), dtype=dtype, device=device)
        for M_ in (A, B, C, D)
    )

    n = A.shape[0] or B.shape[0] or C.shape[1] # try finding non-zero dimensions
    p = B.shape[1] or D.shape[1]
    q = C.shape[0] or D.shape[0]

    # Create zero matrices as needed:
    A = xp.zeros((n, n), dtype=dtype, device=device) if xp_size(A) == 0 else A
    B = xp.zeros((n, p), dtype=dtype, device=device) if xp_size(B) == 0 else B
    C = xp.zeros((q, n), dtype=dtype, device=device) if xp_size(C) == 0 else C
    D = xp.zeros((q, p), dtype=dtype, device=device) if xp_size(D) == 0 else D

    if A.shape != (n, n):
        raise ValueError(f"Parameter A has shape {A.shape} but should be ({n}, {n})!")
    if B.shape != (n, p):
        raise ValueError(f"Parameter B has shape {B.shape} but should be ({n}, {p})!")
    if C.shape != (q, n):
        raise ValueError(f"Parameter C has shape {C.shape} but should be ({q}, {n})!")
    if D.shape != (q, p):
        raise ValueError(f"Parameter D has shape {D.shape} but should be ({q}, {p})!")

    return A, B, C, D


def ss2tf(A, B, C, D, input=0):
    r"""State-space to transfer function.

    A, B, C, D defines a linear state-space system with `p` inputs,
    `q` outputs, and `n` state variables.

    Parameters
    ----------
    A : array_like
        State (or system) matrix of shape ``(n, n)``
    B : array_like
        Input matrix of shape ``(n, p)``
    C : array_like
        Output matrix of shape ``(q, n)``
    D : array_like
        Feedthrough (or feedforward) matrix of shape ``(q, p)``
    input : int, optional
        For multiple-input systems, the index of the input to use.

    Returns
    -------
    num : 2-D ndarray
        Numerator(s) of the resulting transfer function(s). `num` has one row
        for each of the system's outputs. Each row is a sequence representation
        of the numerator polynomial.
    den : 1-D ndarray
        Denominator of the resulting transfer function(s). `den` is a sequence
        representation of the denominator polynomial.

    Notes
    -----
    Before calculating `num` and `den`, the function `abcd_normalize` is called to
    convert the parameters `A`, `B`, `C`, `D` into two-dimensional arrays of the
    same dtype. The resulting dtype will be based on NumPy's dtype promotion rules,
    except in the case where each of `A`, `B`, `C`, and `D` has integer dtype, in which
    case the resulting dtype will be the default floating point dtype of ``float64``.

    The :ref:`tutorial_signal_state_space_representation` section of the
    :ref:`user_guide` presents the corresponding definitions of continuous-time and
    disrcete time state space systems.

    Examples
    --------
    Convert the state-space representation:

    .. math::

        \dot{\textbf{x}}(t) =
        \begin{bmatrix} -2 & -1 \\ 1 & 0 \end{bmatrix} \textbf{x}(t) +
        \begin{bmatrix} 1 \\ 0 \end{bmatrix} \textbf{u}(t) \\

        \textbf{y}(t) = \begin{bmatrix} 1 & 2 \end{bmatrix} \textbf{x}(t) +
        \begin{bmatrix} 1 \end{bmatrix} \textbf{u}(t)

    >>> A = [[-2, -1], [1, 0]]
    >>> B = [[1], [0]]  # 2-D column vector
    >>> C = [[1, 2]]    # 2-D row vector
    >>> D = 1

    to the transfer function:

    .. math:: H(s) = \frac{s^2 + 3s + 3}{s^2 + 2s + 1}

    >>> from scipy.signal import ss2tf
    >>> ss2tf(A, B, C, D)
    (array([[1., 3., 3.]]), array([ 1.,  2.,  1.]))
    """
    # transfer function is C (sI - A)**(-1) B + D

    # Check consistency and make them all floating point 2d arrays:
    A, B, C, D = abcd_normalize(A, B, C, D)

    try:
        xp = array_namespace(A, B, C, D)
    except TypeError:
        # object arrays, test_ltisys.py::TestSS2TF::test_simo_round_trip
        xp = np_compat

    nout, nin = D.shape
    if input >= nin:
        raise ValueError("System does not have the input specified.")

    # make SIMO from possibly MIMO system.
    B = B[:, input:input + 1]
    D = D[:, input:input + 1]

    try:
        den = _char_poly(A, xp=xp)
    except ValueError:
        den = 1

    if (xp_size(B) == 0) and (xp_size(C) == 0):
        num = xp.reshape(D, (-1,))
        if (xp_size(D) == 0) and (xp_size(A) == 0):
            den = []
        return num, den

    num_states = A.shape[0]
    result_dtype = xp_result_type(A, B, C, D, force_floating=True, xp=xp)
    num = xp.empty((nout, num_states + 1), dtype=result_dtype,
                   device=xp_device(A))
    for k in range(nout):
        Ck = C[k:k + 1, :]
        num_k = _char_poly(A - B @ Ck, xp=xp) + (D[k, ...] - 1.0) * den
        num = xpx.at(num)[k, ...].set(xp.astype(num_k, result_dtype))

    return num, den


def zpk2ss(z, p, k):
    """Zero-pole-gain representation to state-space representation.

    Parameters
    ----------
    z, p : sequence
        Zeros and poles.
    k : float
        System gain.

    Returns
    -------
    A, B, C, D : ndarray
        State space representation of the system, in controller canonical
        form.

    """
    return tf2ss(*zpk2tf(z, p, k))


def ss2zpk(A, B, C, D, input=0):
    """State-space representation to zero-pole-gain representation.

    A, B, C, D defines a linear state-space system with `p` inputs,
    `q` outputs, and `n` state variables.

    Parameters
    ----------
    A : array_like
        State (or system) matrix of shape ``(n, n)``
    B : array_like
        Input matrix of shape ``(n, p)``
    C : array_like
        Output matrix of shape ``(q, n)``
    D : array_like
        Feedthrough (or feedforward) matrix of shape ``(q, p)``
    input : int, optional
        For multiple-input systems, the index of the input to use.

    Returns
    -------
    z, p : sequence
        Zeros and poles.
    k : float
        System gain.

    """
    return tf2zpk(*ss2tf(A, B, C, D, input=input))


def cont2discrete(system, dt, method="zoh", alpha=None):
    """
    Transform a continuous to a discrete state-space system.

    Parameters
    ----------
    system : a tuple describing the system or an instance of `lti`
        The following gives the number of elements in the tuple and
        the interpretation:

        * 1: (instance of `lti`)
        * 2: (num, den)
        * 3: (zeros, poles, gain)
        * 4: (A, B, C, D)

    dt : float
        The discretization time step.
    method : str, optional
        Which method to use:

        * gbt: generalized bilinear transformation
        * bilinear: Tustin's approximation ("gbt" with alpha=0.5)
        * euler: Euler (or forward differencing) method ("gbt" with alpha=0)
        * backward_diff: Backwards differencing ("gbt" with alpha=1.0)
        * zoh: zero-order hold (default)
        * foh: first-order hold (*versionadded: 1.3.0*)
        * impulse: equivalent impulse response (*versionadded: 1.3.0*)

    alpha : float within [0, 1], optional
        The generalized bilinear transformation weighting parameter, which
        should only be specified with method="gbt", and is ignored otherwise

    Returns
    -------
    sysd : tuple containing the discrete system
        Based on the input type, the output will be of the form

        * (num, den, dt)   for transfer function input
        * (zeros, poles, gain, dt)   for zeros-poles-gain input
        * (A, B, C, D, dt) for state-space system input

    Notes
    -----
    By default, the routine uses a Zero-Order Hold (zoh) method to perform
    the transformation. Alternatively, a generalized bilinear transformation
    may be used, which includes the common Tustin's bilinear approximation,
    an Euler's method technique, or a backwards differencing technique.

    The Zero-Order Hold (zoh) method is based on [1]_, the generalized bilinear
    approximation is based on [2]_ and [3]_, the First-Order Hold (foh) method
    is based on [4]_.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Discretization#Discretization_of_linear_state_space_models

    .. [2] http://techteach.no/publications/discretetime_signals_systems/discrete.pdf

    .. [3] G. Zhang, X. Chen, and T. Chen, Digital redesign via the generalized
        bilinear transformation, Int. J. Control, vol. 82, no. 4, pp. 741-754,
        2009.
        (https://www.mypolyuweb.hk/~magzhang/Research/ZCC09_IJC.pdf)

    .. [4] G. F. Franklin, J. D. Powell, and M. L. Workman, Digital control
        of dynamic systems, 3rd ed. Menlo Park, Calif: Addison-Wesley,
        pp. 204-206, 1998.

    Examples
    --------
    We can transform a continuous state-space system to a discrete one:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import cont2discrete, lti, dlti, dstep

    Define a continuous state-space system.

    >>> A = np.array([[0, 1],[-10., -3]])
    >>> B = np.array([[0],[10.]])
    >>> C = np.array([[1., 0]])
    >>> D = np.array([[0.]])
    >>> l_system = lti(A, B, C, D)
    >>> t, x = l_system.step(T=np.linspace(0, 5, 100))
    >>> fig, ax = plt.subplots()
    >>> ax.plot(t, x, label='Continuous', linewidth=3)

    Transform it to a discrete state-space system using several methods.

    >>> dt = 0.1
    >>> for method in ['zoh', 'bilinear', 'euler', 'backward_diff', 'foh', 'impulse']:
    ...    d_system = cont2discrete((A, B, C, D), dt, method=method)
    ...    s, x_d = dstep(d_system)
    ...    ax.step(s, np.squeeze(x_d), label=method, where='post')
    >>> ax.axis([t[0], t[-1], x[0], 1.4])
    >>> ax.legend(loc='best')
    >>> fig.tight_layout()
    >>> plt.show()

    """
    if hasattr(system, 'to_discrete') and callable(system.to_discrete):
        return system.to_discrete(dt=dt, method=method, alpha=alpha)

    if len(system) == 2:
        sysd = cont2discrete(tf2ss(system[0], system[1]), dt, method=method,
                             alpha=alpha)
        return ss2tf(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(system) == 3:
        sysd = cont2discrete(zpk2ss(system[0], system[1], system[2]), dt,
                             method=method, alpha=alpha)
        return ss2zpk(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(system) == 4:
        a, b, c, d = system
    else:
        raise ValueError("First argument must either be a tuple of 2 (tf), "
                         "3 (zpk), or 4 (ss) arrays.")

    if method == 'gbt':
        if alpha is None:
            raise ValueError("Alpha parameter must be specified for the "
                             "generalized bilinear transform (gbt) method")
        elif alpha < 0 or alpha > 1:
            raise ValueError("Alpha parameter must be within the interval "
                             "[0,1] for the gbt method")

    if method == 'gbt':
        # This parameter is used repeatedly - compute once here
        ima = np.eye(a.shape[0]) - alpha*dt*a
        ad = linalg.solve(ima, np.eye(a.shape[0]) + (1.0-alpha)*dt*a)
        bd = linalg.solve(ima, dt*b)

        # Similarly solve for the output equation matrices
        cd = linalg.solve(ima.transpose(), c.transpose())
        cd = cd.transpose()
        dd = d + alpha*np.dot(c, bd)

    elif method == 'bilinear' or method == 'tustin':
        return cont2discrete(system, dt, method="gbt", alpha=0.5)

    elif method == 'euler' or method == 'forward_diff':
        return cont2discrete(system, dt, method="gbt", alpha=0.0)

    elif method == 'backward_diff':
        return cont2discrete(system, dt, method="gbt", alpha=1.0)

    elif method == 'zoh':
        # Build an exponential matrix
        em_upper = np.hstack((a, b))

        # Need to stack zeros under the a and b matrices
        em_lower = np.hstack((np.zeros((b.shape[1], a.shape[0])),
                              np.zeros((b.shape[1], b.shape[1]))))

        em = np.vstack((em_upper, em_lower))
        ms = linalg.expm(dt * em)

        # Dispose of the lower rows
        ms = ms[:a.shape[0], :]

        ad = ms[:, 0:a.shape[1]]
        bd = ms[:, a.shape[1]:]

        cd = c
        dd = d

    elif method == 'foh':
        # Size parameters for convenience
        n = a.shape[0]
        m = b.shape[1]

        # Build an exponential matrix similar to 'zoh' method
        em_upper = linalg.block_diag(np.block([a, b]) * dt, np.eye(m))
        em_lower = np.zeros((m, n + 2 * m))
        em = np.block([[em_upper], [em_lower]])

        ms = linalg.expm(em)

        # Get the three blocks from upper rows
        ms11 = ms[:n, 0:n]
        ms12 = ms[:n, n:n + m]
        ms13 = ms[:n, n + m:]

        ad = ms11
        bd = ms12 - ms13 + ms11 @ ms13
        cd = c
        dd = d + c @ ms13

    elif method == 'impulse':
        if not np.allclose(d, 0):
            raise ValueError("Impulse method is only applicable "
                             "to strictly proper systems")

        ad = linalg.expm(a * dt)
        bd = ad @ b * dt
        cd = c
        dd = c @ b * dt

    else:
        raise ValueError(f"Unknown transformation method '{method}'")

    return ad, bd, cd, dd, dt
