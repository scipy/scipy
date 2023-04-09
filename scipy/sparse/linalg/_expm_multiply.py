"""Compute the action of the matrix exponential."""
from warnings import warn

import numpy as np

import scipy.linalg
import scipy.sparse.linalg
from scipy.linalg._decomp_qr import qr
from scipy.sparse._sputils import is_pydata_spmatrix
from scipy.sparse.linalg import aslinearoperator
from scipy.sparse.linalg._interface import IdentityOperator
from scipy.sparse.linalg._onenormest import onenormest

__all__ = ['expm_multiply']


def _exact_inf_norm(A):
    # A compatibility function which should eventually disappear.
    if scipy.sparse.isspmatrix(A):
        return max(abs(A).sum(axis=1).flat)
    elif is_pydata_spmatrix(A):
        return max(abs(A).sum(axis=1))
    else:
        return np.linalg.norm(A, np.inf)


def _exact_1_norm(A):
    # A compatibility function which should eventually disappear.
    if scipy.sparse.isspmatrix(A):
        return max(abs(A).sum(axis=0).flat)
    elif is_pydata_spmatrix(A):
        return max(abs(A).sum(axis=0))
    else:
        return np.linalg.norm(A, 1)


def _trace(A):
    # A compatibility function which should eventually disappear.
    if is_pydata_spmatrix(A):
        return A.to_scipy_sparse().trace()
    else:
        return A.trace()


def traceest(A, m3, seed=None):
    """Estimate `np.trace(A)` using `3*m3` matrix-vector products.

    The result is not deterministic.

    Parameters
    ----------
    A : LinearOperator
        Linear operator whose trace will be estimated. Has to be square.
    m3 : int
        Number of matrix-vector products divided by 3 used to estimate the
        trace.
    seed : optional
        Seed for `numpy.random.default_rng`.
        Can be provided to obtain deterministic results.

    Returns
    -------
    trace : LinearOperator
        Estimate of the trace

    Notes
    -----
    This is the Hutch++ algorithm given in [1]_.

    References
    ----------
    .. [1] Meyer, Raphael A., Cameron Musco, Christopher Musco, and David P.
       Woodruff. "Hutch++: Optimal Stochastic Trace Estimation." In Symposium
       on Simplicity in Algorithms (SOSA), pp. 142-155. Society for Industrial
       and Applied Mathematics, 2021
       https://doi.org/10.1137/1.9781611976496.16

    """
    rng = np.random.default_rng(seed)
    if len(A.shape) != 2 or A.shape[-1] != A.shape[-2]:
        raise ValueError("Expected A to be like a square matrix.")
    n = A.shape[-1]
    S = rng.choice([-1.0, +1.0], [n, m3])
    Q, _ = qr(A.matmat(S), overwrite_a=True, mode='economic')
    trQAQ = np.trace(Q.conj().T @ A.matmat(Q))
    G = rng.choice([-1, +1], [n, m3])
    right = G - Q@(Q.conj().T @ G)
    trGAG = np.trace(right.conj().T @ A.matmat(right))
    return trQAQ + trGAG/m3


def _ident_like(A):
    # A compatibility function which should eventually disappear.
    if scipy.sparse.isspmatrix(A):
        return scipy.sparse._construct.eye(A.shape[0], A.shape[1],
                                           dtype=A.dtype, format=A.format)
    elif is_pydata_spmatrix(A):
        import sparse
        return sparse.eye(A.shape[0], A.shape[1], dtype=A.dtype)
    elif isinstance(A, scipy.sparse.linalg.LinearOperator):
        return IdentityOperator(A.shape, dtype=A.dtype)
    else:
        return np.eye(A.shape[0], A.shape[1], dtype=A.dtype)


def expm_multiply(A, B, start=None, stop=None, num=None,
                  endpoint=None, traceA=None):
    """
    Compute the action of the matrix exponential of A on B.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose exponential is of interest.
    B : ndarray
        The matrix or vector to be multiplied by the matrix exponential of A.
    start : scalar, optional
        The starting time point of the sequence.
    stop : scalar, optional
        The end time point of the sequence, unless `endpoint` is set to False.
        In that case, the sequence consists of all but the last of ``num + 1``
        evenly spaced time points, so that `stop` is excluded.
        Note that the step size changes when `endpoint` is False.
    num : int, optional
        Number of time points to use.
    endpoint : bool, optional
        If True, `stop` is the last time point.  Otherwise, it is not included.
    traceA : scalar, optional
        Trace of `A`. If not given the trace is estimated for linear operators,
        or calculated exactly for sparse matrices. It is used to precondition
        `A`, thus an approximate trace is acceptable.
        For linear operators, `traceA` should be provided to ensure performance
        as the estimation is not guaranteed to be reliable for all cases.

        .. versionadded:: 1.9.0

    Returns
    -------
    expm_A_B : ndarray
         The result of the action :math:`e^{t_k A} B`.

    Warns
    -----
    UserWarning
        If `A` is a linear operator and ``traceA=None`` (default).

    Notes
    -----
    The optional arguments defining the sequence of evenly spaced time points
    are compatible with the arguments of `numpy.linspace`.

    The output ndarray shape is somewhat complicated so I explain it here.
    The ndim of the output could be either 1, 2, or 3.
    It would be 1 if you are computing the expm action on a single vector
    at a single time point.
    It would be 2 if you are computing the expm action on a vector
    at multiple time points, or if you are computing the expm action
    on a matrix at a single time point.
    It would be 3 if you want the action on a matrix with multiple
    columns at multiple time points.
    If multiple time points are requested, expm_A_B[0] will always
    be the action of the expm at the first time point,
    regardless of whether the action is on a vector or a matrix.

    References
    ----------
    .. [1] Awad H. Al-Mohy and Nicholas J. Higham (2011)
           "Computing the Action of the Matrix Exponential,
           with an Application to Exponential Integrators."
           SIAM Journal on Scientific Computing,
           33 (2). pp. 488-511. ISSN 1064-8275
           http://eprints.ma.man.ac.uk/1591/

    .. [2] Nicholas J. Higham and Awad H. Al-Mohy (2010)
           "Computing Matrix Functions."
           Acta Numerica,
           19. 159-208. ISSN 0962-4929
           http://eprints.ma.man.ac.uk/1451/

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import expm, expm_multiply
    >>> A = csc_matrix([[1, 0], [0, 1]])
    >>> A.toarray()
    array([[1, 0],
           [0, 1]], dtype=int64)
    >>> B = np.array([np.exp(-1.), np.exp(-2.)])
    >>> B
    array([ 0.36787944,  0.13533528])
    >>> expm_multiply(A, B, start=1, stop=2, num=3, endpoint=True)
    array([[ 1.        ,  0.36787944],
           [ 1.64872127,  0.60653066],
           [ 2.71828183,  1.        ]])
    >>> expm(A).dot(B)                  # Verify 1st timestep
    array([ 1.        ,  0.36787944])
    >>> expm(1.5*A).dot(B)              # Verify 2nd timestep
    array([ 1.64872127,  0.60653066])
    >>> expm(2*A).dot(B)                # Verify 3rd timestep
    array([ 2.71828183,  1.        ])
    """
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(
            f"Expected A to be like a square matrix, got {A.shape}."
        )
    if A.shape[1] != B.shape[0]:
        raise ValueError(
            f"Shapes of matrices A {A.shape} and B {B.shape} are incompatible."
        )
    if B.ndim not in (1, 2):
        raise ValueError("Expected B to be like a matrix or a vector.")

    is_single_exponent = all(arg is None for arg in (start, stop, num, endpoint))
    if traceA is None:
        is_linear_operator = isinstance(A, scipy.sparse.linalg.LinearOperator)
        if is_linear_operator:
            warn("Trace of LinearOperator not available, it will be estimated."
                 " Provide `traceA` to ensure performance.", stacklevel=2)
        # m3 is bit arbitrary choice, a more accurate trace (larger m3) might
        # speed up exponential calculation, but trace estimation is also costly
        # an educated guess would need to consider the number of time points
        m3 = 1 if is_single_exponent else 5
        traceA = traceest(A, m3=m3) if is_linear_operator else _trace(A)

    mu = traceA / A.shape[0]
    A = A - mu*_ident_like(A)

    if is_single_exponent:
        X = _expm_multiply_simple(A, B, mu)
    else:
        X, status = _expm_multiply_interval(A, B, mu, start, stop, num,
                                            endpoint)
    return X


def _expm_multiply_simple(A, B, mu, t=1.0, balance=False):
    """
    Compute the action of the matrix exponential at a single time point.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose exponential is of interest.
    B : ndarray
        The matrix to be multiplied by the matrix exponential of A.
    t : float
        A time point.
    balance : bool
        Indicates whether or not to apply balancing.

    Returns
    -------
    F : ndarray
        :math:`e^{t A} B`

    Notes
    -----
    This is algorithm (3.2) in Al-Mohy and Higham (2011).

    """
    if balance:
        raise NotImplementedError
    tol = 0.5*np.finfo(np.result_type(A.dtype, B.dtype)).eps
    is_linear_operator = isinstance(A, scipy.sparse.linalg.LinearOperator)
    A_1_norm = onenormest(A) if is_linear_operator else _exact_1_norm(A)
    if t*A_1_norm == 0:
        m_star, s = 0, 1
    else:
        ell = 2
        norm_info = LazyOperatorNormInfo(A, onenorm=A_1_norm, scale=t, t=ell)
        n0 = 1 if B.ndim == 1 else B.shape[1]
        m_star, s = _fragment_3_1(norm_info, n0, tol, ell=ell)
    return _expm_multiply_simple_core(A, B, t, mu, m_star, s, tol, balance)


def _expm_multiply_simple_core(A, B, t, mu, m_star, s, tol=None, balance=False):
    """A helper function."""
    if balance:
        raise NotImplementedError
    if tol is None:
        tol = 0.5*np.finfo(np.result_type(A.dtype, B.dtype)).eps
    F = B
    eta = np.exp(t * mu / s)
    for _ in range(s):
        c1 = _exact_inf_norm(B)
        for j in range(m_star):
            coeff = t / (s * (j + 1))
            B = coeff * A.dot(B)
            c2 = _exact_inf_norm(B)
            F = F + B
            if c1 + c2 <= tol * _exact_inf_norm(F):
                break
            c1 = c2
        F = eta * F
        B = F
    return F


# This table helps to compute bounds.
# They seem to have been difficult to calculate, involving symbolic
# manipulation of equations, followed by numerical root finding.
_theta = {
    # The first 30 values are from table A.3 of Computing Matrix Functions.
    1: 2.29e-16,
    2: 2.58e-8,
    3: 1.39e-5,
    4: 3.40e-4,
    5: 2.40e-3,
    6: 9.07e-3,
    7: 2.38e-2,
    8: 5.00e-2,
    9: 8.96e-2,
    10: 1.44e-1,
    # 11
    11: 2.14e-1,
    12: 3.00e-1,
    13: 4.00e-1,
    14: 5.14e-1,
    15: 6.41e-1,
    16: 7.81e-1,
    17: 9.31e-1,
    18: 1.09,
    19: 1.26,
    20: 1.44,
    # 21
    21: 1.62,
    22: 1.82,
    23: 2.01,
    24: 2.22,
    25: 2.43,
    26: 2.64,
    27: 2.86,
    28: 3.08,
    29: 3.31,
    30: 3.54,
    # The rest are from table 3.1 of
    # Computing the Action of the Matrix Exponential.
    35: 4.7,
    40: 6.0,
    45: 7.2,
    50: 8.5,
    55: 9.9,
}


class LazyOperatorNormInfo:
    """
    Information about an operator is lazily computed.

    The information includes the exact 1-norm of the operator,
    in addition to estimates of 1-norms of powers of the operator.
    This uses the notation of Computing the Action (2011).
    This class is specialized enough to probably not be of general interest
    outside of this module.

    """

    def __init__(self, A, onenorm=None, scale=1, t=2):
        """
        Provide the operator and some norm-related information.

        Parameters
        ----------
        A : linear operator
            The operator of interest.
        onenorm : float, optional
            The exact 1-norm of A.
        scale : scalar, optional
            If specified, return the norms of scale*A instead of A.
        t : int, optional
            A positive parameter controlling the tradeoff between accuracy
            versus time and memory usage. Larger values take longer and use
            more memory but give more accurate output.
        """
        self._A = A
        self._onenorm = onenorm
        self._d = {}
        self.scale = scale
        self.t = t

    @property
    def onenorm(self):
        """Compute the exact 1-norm."""
        if self._onenorm is None:
            self._onenorm = _exact_1_norm(self._A)
        return abs(self.scale)*self._onenorm

    def d(self, p):
        """Estimate d_p(A) ~= || A^p ||^(1/p) where ||.|| is the 1-norm."""
        if p not in self._d:
            est = onenormest(aslinearoperator(self._A)**p, t=self.t)
            self._d[p] = est**(1.0 / p)
        return abs(self.scale)*self._d[p]

    def alpha(self, p):
        """Compute max(d(p), d(p+1))."""
        return max(self.d(p), self.d(p+1))


def _compute_cost_div_m(m, p, norm_info):
    """
    A helper function for computing bounds.

    This is equation (3.10).
    It measures cost in terms of the number of required matrix products.

    Parameters
    ----------
    m : int
        A valid key of _theta.
    p : int
        A matrix power.
    norm_info : LazyOperatorNormInfo
        Information about 1-norms of related operators.

    Returns
    -------
    cost_div_m : int
        Required number of matrix products divided by m.

    """
    return int(np.ceil(norm_info.alpha(p) / _theta[m]))


def _compute_p_max(m_max):
    """
    Compute the largest positive integer p such that p*(p-1) <= m_max + 1.

    Do this in a slightly dumb way, but safe and not too slow.

    Parameters
    ----------
    m_max : int
        A count related to bounds.

    """
    sqrt_m_max = np.sqrt(m_max)
    p_low = int(np.floor(sqrt_m_max))
    p_high = int(np.ceil(sqrt_m_max + 1))
    return max(p for p in range(p_low, p_high+1) if p*(p-1) <= m_max + 1)


def _fragment_3_1(norm_info, n0, tol, m_max=55, ell=2):
    """
    A helper function for the _expm_multiply_* functions.

    Parameters
    ----------
    norm_info : LazyOperatorNormInfo
        Information about norms of certain linear operators of interest.
    n0 : int
        Number of columns in the _expm_multiply_* B matrix.
    tol : float
        Expected to be
        :math:`2^{-24}` for single precision or
        :math:`2^{-53}` for double precision.
    m_max : int
        A value related to a bound.
    ell : int
        The number of columns used in the 1-norm approximation.
        This is usually taken to be small, maybe between 1 and 5.

    Returns
    -------
    best_m : int
        Related to bounds for error control.
    best_s : int
        Amount of scaling.

    Notes
    -----
    This is code fragment (3.1) in Al-Mohy and Higham (2011).
    The discussion of default values for m_max and ell
    is given between the definitions of equation (3.11)
    and the definition of equation (3.12).

    """
    if ell < 1:
        raise ValueError('expected ell to be a positive integer')
    best_m = None
    best_s = None
    if _condition_3_13(norm_info.onenorm, n0, m_max, ell):
        for m, theta in _theta.items():
            s = int(np.ceil(norm_info.onenorm / theta))
            if best_m is None or m * s < best_m * best_s:
                best_m = m
                best_s = s
    else:
        # Equation (3.11).
        for p in range(2, _compute_p_max(m_max) + 1):
            for m in range(p*(p-1)-1, m_max+1):
                if m in _theta:
                    s = _compute_cost_div_m(m, p, norm_info)
                    if best_m is None or m * s < best_m * best_s:
                        best_m = m
                        best_s = s
        best_s = max(best_s, 1)
    return best_m, best_s


def _condition_3_13(A_1_norm, n0, m_max, ell):
    """
    A helper function for the _expm_multiply_* functions.

    Parameters
    ----------
    A_1_norm : float
        The precomputed 1-norm of A.
    n0 : int
        Number of columns in the _expm_multiply_* B matrix.
    m_max : int
        A value related to a bound.
    ell : int
        The number of columns used in the 1-norm approximation.
        This is usually taken to be small, maybe between 1 and 5.

    Returns
    -------
    value : bool
        Indicates whether or not the condition has been met.

    Notes
    -----
    This is condition (3.13) in Al-Mohy and Higham (2011).

    """
    # This is the rhs of equation (3.12).
    p_max = _compute_p_max(m_max)
    a = 2 * ell * p_max * (p_max + 3)

    # Evaluate the condition (3.13).
    b = _theta[m_max] / (n0 * m_max)
    return A_1_norm <= a * b


def _expm_multiply_interval(A, B, mu, start=None, stop=None, num=None,
                            endpoint=None, balance=False,
                            status_only=False):
    """
    Compute the action of the matrix exponential at multiple time points.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose exponential is of interest.
    B : ndarray
        The matrix to be multiplied by the matrix exponential of A.
    start : scalar, optional
        The starting time point of the sequence.
    stop : scalar, optional
        The end time point of the sequence, unless `endpoint` is set to False.
        In that case, the sequence consists of all but the last of ``num + 1``
        evenly spaced time points, so that `stop` is excluded.
        Note that the step size changes when `endpoint` is False.
    num : int, optional
        Number of time points to use.
    endpoint : bool, optional
        If True, `stop` is the last time point. Otherwise, it is not included.
    balance : bool
        Indicates whether or not to apply balancing.
    status_only : bool
        A flag that is set to True for some debugging and testing operations.

    Returns
    -------
    F : ndarray
        :math:`e^{t_k A} B`
    status : int
        An integer status for testing and debugging.

    Notes
    -----
    This is algorithm (5.2) in Al-Mohy and Higham (2011).

    There seems to be a typo, where line 15 of the algorithm should be
    moved to line 6.5 (between lines 6 and 7).

    """
    if balance:
        raise NotImplementedError
    tol = 0.5*np.finfo(np.result_type(A.dtype, B.dtype)).eps

    # Get the linspace samples, attempting to preserve the linspace defaults.
    linspace_kwargs = {'retstep': True}
    if num is not None:
        linspace_kwargs['num'] = num
    if endpoint is not None:
        linspace_kwargs['endpoint'] = endpoint
    samples, step = np.linspace(start, stop, **linspace_kwargs)

    # Convert the linspace output to the notation used by the publication.
    nsamples = len(samples)
    if nsamples < 2:
        raise ValueError('at least two time points are required')
    q = nsamples - 1
    h = step
    t_0 = samples[0]
    t_q = samples[q]

    # Define the output ndarray.
    # Use an ndim=3 shape, such that the last two indices
    # are the ones that may be involved in level 3 BLAS operations.
    X_shape = (nsamples,) + B.shape
    X = np.empty(X_shape, dtype=np.result_type(A.dtype, B.dtype))
    t = t_q - t_0
    is_linear_operator = isinstance(A, scipy.sparse.linalg.LinearOperator)
    A_1_norm = onenormest(A) if is_linear_operator else _exact_1_norm(A)
    ell = 2
    norm_info = LazyOperatorNormInfo(A, onenorm=A_1_norm, scale=t, t=ell)
    n0 = 1 if B.ndim == 1 else B.shape[1]
    if t*A_1_norm == 0:
        m_star, s = 0, 1
    else:
        m_star, s = _fragment_3_1(norm_info, n0, tol, ell=ell)

    # Compute the expm action up to the initial time point.
    X[0] = _expm_multiply_simple_core(A, B, t_0, mu, m_star, s)

    # Compute the expm action at the rest of the time points.
    if q <= s:
        if status_only:
            return 0
        return _expm_multiply_interval_core_0(A, X, h, mu,
                                              q, norm_info, tol, ell, n0)
    elif not (q % s):
        if status_only:
            return 1
        return _expm_multiply_interval_core_1(A, X, h, mu, m_star, s, q, tol)
    elif (q % s):
        if status_only:
            return 2
        return _expm_multiply_interval_core_2(A, X, h, mu, m_star, s, q, tol)

    raise Exception('internal error')


def _expm_multiply_interval_core_0(A, X, h, mu, q, norm_info, tol, ell, n0):
    """Compute the case q <= s."""
    # Compute the new values of m_star and s which should be applied
    # over intervals of size t/q
    if norm_info.onenorm == 0:
        m_star, s = 0, 1
    else:
        norm_info.scale /= q
        m_star, s = _fragment_3_1(norm_info, n0, tol, ell=ell)
        norm_info.scale *= q

    for k in range(q):
        X[k+1] = _expm_multiply_simple_core(A, X[k], h, mu, m_star, s)
    return X, 0


def _expm_multiply_interval_core_1(A, X, h, mu, m_star, s, q, tol):
    """Compute the case q > s and q % s == 0."""
    d = q // s
    K = np.empty((m_star + 1, *X.shape[1:]), dtype=X.dtype)
    for i in range(s):
        _compute_core(A, X, K, d, d, i, h, mu, m_star, tol)
    return X, 1


def _expm_multiply_interval_core_2(A, X, h, mu, m_star, s, q, tol):
    """Compute the case q > s and q % s > 0."""
    d = q // s
    j = q // d
    r = q - d*j
    K = np.empty((m_star + 1, *X.shape[1:]), dtype=X.dtype)
    for i in range(j + 1):
        effective_d = d if i < j else r
        _compute_core(A, X, K, d, effective_d, i, h, mu, m_star, tol)
    return X, 2


def _compute_core(A, X, K, d, effective_d, i, h, mu, m_star, tol):
    """Cumpute common calculations of cores 1 and 2."""
    K[0] = X[i*d]
    high_p = 0  # keep track, how far we calculated the dot products
    for k in range(1, effective_d+1):
        F = K[0]
        c1 = _exact_inf_norm(F)
        for p in range(1, m_star+1):
            if p == high_p + 1:  # haven't pre-calcualted dot
                K[p] = h / p * A.dot(K[p-1])
                high_p = p
            # in case the python integer type is to large we need numpy float
            coeff = K.dtype.type(pow(k, p))
            F = F + coeff*K[p]
            c2 = coeff * _exact_inf_norm(K[p])
            if c1 + c2 <= tol * _exact_inf_norm(F):
                break
            c1 = c2
        X[k + i*d] = np.exp(k*h*mu) * F
