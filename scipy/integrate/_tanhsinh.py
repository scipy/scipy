from dataclasses import dataclass
import numpy as np
from scipy import special

# todo:
#  fix tests broken by jumpstart
#  return lower-level integral estimates for testing?
#  callback
#  respect data types
#  apply np.vectorize as needed?
#  remove maxiter?
#  accept args, kwargs?
#  tests - test minweight
#  decide when output of log-integration is complex vs real
#  debug inflated error estimate with log-integration when complex part is nonzero
#  support complex integration
#  support singularities? interval subdivision? this feature will be added
#    eventually, but do we adjust the interface now?
#  warn (somehow) when invalid function values & weight < minweight
#  vectorize
#  make public?

@dataclass
class QuadratureResult:
    integral: float
    error: float
    feval: int
    success: bool
    status: int
    message: str


_status_messages = {-1: "Iteration in progress.",
                    0: ("The algorithm completed successfully, and the error "
                        "estimate meets the requested tolerance."),
                    1: ("The error estimate does not meet the specified "
                        "tolerance, but performing additional iterations "
                        "cause the function evaluation limit to be exceeded."),
                    2: ("The error estimate does not meet the specified "
                        "tolerance, but performing additional iterations "
                        "would cause the iteration limit to be exceeded."),
                    3: ("An invalid value (e.g. overflow, NaN) was "
                        "encountered within the integration interval. See "
                        "documentation notes for more information.")
                    }


def _tanhsinh(f, a, b, *, maxfun=5000, maxiter=10, atol=None, rtol=None,
              minweight=1e-100, log=False):
    """Evaluate a convergent integral numerically using tanh-sinh quadrature.

    In practice, tanh-sinh quadrature achieves quadratic convergence for
    many integrands: the number of accurate *digits* scales roughly linearly
    with the number of function evaluations [1]_.

    Either or both of the limits of integration may be infinite, and
    singularities at the endpoints are acceptable. Divergent integrals and
    integrands with non-finite derivatives or singularities within an interval
    are out of scope, but may be evaluated be calling `_tanhsinh` on each
    sub-interval separately.

    Parameters
    ----------
    f : callable
        The function to be integrated. `f` must accept a single argument,
        the value at which it is to be evaluated. `f` must accept an array
        as input and evaluate the integrand elementwise.
    a, b : float
        Lower and upper limits of integration.
    maxfun, maxiter: int, default=5000
        The maximum acceptable number of function evaluations (default: 5000)
        and iterations (default: 10), respectively. Note that the number of
        function evaluations is counted as the number of elements at which `f`
        is evaluated; not the number of calls of `f`.
        In the first iteration, `f` is called once (twice for doubly-infinite
        integrals), performing 14 (28) function evaluations. The number of
        function in each subsequent iteration is approximately double the
        number in the previous iteration.
        The function will terminate with nonzero exit status before *either* of
        these limits is reached. For an increase of one of these parameters to
        have an effect, *both* of these values must be increased.
    atol, rtol : float, optional
        Absolute termination tolerance (default: 0) and relative termination
        tolerance (default: 1e-14), respectively. The error estimate is as
        described in [1]_ Section 5; while not theoretically rigorous or
        conservative, it is said to work well in practice. Must be non-negative
        and finite if `log` is False, and must be expressed as the log of a
        non-negative and finite number if `log` is True.
    minweight : float, default: 1e-100
        The minimum nonzero weight to be used in the Euler-Maclaurin Summation
        formula as described in [1]_ Section 4. When evaluating an integral
        with a singularity at an endpoint, values of the integrand will not be
        considered if the weight prescribed by the tanh-sinh quadrature scheme
        falls below this threshold. The default value is 1e-100. For integrals
        without endpoint singularities, smaller values are acceptable. For
        integrals with an endpoint singularity, a larger value will protect
        against overflows of the integrand, but potentially at the expense of
        integral and error estimate accuracy.
    log : bool, default: True
        When set to True, func returns the log of the integrand, `atol` and
        `rtol` must be expressed as the logs of the absolute and relative
        errors, and the result object contains the log of the integral and
        error.

    Returns
    -------
    res : QuadratureResult
        An object with the following attributes.

        integral : float
            An estimate of the integral
        error : float
            An estimate of the error.
        feval : int
            The number of function evaluations, i.e., the number of
            points at which the integrand was evaluated.
        success : bool
            Whether the algorithm terminated successfully with an error
            estimate satisfying the specified `atol` or `rtol`.
        status : int
            A numerical code indicating the status of the algorithm at
            termination. ``0`` indicates successful termination. ``1`` and
            ``2`` indicate termination due to function or iteration limit
            (respectively), and ``3`` indicates that an invalid value
            was encountered within the integration interval.
        message : str
            A description of the termination status.

    See Also
    --------
    quad, quadrature

    Notes
    -----
    Implements the algorithm as described in [1]_ with minor adaptations for
    fixed-precision arithmetic. The tanh-sinh scheme was originally introduced
    in [2]_.

    Due to floating-point error in the abscissae, the function may be evaluated
    at the endpoints of the interval, but the value returned will be ignored.

    If the function returns an invalid value (e.g. infinity, NaN) within the
    integral and the corresponding weight is greater than `minweight`,
    `_tanhsinh` will terminate with exit status ``3``. If this is caused by
    an interior singularity, break the interval into sub-intervals with
    the singularity at an endpoint, and call `_tanhsinh` on each.

    References
    ----------
    [1] Bailey, David H., Karthik Jeyabalan, and Xiaoye S. Li. "A comparison of
        three high-precision quadrature schemes." Experimental Mathematics 14.3
        (2005): 317-329.
    [2] Takahasi, Hidetosi, and Masatake Mori. "Double exponential formulas for
        numerical integration." Publications of the Research Institute for
        Mathematical Sciences 9.3 (1974): 721-741.

    Example
    -------
    Evaluate the Gaussian integral:

    >>> import numpy as np
    >>> from scipy import integrate
    >>> def f(x):
    ...     return np.exp(-x**2)
    >>> res = _tanhsinh(f, -np.inf, np.inf)
    >>> res.integral  # true value is np.sqrt(np.pi), 1.7724538509055159
    1.772453850905516
    >>> res.error  # actual error is ~2.220446049250313e-16
    1.062283726062058e-15

    The value of the Gaussian function (bell curve) is nearly zero for
    arguments sufficiently far from zero, so the value of the integral
    over a finite interval is nearly the same.

    >>> _tanhsinh(f, -20, 20).integral
    1.7724538509055152

    However, with unfavorable integration limits, the integration scheme
    may not be able to find the important region.

    >>> _tanhsinh(f, -np.inf, 20).integral
    1.967881377548732e-19

    In such cases, or when there are singularities within the interval,
    break the integral into parts with endpoints at the important points.

    >>> _tanhsinh(f, -np.inf, 0).integral + _tanhsinh(f, 0, 20).integral
    1.7724538509055163

    """
    # Input validation and standardization
    res = _tanhsinh_iv(f, a, b, maxfun, maxiter, atol, rtol, minweight, log)
    f, a, b, maxfun, maxiter, atol, rtol, minweight, log, feval_factor = res

    # Initialization
    Sk = []  # sequence of integral estimates for error estimation
    Sn = aerr = np.nan  # integral and error are NaN until determined otherwise
    status = -1  # "Iteration in progress." for callback
    feval = 0  # function evaluation counter

    for n in range(0, maxiter):  # for each "level"
        h = 1 / 2**n  # step size

        # retrieve abscissa-weight pairs from cache
        xjc, wj = _get_pairs(n, inclusive=(n==0))

        # Determine whether evaluating the function at the abscissae will
        # cause the function evaluation limit to be exceeded. If so, break.
        # Otherwise, increment the function eval counter.
        if feval + len(xjc) * feval_factor > maxfun:
            status = 1
            break

        # Transform the abscissae as required by the limits of integration
        xj, wj = _transform_to_limits(a, b, xjc, wj)

        # Perform the Euler-Maclaurin sum
        feval += 2 * len(xjc) * feval_factor  # function evals happen next
        fjwj, Sn = _euler_maclaurin_sum(f, h, xj, wj, minweight, log)

        # Check for infinities / NaNs. If encountered, break.
        not_finite = ((log and np.any(np.isposinf(np.real(fjwj)) | np.isnan(fjwj))) or
                      (not log and not np.all(np.isfinite(fjwj))))
        if not_finite:
            status = 3  # function evaluation limit
            break

        # If we have integral estimates from a lower level, update it.
        if Sk:
            Snm1 = Sk[-1]
            Sn = (special.logsumexp([Snm1 - np.log(2), Sn]) if log
                  else Snm1 / 2 + Sn)
        # Otherwise, the integral estimate is just Sn.

        # Check error estimate (see [1] 5. Error Estimation, page 11)
        rerr, aerr = _error_estimate(h, n, Sn, Sk, fjwj, log)
        success = (rerr < rtol or (rerr + Sn < atol) if log
                   else rerr < rtol or rerr*abs(Sn) < atol)
        if success:
            status = 0  # Success
            break

        # Store integral estimate.
        Sk.append(Sn)
    else:
        status = 2  # Iteration limit

    message = _status_messages[status]
    return QuadratureResult(integral=Sn, error=aerr, feval=feval,
                            success=status==0, status=status, message=message)


def _compute_pair(k):
    # Compute the abscissa-weight pairs for each level m. See [1] page 9.

    # "....each level k of abscissa-weight pairs uses h = 2 **-k"
    h = 1 / (2 ** k)

    # "We find that roughly 3.6 * 2^k abscissa-weight pairs are generated at
    # "level k." The actual number per level can be generated like:
    # for i in range(10):
    #     _, xjc, wj = _compute_pair(i)
    #     # don't want infinite weights or to evaluate f at endpoints
    #     valid = (xjc > 0) & (wj > 0) & np.isfinite(wj)
    #     print(np.sum(valid))
    # Running this code, I'm finding that the maximum index value w/ 64-bit is:
    max = int(np.ceil(6.115 * 2**k))
    # This reproduces all the integers produced by the loop above for k <= 10.
    # Note that the actual number of pairs is *half* of this (see below).

    # For iterations after the first, "....the integrand function needs to be
    # evaluated only at the odd-indexed abscissas at each level."
    j = np.arange(max) if k == 0 else np.arange(1, max, 2)
    jh = j * h

    # "In this case... the weights wj = u1/cosh(u2)^2, where..."
    pi_2 = np.pi / 2
    u1 = pi_2*np.cosh(jh)
    u2 = pi_2*np.sinh(jh)
    wj = u1 / np.cosh(u2)**2

    # "We actually store 1-xj = 1/(...)."
    xjc = 1 / (np.exp(u2) * np.cosh(u2))  # complement of xj = np.tanh(u2)

    # When level k == 0, the zeroth xj corresponds with xj = 0. To simplify
    # code, the function will be evaluated there twice; each gets half weight.
    wj[0] = wj[0] / 2 if k == 0 else wj[0]

    return xjc, wj


def _pair_cache(max_level):
    # Cache the ascissa-weight pairs up to a specified level
    # All abscissae (weights) are stored concatenated in a 1D array;
    # `index` notes the last index of each level.
    pairs = [_compute_pair(k) for k in range(max_level + 1)]
    lengths = [len(pair[0]) for pair in pairs]
    _pair_cache.xjc, _pair_cache.wj = np.concatenate(pairs, axis=-1)
    _pair_cache.indices = np.cumsum(lengths)
_pair_cache.xjc = None
_pair_cache.wj = None
_pair_cache.indices = []


def _get_pairs(k, inclusive=False):
    # Retrieve the specified abscissa-weight pairs from the cache
    # If `inclusive`, return all up to and including the specified level
    indices = _pair_cache.indices
    if len(indices) <= k:
        _pair_cache(max(10, k))
        indices = _pair_cache.indices
    xjc = _pair_cache.xjc
    wj = _pair_cache.wj

    start = 0 if (k == 0 or inclusive) else indices[k-1] - 1
    end = indices[k]

    return xjc[start:end].copy(), wj[start:end].copy()


def _transform_to_limits(a, b, xjc, wj):
    # Transform integral according to user-specified limits. This is just
    # math that follows from the fact that the standard limits are (-1, 1).
    # Note: If we had stored xj instead of xjc, we would have
    # xj = alpha * xj + beta, where beta = (a + b)/2
    alpha = (b - a) / 2
    xj = np.concatenate((-alpha * xjc + b, alpha * xjc + a), axis=-1)
    wj = wj*alpha  # these need to get broadcasted, so no *=
    wj = np.concatenate((wj, wj), axis=-1)

    # Points at the boundaries can be generated due to finite precision
    # arithmetic. Ideally we wouldn't evaluate the function at these points
    # or when weights are zero; however, it would be tricky to filter out
    # points when this function is vectorized. Instead, zero the weights.
    invalid = (xj <= a) | (xj >= b)
    wj[invalid] = 0
    return xj, wj


def _euler_maclaurin_sum(f, h, xj, wj, minweight, log):
    # Perform the Euler-Maclaurin Sum, [1] Section 4
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        # Warnings due to function evaluation at endpoints are not cause for
        # concern; the resulting values will be given 0 weight.
        fj = f(xj)
        # Zero weight multiplied by infinity results in NaN. These will be
        # removed below.
        fjwj = fj + np.log(wj) if log else fj * wj

    # Points on the boundaries can be generated due to finite precision
    # arithmetic. Ideally we wouldn't evaluate the function at these points
    # or when weights are zero; however, it would be tricky to filter out
    # points when this function is vectorized. Set the results to zero.
    invalid = (wj <= minweight)
    replacement = -np.inf if log else 0
    fjwj[invalid] = replacement

    # update integral estimate
    Sn = (special.logsumexp(fjwj + np.log(h), axis=-1) if log
          else np.sum(fjwj, axis=-1) * h)

    return fjwj, Sn


def _error_estimate(h, n, Sn, Sk, fjwj, log):
    # Estimate the error according to [1] Section 5

    # The paper says error 1 for n<=2, but I disagree. We start at level 0, so
    # there is no reason not to compute the error beginning in level 2. Before
    # then, the error cannot be calculated, so use NaN, not 1.
    if n <= 1:
        return np.nan, np.nan

    # With a jump start (starting at level higher than 0), we haven't
    # explicitly calculated the integral estimate at lower levels. But we have
    # all the function value-weight products, so we can compute the
    # lower-level estimates.
    if len(Sk) < 2:
        indices = _pair_cache.indices
        hm1 = 2 * h
        hm2 = 4 * h
        z = fjwj.reshape(2, -1)
        fjwjm1 = z[..., :indices[n-1]]
        fjwjm2 = z[..., :indices[n-2]]
        Snm1 = (special.logsumexp(fjwjm1 + np.log(hm1)) if log
                else np.sum(fjwjm1) * hm1)
        Snm2 = (special.logsumexp(fjwjm2 + np.log(hm2)) if log
                else np.sum(fjwjm2) * hm2)
    else:
        Snm2, Snm1 = Sk[-2:]

    e1 = np.finfo(np.float64).eps

    if log:
        log_e1 = np.log(e1)
        # Currently, only real integrals are supported. All complex values
        # have imaginary part in increments of pi*j, which just carries sign
        # information. Use of `np.real` here is equivalent to absolute value
        # in real scale.
        d1 = np.real(special.logsumexp([Sn, Snm1 + np.pi*1j]))
        d2 = np.real(special.logsumexp([Sn, Snm2 + np.pi*1j]))
        d3 = log_e1 + np.max(np.real(fjwj))
        d4 = np.max(np.real(np.reshape(fjwj, (2, -1))[:, -1]))
        rerr = np.max([d1 ** 2 / d2, 2 * d1, d3, d4])
        aerr = max(log_e1, rerr) + np.real(Sn)
    else:
        # Note: explicit computation of log10 of each of these is unnecessary.
        d1 = np.abs(Sn - Snm1)
        d2 = np.abs(Sn - Snm2)
        d3 = e1 * np.max(np.abs(fjwj))
        d4 = np.max(np.reshape(np.abs(fjwj), (2, -1))[:, -1])
        rerr = np.max([d1**((np.log(d1)/np.log(d2))), d1**2, d3, d4])
        aerr = max(e1, rerr) * np.abs(Sn)
    return rerr, aerr


def _tanhsinh_iv(f, a, b, maxfun, maxiter, atol, rtol, minweight, log):
    # Input validation and standardization

    message = '`f` must be callable.'
    if not callable(f):
        raise ValueError(message)

    message = '`log` must be True or False.'
    if log not in {True, False}:
        raise ValueError(message)
    log = bool(log)

    if atol is None:
        atol = -np.inf if log else 0

    if rtol is None:
        rtol = np.log(1e-14) if log else 1e-14

    message = ('Integration limits `a` and `b`, tolerances `atol` and '
               '`rtol`, and `minweight` must be reals.')
    params = np.asarray([a, b, atol, rtol, minweight])
    if not np.issubdtype(params.dtype, np.floating) or np.any(np.isnan(params)):
        raise ValueError(message)
    if log:
        message = '`atol` and `rtol` may not be positive infinity.'
        if np.any(np.isposinf(params[-3:-1])):
            raise ValueError(message)
    else:
        message = '`atol` and `rtol` must be non-negative and finite.'
        if np.any(params[-3:-1] < 0) or np.any(np.isinf(params[-3:-1])):
            raise ValueError(message)
    message = '`minweight` must be positive and finite.'
    if np.any(params[-1] <= 0) or np.any(np.isinf(params[-1])):
        raise ValueError(message)
    a, b, atol, rtol, minweight = params

    message = '`maxfun` and `maxiter` must be integers.'
    params = np.asarray([maxfun, maxiter])
    if not np.issubdtype(params.dtype, np.integer):
        raise ValueError(message)
    message = '`maxfun` and `maxiter` must be positive.'
    if np.any(params <= 0):
        raise ValueError(message)
    maxfun, maxiter = params

    # Transform integrals as needed for infinite limits
    # There are more efficient ways of doing this to avoid function call
    # overhead, but let's stick with the simplest for now.
    if b < a:
        def f(x, f=f):
            return (f(x) + np.pi*1j if log
                    else -f(x))
        a, b = b, a

    if np.isinf(a) and np.isinf(b):
        def f(x, f=f):
            return (special.logsumexp([f(x), f(-x)], axis=0) if log
                    else f(x) + f(-x))
        a, b = 0, np.inf
        feval_factor = 2  # user function evaluated twice each call
    elif np.isinf(a):
        def f(x, f=f):
            return f(-x)
        a, b = -b, -a
        feval_factor = 1
    else:
        feval_factor = 1

    if np.isinf(b):
        def f(x, f=f, a=a):
            return (f(1/x - 1 + a) - 2*np.log(x+0j) if log
                    else f(1/x - 1 + a)*x**-2)
        a, b = 0, 1

    return f, a, b, maxfun, maxiter, atol, rtol, minweight, log, feval_factor
