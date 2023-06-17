# mypy: disable-error-code="attr-defined"
from dataclasses import dataclass
import numpy as np
from scipy import special

# todo:
#  without `minweight`, we are also suppressing infinities within the interval.
#    Is that OK? If so, we can probably get rid of `status=3`.
#  Add heuristic to stop when improvement is too slow
#  special case for equal limits of integration?
#  callback - test results at lower levels with higher maxlevel against
#             final results at lower maxlevel
#  accept args, kwargs?
#  support singularities? interval subdivision? this feature will be added
#    eventually, but do we adjust the interface now?
#  vectorize, respecting data type
#  When doing log-integration, should the tolerances control the error of the
#    log-integral or the error of the integral?  The trouble is that `log`
#    inherently looses some precision so it may not be possible to refine
#    the integral further. Example: 7th moment of stats.f(15, 20)
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
                        "tolerance, but performing an additional iteration "
                        "would cause the maximum level to be exceeded."),
                    2: ("The error estimate does not meet the specified "
                        "tolerance, but performing an additional iteration "
                        "cause the function evaluation limit to be exceeded."),
                    3: ("An invalid value (e.g. overflow, NaN) was "
                        "encountered within the integration interval. See "
                        "documentation notes for more information.")
                    }

def _tanhsinh(f, a, b, *, log=False, maxfun=None, maxlevel=None, minlevel=2,
              atol=None, rtol=None):
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
    log : bool, default: False
        Setting to True indicates that `func` returns the log of the integrand
        and that `atol` and `rtol` are expressed as the logs of the absolute
        and relative errors. In this case, the result object will contain the
        log of the integral and error. This is useful for integrands for which
        numerical underflow or overflow would lead to inaccuracies.
    maxfun, maxlevel : int or np.inf, optional
        The maximum acceptable number of function evaluations and refinement
        level, respectively. Note that the number of function evaluations is
        counted as the number of elements at which `f` is evaluated; not the
        number of calls of `f`.

        At the zeroth level, `f` is called once (twice for doubly-infinite
        integrals), performing 14 (28) function evaluations. Each subsequent
        level approximately doubles the number of function evaluations
        performed. Accordingly, for many integrands, each successive level will
        double the number of accurate digits in the result (up to the limits of
        floating point precision).

        - (default) If neither `maxfun` nor `maxlevel` is provided, the
          algorithm will terminate after completing refinement level 10.
        - If only `maxfun` (`maxlevel`) is provided, the function will terminate
          with exit status ``1`` (``2``) if proceeding to the next refinement
          level would cause the specified limit to be exceeded.
        - If both `maxfun` and `maxlevel` are provided, the function will
          terminate with nonzero exit status before *either* of these limits is
          reached.

    minlevel : int, optional
        The level at which to begin iteration (default: 0). This does not
        change the total number of function evaluations or the points at
        which the function is evaluated; it changes only the *number of times*
        `f` is called. If ``minlevel=k``, then the integrand is evaluated at
        all points from levels ``0`` through ``k`` in the first call.
    atol, rtol : float, optional
        Absolute termination tolerance (default: 0) and relative termination
        tolerance (default: 1e-12), respectively. The error estimate is as
        described in [1]_ Section 5; while not theoretically rigorous or
        conservative, it is said to work well in practice. Must be non-negative
        and finite if `log` is False, and must be expressed as the log of a
        non-negative and finite number if `log` is True.

    Returns
    -------
    res : QuadratureResult
        An object with the following attributes.

        integral : float
            An estimate of the integral
        error : float
            An estimate of the error. Only available if level two or higher
            has been evaluated; otherwise NaN.
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
    fixed-precision arithmetic, including some described by [2]_ and [3]_. The
    tanh-sinh scheme was originally introduced in [4]_.

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
    [2] Vanherck, Joren, Bart Sorée, and Wim Magnus. "Tanh-sinh quadrature for
        single and multiple integration using floating-point arithmetic."
        arXiv preprint arXiv:2007.15057 (2020).
    [3] van Engelen, Robert A.  "Improving the Double Exponential Quadrature
        Tanh-Sinh, Sinh-Sinh and Exp-Sinh Formulas."
        https://www.genivia.com/files/qthsh.pdf
    [4] Takahasi, Hidetosi, and Masatake Mori. "Double exponential formulas for
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
    res = _tanhsinh_iv(f, a, b, log, maxfun, maxlevel, minlevel, atol, rtol)
    (f, a, b, log, maxfun, maxlevel, minlevel, atol, rtol) = res

    return _tanhsinh_noiv(f, a, b, log, maxfun, maxlevel, minlevel, atol, rtol)

def _tanhsinh_noiv(f, a, b, log, maxfun, maxlevel, minlevel, atol, rtol):
    # Transform improper integrals
    f, a, b, feval_factor = _transform_integrals(f, a, b, log)

    # Initialization
    Sk = []  # sequence of integral estimates for error estimation
    # Most extreme abscissae and corresponding `fj`s, `wj`s in Euler-Maclaurin
    # sum. These are dummy values that will be replaced in the first iteration.
    last_terms = (np.inf, np.nan, 0, -np.inf, np.nan, 0, 0)
    Sn = aerr = np.nan  # integral and error are NaN until determined otherwise
    status = -1  # "Iteration in progress." for callback
    feval = 0  # function evaluation counter
    h0 = _get_base_step()

    for n in range(minlevel, maxlevel+1):  # for each "level"
        h = h0 / 2**n  # step size

        # retrieve abscissa-weight pairs from cache
        xjc, wj = _get_pairs(n, h0, inclusive=(n==minlevel))

        # Determine whether evaluating the function at the abscissae will
        # cause the function evaluation limit to be exceeded. If so, break.
        # Otherwise, increment the function eval counter. `feval_factor` is
        # needed for doubly-infinite integrals (see `_transform_integral`).
        if feval + 2*len(xjc) * feval_factor > maxfun:
            status = 2  # function evaluation limit
            break

        # Transform the abscissae as required by the limits of integration
        xj, wj = _transform_to_limits(xjc, wj, a, b)

        # Perform the Euler-Maclaurin sum
        feval += len(xj) * feval_factor  # function evals happen next
        fjwj, Sn, last_terms = _euler_maclaurin_sum(f, xj, wj, h, last_terms, log)

        # Check for infinities / NaNs. If encountered, break.
        not_finite = ((log and (np.isposinf(np.real(Sn)) or np.isnan(Sn))) or
                      (not log and not np.isfinite(Sn)))
        if not_finite:
            status = 3  # invalid value encountered
            break

        # If we have an integral estimate from a lower level, `Sn` calculated
        # by _euler_maclaurin_sum is an update.
        if Sk:
            Snm1 = Sk[-1]
            Sn = (special.logsumexp([Snm1 - np.log(2), Sn]) if log
                  else Snm1 / 2 + Sn)
        # Otherwise, the integral estimate is just Sn.

        # Check error estimate (see [1] 5. Error Estimation, page 11)
        rerr, aerr, Sk = _estimate_error(n, Sn, Sk, fjwj, h,
                                         last_terms, log)
        success = (rerr < rtol or (rerr + np.real(Sn) < atol) if log
                   else rerr < rtol or rerr*abs(Sn) < atol)
        if success:
            status = 0  # Success
            break

        # Store integral estimate.
        Sk.append(Sn)
    else:
        status = 1  # Iteration limit

    message = _status_messages[status]
    return QuadratureResult(integral=Sn, error=aerr, feval=feval,
                            success=status==0, status=status, message=message)


# argument `dtype` currently unused, but will be used
# when I vectorize and respect dtypes.
def _get_base_step(dtype=np.float64):
    # Compute the base step length for the provided dtype. Theoretically, the
    # Euler-Maclaurin sum is infinite, but it gets cut off when either the
    # weights underflow or the abscissae cannot be distinguished from the
    # limits of integration. The latter happens to occur first for float32 and
    # float64, and it occurs when`xjc` (the abscissa complement)
    # in `_compute_pair` underflows. We can solve for the argument at which
    # it will underflow `tmax` using [2] Eq. 13.
    fmin = np.finfo(dtype).tiny*4  # stay a little away from the limit
    tmax = np.arcsinh(np.log(2/fmin - 1)/np.pi)

    # Based on this, we can choose a base step size `h` for level 0.
    # The number of function evaluations will be `2 + m*2^(k+1)`, where `k` is
    # the level and `m` is an integer we get to choose. I choose `8` somewhat
    # arbitrarily, but a rationale is that a power of 2 makes floating point
    # arithmetic more predictable. It also results in a base step size close
    # to `1`, which is what [1] uses (and I used here until I found [2] and
    # these ideas settled).
    h0 = tmax / _N_BASE_STEPS
    return h0


_N_BASE_STEPS = 8


def _compute_pair(k, h0):
    # Compute the abscissa-weight pairs for each level k. See [1] page 9.

    # "....each level k of abscissa-weight pairs uses h = 2 **-k"
    # We adapt to floating point arithmetic using ideas of [2].
    h = h0 / 2**k
    max = _N_BASE_STEPS * 2**k

    # For iterations after the first, "....the integrand function needs to be
    # evaluated only at the odd-indexed abscissas at each level."
    j = np.arange(max+1) if k == 0 else np.arange(1, max+1, 2)
    jh = j * h

    # "In this case... the weights wj = u1/cosh(u2)^2, where..."
    pi_2 = np.pi / 2
    u1 = pi_2*np.cosh(jh)
    u2 = pi_2*np.sinh(jh)
    # Denominators get big here. Overflow then underflow doesn't need warning.
    with np.errstate(under='ignore', over='ignore'):
        wj = u1 / np.cosh(u2)**2
        # "We actually store 1-xj = 1/(...)."
        xjc = 1 / (np.exp(u2) * np.cosh(u2))  # complement of xj = np.tanh(u2)

    # When level k == 0, the zeroth xj corresponds with xj = 0. To simplify
    # code, the function will be evaluated there twice; each gets half weight.
    wj[0] = wj[0] / 2 if k == 0 else wj[0]

    return xjc, wj


def _pair_cache(max_level, h0):
    # Cache the ascissa-weight pairs up to a specified level
    # All abscissae (weights) are stored concatenated in a 1D array;
    # `index` records the first index of each level after level 0.
    pairs = [_compute_pair(k, h0) for k in range(max_level + 1)]
    _pair_cache.xjc, _pair_cache.wj = np.concatenate(pairs, axis=-1)
    lengths = [len(pair[0]) for pair in pairs]
    _pair_cache.indices = np.cumsum(lengths)
_pair_cache.xjc = None
_pair_cache.wj = None
_pair_cache.indices = []


def _get_pairs(k, h0, inclusive=False):
    # Retrieve the specified abscissa-weight pairs from the cache
    # If `inclusive`, return all up to and including the specified level
    indices = _pair_cache.indices
    if len(indices) <= k:
         # rarely are more than 6 levels needed; 10 is plenty to start with
        _pair_cache(max(10, k), h0)
        indices = _pair_cache.indices
    xjc = _pair_cache.xjc
    wj = _pair_cache.wj

    start = 0 if (k == 0 or inclusive) else indices[k-1]
    end = indices[k]

    return xjc[start:end], wj[start:end]


def _transform_to_limits(xjc, wj, a, b):
    # Transform integral according to user-specified limits. This is just
    # math that follows from the fact that the standard limits are (-1, 1).
    # Note: If we had stored xj instead of xjc, we would have
    # xj = alpha * xj + beta, where beta = (a + b)/2
    alpha = (b - a) / 2
    xj = np.concatenate((-alpha * xjc + b, alpha * xjc + a), axis=-1)
    wj = wj*alpha  # these need to get broadcasted, so no *=
    wj = np.concatenate((wj, wj), axis=-1)

    # Points at the boundaries can be generated due to finite precision
    # arithmetic, but these function values aren't supposed to be included in
    # the Euler-Maclaurin sum. Ideally we wouldn't evaluate the function at
    # these points; however, we can't easily filter out points when this
    # function is vectorized. Instead, zero the weights.
    invalid = (xj <= a) | (xj >= b)
    wj[invalid] = 0
    return xj, wj


def _euler_maclaurin_sum(f, xj, wj, h, last_terms, log):
    # Perform the Euler-Maclaurin Sum, [1] Section 4
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        # Warnings due to any endpoints singularities are not cause for
        # concern; the resulting values will be replaced.
        fj = f(xj)

    # The error estimate needs to know the magnitude of the last term
    # omitted from the Euler-Maclaurin sum. This is a bit involved because
    # it may have been computed at a previous level. I sure hope it's worth
    # all the trouble.
    xl0, fl0, wr0, xr0, fr0, wl0, d4 = last_terms  # incumbent last terms

    # Find the most extreme abscissae corresponding with terms that are
    # included in the sum in this level.
    invalid = ~np.isfinite(fj) | (wj == 0)  # these terms aren't in the sum
    invalid_r, invalid_l = invalid.reshape(2, -1)
    xr, xl = xj.reshape(2, -1).copy()  # this gets modified
    xr[invalid_r], xl[invalid_l] = -np.inf, np.inf
    ir, il = np.argmax(xr), np.argmin(xl)

    # Determine whether the most extreme abscissae are from this level or
    # a previous level. Update the record of the corresponding weights and
    # function values.
    fr, fl = fj.reshape(2, -1)
    wr, wl = wj.reshape(2, -1)
    xr0, fr0, wr0 = ((xr[ir], fr[ir], wr[ir]) if xr[ir] > xr0
                     else (xr0, fr0, wr0))
    xl0, fl0, wl0 = ((xl[il], fl[il], wl[il]) if xl[il] < xl0
                     else (xl0, fl0, wl0))

    # Compute the error estimate `d4` - the magnitude of the leftmost or
    # rightmost term, whichever is greater.
    flwl0 = fl0 + np.log(wl0) if log else fl0 * wl0  # leftmost term
    frwr0 = fr0 + np.log(wr0) if log else fr0 * wr0  # rightmost term
    magnitude = np.real if log else np.abs
    d4 = np.maximum(magnitude(flwl0), magnitude(frwr0))
    last_terms = xl0, fl0, wl0, xr0, fr0, wr0, d4

    # There are two approaches to dealing with function values that are
    # numerically infinite due to approaching a singularity - zero them, or
    # replace them with the function value at the nearest non-infinite point.
    # [3] pg. 22 suggests the latter, so let's do that given that we have the
    # information.
    fr[invalid_r] = fr0
    fl[invalid_l] = fl0

    # When wj is zero, log emits a warning
    with np.errstate(divide='ignore'):
        fjwj = fj + np.log(wj) if log else fj * wj

    # update integral estimate
    Sn = (special.logsumexp(fjwj + np.log(h), axis=-1) if log
          else np.sum(fjwj, axis=-1) * h)

    return fjwj, Sn, last_terms


def _estimate_error(n, Sn, Sk, fjwj, h, last_terms, log):
    # Estimate the error according to [1] Section 5

    if n == 0:
        # The paper says to use "one" as the error before it can be calculated.
        # NaN seems to be more appropriate.
        return np.nan, np.nan, Sk

    indices = _pair_cache.indices
    z = fjwj.reshape(2, -1)

    # With a jump start (starting at level higher than 0), we haven't
    # explicitly calculated the integral estimate at lower levels. But we have
    # all the function value-weight products, so we can compute the
    # lower-level estimates.
    if len(Sk) == 0:
        hm1 = 2 * h
        fjwjm1_rl = fjwj_rl[:, :indices[n-1]]
        Snm1 = (special.logsumexp(fjwjm1_rl + np.log(hm1)) if log
                else np.sum(fjwjm1_rl) * hm1)
        Sk.append(Snm1)

    if n == 1:
        return np.nan, np.nan, Sk

    # The paper says not to calculate the error for n<=2, but it's not clear
    # about whether it starts at level 0 or level 1. We start at level 0, so
    # why not compute the error beginning in level 2?
    if len(Sk) < 2:
        hm2 = 4 * h
        fjwjm2_rl = fjwj_rl[:, :indices[n-2]]
        Snm2 = (special.logsumexp(fjwjm2_rl + np.log(hm2)) if log
                else np.sum(fjwjm2_rl) * hm2)
        Sk.insert(0, Snm2)

    Snm2, Snm1 = Sk[-2:]

    e1 = np.finfo(np.float64).eps

    if log:
        log_e1 = np.log(e1)
        # Currently, only real integrals are supported in log-scale. All
        # complex values have imaginary part in increments of pi*j, which just
        # carries sign information of the original integral, so use of
        # `np.real` here is equivalent to absolute value in real scale.
        fjwj = np.real(fjwj)
        d1 = np.real(special.logsumexp([Sn, Snm1 + np.pi*1j]))
        d2 = np.real(special.logsumexp([Sn, Snm2 + np.pi*1j]))
        d3 = log_e1 + np.max(fjwj)
        d4 = last_terms[-1]
        aerr = np.max([d1 ** 2 / d2, 2 * d1, d3, d4])
        rerr = max(log_e1, aerr - np.real(Sn))
    else:
        # Note: explicit computation of log10 of each of these is unnecessary.
        fjwj = np.abs(fjwj)
        d1 = np.abs(Sn - Snm1)
        d2 = np.abs(Sn - Snm2)
        d3 = e1 * np.max(fjwj)
        d4 = last_terms[-1]
        # If `d1` is 0, no need to warn. This does the right thing.
        with np.errstate(divide='ignore'):
            aerr = np.max([d1**(np.log(d1)/np.log(d2)), d1**2, d3, d4])
        rerr = max(e1, aerr/np.abs(Sn))
    return rerr, aerr, Sk

def _tanhsinh_iv(f, a, b, log, maxfun, maxlevel, minlevel, atol, rtol):
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
        rtol = np.log(1e-12) if log else 1e-12

    message = ('Integration limits `a` and `b` and tolerances `atol` and '
               '`rtol` must be reals.')
    params = np.asarray([a, b, atol, rtol, 0.])
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
    a, b, atol, rtol, _ = params

    BIGINT = int(2**63-2)  # avoid overflow when this is incremented
    if maxfun is None and maxlevel is None:
        maxlevel = 10
    # `np.isposinf` errors for objects other than numbers. This seems to work
    # for any common representation of infinity.
    maxfun = BIGINT if (maxfun is None or maxfun == np.inf) else maxfun
    maxlevel = BIGINT if (maxlevel is None or maxlevel == np.inf) else maxlevel

    message = '`maxfun`, `maxlevel`, and `minlevel` must be integers.'
    params = np.asarray([maxfun, maxlevel, minlevel])
    if not np.issubdtype(params.dtype, np.integer):
        raise ValueError(message)
    message = '`maxfun`, `maxlevel`, and `minlevel` must be non-negative.'
    if np.any(params < 0):
        raise ValueError(message)
    maxfun, maxlevel, minlevel = params

    return (f, a, b, log, maxfun, maxlevel, minlevel, atol, rtol)

def _transform_integrals(f, a, b, log):
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
            return (f(1/x - 1 + a) - 2*np.log(abs(x)) if log
                    else f(1/x - 1 + a)*x**-2)
        a, b = 0, 1

    return f, a, b, feval_factor
