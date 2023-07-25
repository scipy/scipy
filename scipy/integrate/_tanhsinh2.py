# mypy: disable-error-code="attr-defined"
from dataclasses import dataclass
import numpy as np
from scipy import special
from scipy.optimize import OptimizeResult
from scipy.optimize._zeros_py import (_scalar_optimization_loop,
                                      _scalar_optimization_initialize)

# todo:
#  improper integral transformations
#  update input validation
#  check function evaluation limit
#  update tests
#  address https://github.com/scipy/scipy/pull/18650#discussion_r1232935669
#  address https://github.com/scipy/scipy/pull/18650#discussion_r1233032521
#  without `minweight`, we are also suppressing infinities within the interval.
#    Is that OK? If so, we can probably get rid of `status=3`.
#  Add heuristic to stop when improvement is too slow
#  callback - test results at lower levels with higher maxlevel against
#             final results at lower maxlevel
#  support singularities? interval subdivision? this feature will be added
#    eventually, but do we adjust the interface now?
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


def _pair_cache(k, h0):
    # Cache the ascissa-weight pairs up to a specified level.
    # Abscissae and weights of consecutive levels are concatenated in separate
    # 1D arrays. `index` records the indices that correspond with each level;
    # `xjc[index[k]:index[k+1]` extracts the level `k` abscissae.
    xjcs = [_pair_cache.xjc]
    wjs = [_pair_cache.wj]

    for i in range(len(_pair_cache.indices)-1, k + 1):
        xjc, wj = _compute_pair(i, h0)
        xjcs.append(xjc)
        wjs.append(wj)
        _pair_cache.indices.append(_pair_cache.indices[-1] + len(xjc))

    _pair_cache.xjc = np.concatenate(xjcs)
    _pair_cache.wj = np.concatenate(wjs)

_pair_cache.xjc = np.empty(0)
_pair_cache.wj = np.empty(0)
_pair_cache.indices = [0]


def _get_pairs(k, h0, inclusive=False):
    # Retrieve the specified abscissa-weight pairs from the cache
    # If `inclusive`, return all up to and including the specified level
    if len(_pair_cache.indices) <= k+2:
        _pair_cache(k, h0)

    xjc = _pair_cache.xjc
    wj = _pair_cache.wj
    indices = _pair_cache.indices

    start = 0 if inclusive else indices[k]
    end = indices[k+1]

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


def _euler_maclaurin_sum(fj, work):
    # Perform the Euler-Maclaurin Sum, [1] Section 4

    # The error estimate needs to know the magnitude of the last term
    # omitted from the Euler-Maclaurin sum. This is a bit involved because
    # it may have been computed at a previous level. I sure hope it's worth
    # all the trouble.
    xl0, fl0, wr0, xr0, fr0, wl0, d4 = work.xl0, work.fl0, work.wr0, work.xr0, work.fr0, work.wl0, work.d4  # incumbent last terms

    # Find the most extreme abscissae corresponding with terms that are
    # included in the sum in this level.
    invalid = ~np.isfinite(fj) | (work.wj == 0)  # these terms aren't in the sum
    invalid_rl = invalid.reshape(len(fj), 2, -1)
    invalid_r, invalid_l = invalid_rl[:, 0], invalid_rl[:, 1]
    xrl = work.xj.reshape(len(fj), 2, -1).copy()  # this gets modified
    xr, xl = xrl[:, 0], xrl[:, 1]
    xr[invalid_r], xl[invalid_l] = -np.inf, np.inf
    ir, il = np.argmax(xr, axis=-1), np.argmin(xl, axis=-1)

    # Determine whether the most extreme abscissae are from this level or
    # a previous level. Update the record of the corresponding weights and
    # function values.
    frl = fj.reshape(len(fj), 2, -1)
    fr, fl = frl[:, 0], frl[:, 1]
    wrl = work.wj.reshape(len(work.wj), 2, -1)
    wr, wl = wrl[:, 0], wrl[:, 1]
    j = np.take_along_axis(xr, ir[:, np.newaxis], axis=-1).squeeze() > xr0  # xr[:, ir] > xr0
    xr0[j] = np.take_along_axis(xr[j], ir[j, np.newaxis], axis=-1).squeeze()
    fr0[j] = np.take_along_axis(fr[j], ir[j, np.newaxis], axis=-1).squeeze()
    wr0[j] = np.take_along_axis(wr[j], ir[j, np.newaxis], axis=-1).squeeze()
    j = np.take_along_axis(xl, il[:, np.newaxis], axis=-1).squeeze() < xl0  # xl[:, il] < xl0
    xl0[j] = np.take_along_axis(xl[j], il[j, np.newaxis], axis=-1).squeeze()
    fl0[j] = np.take_along_axis(fl[j], il[j, np.newaxis], axis=-1).squeeze()
    wl0[j] = np.take_along_axis(wl[j], il[j, np.newaxis], axis=-1).squeeze()

    # Compute the error estimate `d4` - the magnitude of the leftmost or
    # rightmost term, whichever is greater.
    flwl0 = fl0 + np.log(wl0) if work.log else fl0 * wl0  # leftmost term
    frwr0 = fr0 + np.log(wr0) if work.log else fr0 * wr0  # rightmost term
    magnitude = np.real if work.log else np.abs
    d4 = np.maximum(magnitude(flwl0), magnitude(frwr0))
    last_terms = xl0, fl0, wl0, xr0, fr0, wr0, d4

    # There are two approaches to dealing with function values that are
    # numerically infinite due to approaching a singularity - zero them, or
    # replace them with the function value at the nearest non-infinite point.
    # [3] pg. 22 suggests the latter, so let's do that given that we have the
    # information.
    fr0b = np.broadcast_to(fr0[:, np.newaxis], fr.shape)
    fl0b = np.broadcast_to(fl0[:, np.newaxis], fl.shape)
    fr[invalid_r] = fr0b[invalid_r]
    fl[invalid_l] = fl0b[invalid_l]

    # When wj is zero, log emits a warning
    with np.errstate(divide='ignore'):
        fjwj = fj + np.log(work.wj) if work.log else fj * work.wj

    # update integral estimate
    Sn = (special.logsumexp(fjwj + np.log(work.h), axis=-1) if work.log
          else np.sum(fjwj, axis=-1) * work.h)

    work.xl0, work.fl0, work.wr0, work.xr0, work.fr0, work.wl0, work.d4 = xl0, fl0, wr0, xr0, fr0, wl0, d4

    return fjwj, Sn


def _estimate_error(work):
    n = work.n
    Sn = work.Sn
    Sk = work.Sk
    h = work.h
    log = work.log
    last_terms = None
    # Estimate the error according to [1] Section 5

    if work.n == 0 or work.nit == 0:
        # The paper says to use "one" as the error before it can be calculated.
        # NaN seems to be more appropriate.
        nan = np.full_like(Sn, np.nan)
        return nan, nan, Sk

    fjwj = work.fjwj

    indices = _pair_cache.indices
    fjwj_rl = fjwj.reshape(len(work.Sn), 2, -1)

    # With a jump start (starting at level higher than 0), we haven't
    # explicitly calculated the integral estimate at lower levels. But we have
    # all the function value-weight products, so we can compute the
    # lower-level estimates.
    if Sk.shape[-1] == 0:
        hm1 = 2 * work.h
        fjwjm1_rl = fjwj_rl[..., :indices[n]]
        Snm1 = (special.logsumexp(fjwjm1_rl + np.log(hm1), axis=(-1, -2)) if log
                else np.sum(fjwjm1_rl, axis=(-1, -2)) * hm1)
        Sk = np.concatenate((Snm1[:, np.newaxis], Sk), axis=-1)

    if n == 1:
        nan = np.full_like(Sn, np.nan)
        return nan, nan, Sk

    # The paper says not to calculate the error for n<=2, but it's not clear
    # about whether it starts at level 0 or level 1. We start at level 0, so
    # why not compute the error beginning in level 2?
    if Sk.shape[-1] < 2:
        hm2 = 4 * work.h
        fjwjm2_rl = fjwj_rl[..., :indices[n-1]]
        Snm2 = (special.logsumexp(fjwjm2_rl + np.log(hm2), axis=(-1, -2)) if log
                else np.sum(fjwjm2_rl, axis=(-1, -2)) * hm2)
        Sk = np.concatenate((Snm2[:, np.newaxis], Sk), axis=-1)

    Snm2 = Sk[..., -2]
    Snm1 = Sk[..., -1]

    e1 = np.finfo(np.float64).eps

    if log:
        log_e1 = np.log(e1)
        # Currently, only real integrals are supported in log-scale. All
        # complex values have imaginary part in increments of pi*j, which just
        # carries sign information of the original integral, so use of
        # `np.real` here is equivalent to absolute value in real scale.
        fjwj = np.real(fjwj)
        d1 = np.real(special.logsumexp([Sn, Snm1 + np.pi*1j], axis=0))
        d2 = np.real(special.logsumexp([Sn, Snm2 + np.pi*1j], axis=0))
        d3 = log_e1 + np.max(fjwj, axis=-1)
        d4 = work.d4
        aerr = np.max([d1 ** 2 / d2, 2 * d1, d3, d4], axis=0)
        rerr = np.maximum(log_e1, aerr - np.real(Sn))
    else:
        # Note: explicit computation of log10 of each of these is unnecessary.
        fjwj = np.abs(fjwj)
        d1 = np.abs(Sn - Snm1)
        d2 = np.abs(Sn - Snm2)
        d3 = e1 * np.max(fjwj, axis=-1)
        d4 = work.d4
        # If `d1` is 0, no need to warn. This does the right thing.
        with np.errstate(divide='ignore'):
            aerr = np.max([d1**(np.log(d1)/np.log(d2)), d1**2, d3, d4], axis=0)
        rerr = np.maximum(e1, aerr/np.abs(Sn))
    return rerr, aerr.reshape(Sn.shape), Sk


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

    # these need input validation
    a, b = np.broadcast_arrays(a, b)

    params = np.asarray([atol, rtol, 0.])
    if log:
        message = '`atol` and `rtol` may not be positive infinity.'
        if np.any(np.isposinf(params[-3:-1])):
            raise ValueError(message)
    else:
        message = '`atol` and `rtol` must be non-negative and finite.'
        if np.any(params[-3:-1] < 0) or np.any(np.isinf(params[-3:-1])):
            raise ValueError(message)
    atol, rtol, _ = params

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


def _tanhsinh2(f, a, b, *, log=False, maxfun=None, maxlevel=None,
               minlevel=2, atol=None, rtol=None, args=(), callback=None):

    # need IV for callback, args
    res = _tanhsinh_iv(f, a, b, log, maxfun, maxlevel, minlevel, atol, rtol)
    (f, a, b, log, maxfun, maxlevel, minlevel, atol, rtol) = res

    # Transform improper integrals
    # func, a, b, feval_factor = _transform_integrals(f, a, b, log)
    func = f

    # Initialization
    # func(a/b) is not used. Can we get avoid doing this?
    xs, fs, args, shape, dtype = _scalar_optimization_initialize(func, (a, b), args)
    a, b = xs
    zero = -np.inf if log else 0
    Sn = np.full(shape, zero, dtype=dtype).ravel()
    Sk = np.empty_like(Sn)[:, np.newaxis][:, 0:0]  # add zero length new axis
    aerr = np.full(shape, np.nan, dtype=dtype).ravel()
    status = np.full(shape, 1, dtype=int).ravel()  # in progress - STANDARDIZE
    h0 = _get_base_step(dtype=dtype)
    maxiter = maxlevel - minlevel + 1

    xl0 = np.full(shape, np.inf, dtype=dtype).ravel()
    fl0 = np.full(shape, np.nan, dtype=dtype).ravel()
    wl0 = np.zeros(shape, dtype=dtype).ravel()
    xr0 = np.full(shape, -np.inf, dtype=dtype).ravel()
    fr0 = np.full(shape, np.nan, dtype=dtype).ravel()
    wr0 = np.zeros(shape, dtype=dtype).ravel()
    d4 = np.zeros(shape, dtype=dtype).ravel()

    nit, nfev = 0, 2  # one function evaluations performed above

    work = OptimizeResult(Sn=Sn, Sk=Sk, aerr=aerr, h0=h0, h=h0,
                          atol=atol, rtol=rtol, nit=nit, nfev=nfev,
                          status=status, dtype=dtype, minlevel=minlevel,
                          a=a[:, np.newaxis], b=b[:, np.newaxis], log=log,
                          n = minlevel, xl0=xl0, fl0=fl0, wl0=wl0, xr0=xr0,
                          fr0=fr0, wr0=wr0, d4=d4)
    # This is the correspondence between terms in the `work` object and the
    # final result. In this case, the mapping is trivial. Note that `success`
    # is prepanded automatically.
    res_work_pairs = [('status', 'status'), ('integral', 'Sn'),
                      ('error', 'aerr'), ('nit', 'nit'), ('nfev', 'nfev')]

    def pre_func_eval(work):
        work.h = work.h0 / 2**work.n
        xjc, wj = _get_pairs(work.n, work.h0, inclusive=(work.n == work.minlevel))
        work.xj, work.wj = _transform_to_limits(xjc, wj, work.a, work.b)
        return work.xj

    def post_func_eval(x, fj, work):
        fjwj, Sn = _euler_maclaurin_sum(fj, work)
        if work.Sk.shape[-1]:
            Snm1 = work.Sk[:, -1]
            Sn = (special.logsumexp([Snm1 - np.log(2), Sn], axis=0) if log
                  else Snm1 / 2 + Sn)
        work.fjwj = fjwj
        work.Sn = Sn

    def check_termination(work):
        """Terminate due to convergence, non-finite values, or error increase"""
        stop = np.zeros_like(work.Sn).astype(bool)

        if work.nit == 0:
            # The only way we can terminate on the zeroth iteration is if
            # the integration limits are equal.
            i = (work.a == work.b).ravel()  # these are guaranteed to be 1d
            zero = -np.inf if log else 0
            work.Sn[i] = zero
            work.aerr[i] = zero
            work.status[i] = 0
            stop[i] = True
            return stop

        work.rerr, work.aerr, work.Sk = _estimate_error(work)
        i = ((work.rerr < work.rtol) | (work.rerr + np.real(work.Sn) < work.atol) if log
             else (work.rerr < work.rtol) | (work.rerr * abs(work.Sn) < work.atol))
        work.status[i] = 0
        stop[i] = True

        i = ~np.isfinite(work.Sn) & ~stop
        work.status[i] = 3
        stop[i] = True

        return stop

    def post_termination_check(work):
        work.n += 1
        work.Sk = np.concatenate((work.Sk, work.Sn[:, np.newaxis]), axis=-1)
        return

    def customize_result(res):
        return

    return _scalar_optimization_loop(work, callback, shape,
                                     maxiter, func, args, dtype,
                                     pre_func_eval, post_func_eval,
                                     check_termination, post_termination_check,
                                     customize_result, res_work_pairs)