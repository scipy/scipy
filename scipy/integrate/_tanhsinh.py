from dataclasses import dataclass
import numpy as np

# todo:
#  documentation
#  cache abscissa/weight pairs
#  tests - test rtol, maxfun, and minweight
#  respect data types
#  callback
#  log integration
#  warn (somehow) when invalid function values & weight < minweight
#  vectorize

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

def _compute_pairs(k):
    # Compute the abscissa-weight pairs for each level m. See [1] page 9.

    "....each level k of abscissa-weight pairs uses h = 2 **-k"
    h = 1 / (2 ** k)

    # "We find that roughly 3.6 * 2^k abscissa-weight pairs are generated at
    # "level k." The actual number per level can be generated like:
    # for i in range(10):
    #     _, xjc, wj = _compute_pairs(i)
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
    return h, xjc, wj

def _tanhsinh(f, a, b, *, maxfun=5000, maxiter=10, atol=0, rtol=1e-14,
              minweight=1e-100):

    # Input validation and standardization
    res = _quadts_iv(f, a, b, maxfun, maxiter, atol, rtol, minweight)
    f, a, b, maxfun, maxiter, atol, rtol, minweight, feval_factor = res

    # Initialization
    Sk = []  # sequence of integral estimates for error estimation
    Sn = aerr = np.nan  # integral and error are NaN until determined otherwise
    status = -1  # "Iteration in progress."
    feval = 0  # function evaluation counter

    for n in range(maxiter):
        h, xjc, wj = _compute_pairs(n)

        # Transform integral according to user-specified limits. This is just
        # math that follows from the fact that the standard limits are (-1, 1).
        # Note: If we had stored xj instead of xjc, we would have
        # xj = alpha * xj + beta, where beta = (a + b)/2
        alpha = (b-a)/2
        xj = np.concatenate((-alpha * xjc + b, alpha * xjc + a))
        wj *= alpha
        wj = np.concatenate((wj, wj))

        if feval + len(xj) * feval_factor > maxfun:
            status = 1
            break
        with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
            fj = f(xj)
        feval += len(xj) * feval_factor
        fjwj = fj * wj  # I'd use @, but we'll need the individual products

        # Points on the boundaries can be generated due to finite precision
        # arithmetic. Ideally we wouldn't evaluate the function at these points
        # or when weights are zero; however, it would be tricky to filter out
        # points when this function is vectorized. Set the results to zero.
        invalid = (xj <= a) | (xj >= b) | (wj <= minweight)
        fjwj[invalid] = 0

        # Check for infinities
        if not np.all(np.isfinite(fjwj)):
            status = 3
            break

        # update integral estimate
        Snm1 = 0 if not Sk else Sk[-1]
        Sn = Snm1/2 + np.sum(fjwj) * h

        # Check error estimate (see "5. Error Estimation, page 11")
        if n >= 2:
            Snm2, Snm1 = Sk[-2:]
            with np.errstate(divide='ignore'):  # when values are zero
                d1 = np.log10(abs(Sn - Snm1))
                d2 = np.log10(abs(Sn - Snm2))
                e1 = np.finfo(np.float64).eps
                d3 = np.log10(e1 * np.max(np.abs(fjwj)))
                d4 = np.log10(np.max(np.reshape(np.abs(fjwj), (2, -1))[:, -1]))
            d = np.max([d1**2/d2, 2*d1, d3, d4])
            rerr = 10**d
            aerr = max(np.finfo(np.float64).eps, rerr) * abs(Sn)
            if rerr < rtol or rerr*abs(Sn) < atol:
                status = 0
                break

        Sk.append(Sn)
    else:
        status = 2

    message = _status_messages[status]
    return QuadratureResult(integral=Sn, error=aerr, feval=feval,
                            success=status==0, status=status, message=message)


def _quadts_iv(f, a, b, maxfun, maxiter, atol, rtol, minweight):

    message = '`f` must be callable.'
    if not callable(f):
        raise ValueError(message)

    message = ('Integration limits `a` and `b`, tolerances `atol` and '
               '`rtol`, and `minweight` must be reals.')
    params = np.asarray([a, b, atol, rtol, minweight])
    if not np.issubdtype(params.dtype, np.floating) or np.any(np.isnan(params)):
        raise ValueError(message)
    message = '`rtol` and `minweight` must be positive and finite.'
    if np.any(params[-2:] <= 0) or  np.any(np.isinf(params[-2:])):
        raise ValueError(message)
    a, b, atol, rtol, minweight = params
    message = '`atol` must be non-negative and finite.'
    if atol < 0 or np.isinf(atol):
        raise ValueError(message)

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
            return -f(x)
        a, b = b, a

    if np.isinf(a) and np.isinf(b):
        def f(x, f=f):
            return f(x) + f(-x)
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
            return f(1/x - 1 + a)*x**-2
        a, b = 0, 1

    return f, a, b, maxfun, maxiter, atol, rtol, minweight, feval_factor
