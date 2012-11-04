# Contributed by Thouis Jones, 2009.
from numpy import log, array, isfinite, inf
from quadrature import vectorize1, AccuracyWarning

__all__ = ["quad_de"]

def quad_de(func, a, b, args=(), tol=1e-10, vec_func=True, max_level=None):
    """
    Integrate a function over an interval using double exponential
    quadrature.

    This method is well-suited to integrating analytic functions.  It
    should not be used with functions with discontinuities, although
    singularities at the endpoints are tolerated.

    If quad_de() encounters a singularity (+/-infinity or NaN) within
    the integration limits, it will raise a ValueError.

    Parameters
    ----------
    func : callable
        The integrand
    a : float
        Lower integration limit.
    b : float
        Upper integration limit.
    args : tuple, optional
        Additional arguments to the function.
    tol : float, optional
        Desired absolute error tolerance.
    vec_func : bool, optional
        Whether the function is vectorized.
    max_level : int, optional
        Maximum number of subdivisions.

    Returns
    -------
    val : float
        Estimated value of the integral
    err : float
        Estimated error

    Raises
    ------
    ValueError
        If the function returns non-finite values.

    See Also
    --------
    quad

    Notes
    -----
    Double exponential quadrature (AKA double exponential integration
    and tanh-sinh quadrature [TS]) is an integration method introduced
    by Takahashi and Mori in 1974 [TM].  The method applies a change
    of variables to map the integration limits to [-inf, inf], with a
    doubly exponential falloff in either direction, then uses
    quadrature with uniform intervals to compute the integral.

    This is based on John D. Cook's implementation [FN, DE], used with
    permission.

    Examples
    --------
    >>> from scipy.integrate import quad_de
    >>> quad_de(lambda x: np.cos(x)/(1 + x**4), 0, np.inf)
    (0.77218652541156108, 7.0590915539313116e-12)

    Weak singularity:

    >>> def f(x):
    ...     return 1/((1 - x)**.5 * (1 + x)**.25)
    >>> quad_de(f, -1, 1)
    (2.8496737568265096, 1.8710389421139337e-09)

    Strong singularity:

    >>> def f(x):
    ...     return 1/((1 - x)**.5 * (1 + x)**.9)
    >>> quad_de(f, -1, 1)
    (8.3860477551851123, 0.0030324420851587171)

    Note that in this case, `quad_de` fails to achieve the desired
    accuracy.  The error estimate is also only an estimate -- the
    exact answer here is 8.581295... and the actual error is larger
    than estimated. In these cases, you can use `quad` which may get
    better results.

    References
    ----------
    .. [DE] http://www.johndcook.com/double_exponential_integration.html
    .. [FN] http://www.codeproject.com/KB/recipes/FastNumericalIntegration.aspx
    .. [TH] http://en.wikipedia.org/wiki/Tanh-sinh_quadrature
    .. [TM] Double exponential formulas for numerical integration, 
            Hidetosi Takahasi and Masatake Mori 
            Publ. Res. Inst. Math. Sci. Volume 9, Number 3 (1973), 721-741.  
            http://dx.doi.org/10.2977%2Fprims%2F1195192451 
    """

    # vectorize func, and fold in args
    vfunc = vectorize1(func, args, vec_func=vec_func)
    
    # linear change of variables from (a,b) to (-1, 1), (0, inf) or (-inf, inf)
    # Int[f(x), {x,a,b}] = c*Int[f(ct+d), {t,-1,1}]

    if isfinite(a) and isfinite(b):
        c = float(b - a) / 2
        d = float(a + b) / 2
        _weights = _weights_m1_1
        _abscissas = _abscissas_m1_1
        _initial_step_size = _initial_step_size_m1_1
    else:
        _initial_step_size = 1
        if b < a:
            _initial_step_size *= -1
            a, b = b, a

        if a == -inf and b == inf:
            c = 1
            d = 0
            _weights = _weights_minf_inf
            _abscissas = _abscissas_minf_inf
            _initial_step_size *= _initial_step_size_minf_inf
        else:
            _weights = _weights_0_inf
            _abscissas = _abscissas_0_inf
            _initial_step_size *= _initial_step_size_0_inf
            if a == -inf and isfinite(b):
                c = -1
                d = b
            elif isfinite(a) and b == inf:
                c = 1
                d = a
            else:
                raise ValueError("Invalid bounds (%g, %g)" % (a, b))

    def contribution_at_level(l):
        # _abscissas and _weights are defined below
        return (2.0**(-l) * _weights[l] * vfunc(c * _abscissas[l] + d)).sum()

    def raise_bad_value(l):
        vals = vfunc(c * _abscissas[l] + d)
        raise ValueError('Singularities encountered in quad_de evaluating %s([%s]) = [%s]', 
                         func.__name__, 
                         ", ".join(['%f'%v for v in (c * _abscissas[l] + d)[~isfinite(vals)]]),
                         ", ".join(['%f'%v for v in vals[~isfinite(vals)]]))

    integral = contribution_at_level(0)
    if not isfinite(integral):
        raise_bad_value(0)
    current_delta = inf

    if max_level is None:
        max_level = inf

    error_estimate = inf

    for level in range(1, min(len(_abscissas), max_level)):
        new_contribution = contribution_at_level(level)

        previous_delta = current_delta
        current_delta = abs(integral / 2 - new_contribution)

        integral = integral / 2 + new_contribution

        if not isfinite(integral):
            raise_bad_value(level)

        if current_delta == 0.0:
            error_estimate = 0.0
            break

        if level == 1:
            continue

        # both of these are nonzero, or we would have already exited.
        r = log(current_delta) / log(previous_delta)
        
        if 1.9 < r < 2.1:
            # If convergence theory applied perfectly, r would be 2 in
            # the convergence region.  r close to 2 is good enough. We
            # expect the difference between this integral estimate and
            # the next one to be roughly delta^2.
            error_estimate = c * current_delta ** 2
        else:
            # not yet converging
            error_estimate = c * current_delta * 2

        if error_estimate < 0.1 * tol:
            break

    if error_estimate > tol:
        import warnings
        warnings.warn("Requested tolerance not satisfied", AccuracyWarning)

    value = c * integral * _initial_step_size
    error = abs(_initial_step_size * error_estimate)
    return value, error

# These are computed from _generate_doubleexp.abw_*() for the
# transformed range [-z, z] (where the abscissa `z` is chosen
# according to numerical limits).
#
# This points out one way that the integrator's accuracy might be
# improved: accept two versions of f(), one for points near the lower
# boundary, one for points nearer the upper boundary, with the point
# expressed as the distance from that boundary.  Another improvement
# would involve passing in the weight to be incorporated directly into
# the computation of f().  For functions that are more easily computed
# in the log domain, this could help preserve accuracy.

from _doubleexp_coefs import \
     _initial_step_size_m1_1, \
     _abscissas_raw_m1_1, \
     _weights_raw_m1_1, \
     _initial_step_size_0_inf, \
     _abscissas_raw_0_inf, \
     _weights_raw_0_inf, \
     _initial_step_size_minf_inf, \
     _abscissas_raw_minf_inf, \
     _weights_raw_minf_inf

def _mirror_abscissas(a):
    return a + [-v for v in a if v > 0]

def _mirror_weights(w, a):
    # special case for abscissa == 0, to prevent duplicate evaluations
    return w + [v for idx, v in enumerate(w) if a[idx] > 0]

_abscissas_m1_1 = [array(_mirror_abscissas(a)) for a in _abscissas_raw_m1_1]
_weights_m1_1 = [array(_mirror_weights(w, a))
                 for w, a in zip(_weights_raw_m1_1, _abscissas_raw_m1_1)]

_abscissas_0_inf = map(array, _abscissas_raw_0_inf)
_weights_0_inf = map(array, _weights_raw_0_inf)

_abscissas_minf_inf = [array(_mirror_abscissas(a)) for a in _abscissas_raw_minf_inf]
_weights_minf_inf = [array(_mirror_weights(w, a))
                     for w, a in zip(_weights_raw_minf_inf, _abscissas_raw_minf_inf)]
