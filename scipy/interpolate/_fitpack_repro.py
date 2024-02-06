""" Replicate FITPACK's logic, in bits and pieces.
"""


import operator
import numpy as np

from scipy.interpolate._bsplines import (
    _not_a_knot,  BSpline,
    _lsq_solve_qr
)

#    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#    c  part 1: determination of the number of knots and their position     c
#    c  **************************************************************      c
#
# https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fpcurf.f#L31

# Hardcoded in curfit.f
TOL = 0.001
MAXIT = 20


def _get_residuals(x, y, t, k, w):
 #   from scipy.interpolate._bsplines import make_lsq_spline
    # FITPACK has (w*(spl(x)-y))**2; make_lsq_spline has w*(spl(x)-y)**2
    w2 = w**2

    # inline the relevant part of >>> spl = make_lsq_spline(x, y, w=w2, t=t, k=k)
    assert y.ndim == 1
    yy = y[:, None]
    _, _, c = _lsq_solve_qr(x, yy, t, k, w)
    c = c.reshape(c.shape[0])
    c = np.ascontiguousarray(c)
    spl = BSpline.construct_fast(t, c, k)

    residuals = w2 * (spl(x) - y)**2
    return residuals


def _add_knot(x, t, k, residuals):
    """Add a new knot.

    (Approximately) replicate FITPACK's logic:
      1. split the `x` array into knot intervals, ``t(j+k) <= x(i) <= t(j+k+1)``
      2. find the interval with the maximum sum of residuals
      3. insert a new knot into the middle of that interval.

    NB: a new knot is in fact an `x` value at the middle of the interval.
    So *the knots are a subset of `x`*.

    This routine is an analog of
    https://github.com/scipy/scipy/blob/v1.11.4/scipy/interpolate/fitpack/fpcurf.f#L190-L215
    (cf _split function)

    and https://github.com/scipy/scipy/blob/v1.11.4/scipy/interpolate/fitpack/fpknot.f
    """
    fparts, ix = _split(x, t, k, residuals)

    ### redo
    assert all(ix == np.searchsorted(x, t[k:-k]))

    # find the interval with max fparts and non-zero number of x values inside
    idx_max = -101
    fpart_max = -1e100
    for i in range(len(fparts)):
        if ix[i+1] - ix[i] > 1 and fparts[i] > fpart_max:
            idx_max = i
            fpart_max = fparts[i]

    if idx_max == -101:
        raise ValueError("Internal error, please report it to SciPy developers.")

    # round up, like Dierckx does? This is really arbitrary though.
    idx_newknot = (ix[idx_max] + ix[idx_max+1] + 1) // 2
    new_knot = x[idx_newknot]

    idx_t = np.searchsorted(t, new_knot)
    t_new = np.r_[t[:idx_t], new_knot, t[idx_t:]]

    return t_new


def _split(x, t, k, residuals):
    """Split the `x` array into knot intervals and compute the residuals.

    The intervals are `t(j+k) <= x(i) <= t(j+k+1)`. Return the arrays of
    the start and end `x` indices of the intervals and the sums of residuals:
    `fparts[i]` corresponds to the interval
    `x_intervals[i] <= xvalue <= x_intervals[i+1]]`.

    This routine is a (best-effort) translation of
    https://github.com/scipy/scipy/blob/v1.11.4/scipy/interpolate/fitpack/fpcurf.f#L190-L215
    """
# c  search for knot interval t(number+k) <= x <= t(number+k+1) where
# c  fpint(number) is maximal on the condition that nrdata(number)
# c  not equals zero.
   # residuals = _get_residuals(x, y, t, k, w)

    interval = k+1
    x_intervals = [0]
    fparts = []
    fpart = 0.0
    for it in range(len(x)):
        xv, rv = x[it], residuals[it]
        fpart += rv

        if (xv >= t[interval]) and interval < len(t) - k - 1:
            # end of the current t interval: split the weight at xv by 1/2
            # between two intervals
            carry = rv / 2.0
            fpart -= carry
            fparts.append(fpart)

            fpart = carry
            interval += 1

            x_intervals.append(it)

    x_intervals.append(len(x)-1)
    fparts.append(fpart)

    return fparts, x_intervals

    '''
    The whole _split routine is basically this:

    ix = np.searchsorted(x, t[k:-k])
    # sum half-open intervals
    fparts = [residuals[ix[i]:ix[i+1]].sum() for i in range(len(ix)-1)]
    carries = residuals[x[ix[1:-1]]]

    for i in range(len(carries)):     # split residuals at internal knots
        carry = carries[i] / 2
        fparts[i] += carry
        fparts[i+1] -= carry

    fparts[-1] += residuals[-1]       # add the contribution of the last knot

    assert sum(fparts) == sum(residuals)
    '''



def _validate_inputs(x, y, w, k, s, xb, xe):
    """Common input validations for generate_knots and make_splrep.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if y.ndim != 1:
        raise NotImplementedError(f"{y.ndim = } not implemented yet.")

    if w is None:
        w = np.ones_like(y, dtype=float)
    else:
        w = np.asarray(w, dtype=float)
        if (w < 0).any():
            raise ValueError("Weights must be non-negative")
    if w.ndim != 1 or w.shape[0] != x.shape[0]:
        raise ValueError(f"Weights is incompatible: {w.shape =} != {x.shape}.")

    if x.shape[0] != y.shape[0]:
        raise ValueError(f"Data is incompatible: {x.shape = } and {y.shape = }.")
    if (x[1:] < x[:-1]).any():
        raise ValueError("Expect `x` to be an ordered 1D sequence.")

    k = operator.index(k)

    if s < 0:
        raise ValueError(f"`s` must be non-negative. Got {s = }")

    if xb is None:
        xb = min(x)
    if xe is None:
        xe = max(x)

    return x, y, w, k, s, xb, xe


def generate_knots(x, y, k=3, *, s=0, w=None, nest=None, xb=None, xe=None):
    """Replicate FITPACK's constructing the knot vector.

    Parameters
    ----------
    x, y : array_like
        The data points defining the curve ``y = f(x)``.
    k : int, optional
        The spline degree. Default is cubic, ``k = 3``.
    s : float, optional
        The smoothing factor.
    w : array_like, optional
        Weights.
    nest : int, optional
        Stop when at least this many knots are placed.
    xb : float, optional
        The boundary of the approximation interval. If None (default),
        is set to ``x[0]``.
    xe : float, optional
        The boundary of the approximation interval. If None (default),
        is set to ``x[-1]``.

    Yields
    ------
    t : ndarray
        Knot vectors with an increasing number of knots.
        The generator is finite: it stops when the smoothing critetion is
        satisfied, or when then number of knots exceeds the maximum value:
        the user-provided `nest` or `x.size + k + 1` --- which is the knot vector
        for the interpolating spline.

    Examples
    --------
    Generate some noisy data and fit a sequence of LSQ splines:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.interpolate import make_lsq_spline, generate_knots
    >>> rng = np.random.default_rng(12345)
    >>> x = np.linspace(-3, 3, 50)
    >>> y = np.exp(-x**2) + 0.1 * rng.standard_normal(size=50)

    >>> knots = list(generate_knots(x, y, s=1e-10))
    >>> for t in knots[::3]:
    ...     spl = make_lsq_spline(x, y, t)
    ...     xs = xs = np.linspace(-3, 3, 201)
    ...     plt.plot(xs, spl(xs), '-', label=f'n = {len(t)}', lw=3, alpha=0.7)
    >>> plt.plot(x, y, 'o', label='data')
    >>> plt.plot(xs, np.exp(-xs**2), '--')
    >>> plt.legend()

    Note that increasing the number of knots make the result follow the data
    more and more closely.

    Also note that a step of the generator may add multiple knots:

    >>> [len(t) for t in knots]
    [8, 9, 10, 12, 16, 24, 40, 48, 52]

    Notes
    -----
    The routine generates successive knots vectors of increasing length, starting
    from ``2*(k+1)`` to ``len(x) + k + 1``, trying to make knots more dense
    in the regions where the deviation of the LSQ spline from data is large.

    When the maximum number of knots, ``len(x) + k + 1`` is reached
    (this happens when ``s`` is small and ``nest`` is large), the generator
    stops, and the last output is the knots for interpolation with the
    not-a-knot boundary condition.

    Knots are located at data sites, unless ``k`` is even and the number of knots
    is ``len(x) + k + 1``. In that case, the last output of the generator
    has internal knots at Greville sites, ``(x[1:] + x[:-1]) / 2``.
    """
    if s == 0:
        if nest is not None or w is not None:
            raise ValueError("s == 0 is interpolation only")
        t = _not_a_knot(x, k)
        yield t
        return

    x, y, w, k, s, xb, xe = _validate_inputs(x, y, w, k, s, xb, xe)

    acc = s * TOL
    m = x.size    # the number of data points

    if nest is None:
        # the max number of knots. This is set in _fitpack_impl.py line 274
        # and fitpack.pyf line 198
        nest = max(m + k + 1, 2*k + 3)
    else:
        if nest < 2*(k+1):
            raise ValueError(f"`nest` too small: {nest = } < 2*(k+1) = {2*(k+1)}.")

    nmin = 2*(k+1)    # the number of knots for an LSQ polynomial approximation
    nmax = m + k + 1  # the number of knots for the spline interpolation

    # start from no internal knots
    t = np.asarray([xb]*(k+1) + [xe]*(k+1), dtype=float)
    n = t.shape[0]
    fp = 0.0
    fpold = 0.0

    # c  main loop for the different sets of knots. m is a safe upper bound
    # c  for the number of trials.
    for _ in range(m):
        yield t

        # construct the LSQ spline with this set of knots
        fpold = fp
        residuals = _get_residuals(x, y, t, k, w=w)
        fp = residuals.sum()
        fpms = fp - s

        # c  test whether the approximation sinf(x) is an acceptable solution.
        # c  if f(p=inf) < s accept the choice of knots.
        if (abs(fpms) < acc) or (fpms < 0):
            return

        # ### c  increase the number of knots. ###

        # c  determine the number of knots nplus we are going to add.
        if n == nmin:
            # the first iteration
            nplus = 1
        else:
            delta = fpold - fp
            npl1 = int(nplus * fpms / delta) if delta > acc else nplus*2
            nplus = min(nplus*2, max(npl1, nplus//2, 1))

        # actually add knots
        for j in range(nplus):
            t = _add_knot(x, t, k, residuals)

            # check if we have enough knots already

            n = t.shape[0]
            # c  if n = nmax, sinf(x) is an interpolating spline.
            # c  if n=nmax we locate the knots as for interpolation.
            if (n >= nmax):
                t = _not_a_knot(x, k)
                yield t
                return

            # c  if n=nest we cannot increase the number of knots because of
            # c  the storage capacity limitation.
            if (n >= nest):
                yield t
                return

            # recompute if needed
            if j < nplus - 1:
                residuals = _get_residuals(x, y, t, k, w=w)

    # this should never be reached
    return


def construct_knot_vector(x, y, *, s, k=3, w=None, nest=None, xb=None, xe=None):
    # return the last value generated
    # XXX: needed? vs list(generate_knots(...))[-1]
    for t in generate_knots(x, y, k=k, s=s, w=w, nest=nest, xb=xb, xe=xe):
        pass
    return t

