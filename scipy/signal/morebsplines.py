import numpy as np


__all__ = [ 'cubic', 'quadratic','cspline1d', 'qspline1d', 
            'cspline1d_eval', 'qspline1d_eval']


#pythran export cubic(int[])
#pythran export cubic(float[])
def cubic(x):
    """A cubic B-spline.

    This is a special case of `bspline`, and equivalent to ``bspline(x, 3)``.

    Parameters
    ----------
    x : array_like
        a knot vector

    Returns
    -------
    res : ndarray
        Cubic B-spline basis function values

    See Also
    --------
    bspline : B-spline basis function of order n
    quadratic : A quadratic B-spline.

    Examples
    --------
    We can calculate B-Spline basis function of several orders:

    >>> from scipy.signal import bspline, cubic, quadratic
    >>> bspline(0.0, 1)
    1

    >>> knots = [-1.0, 0.0, -1.0]
    >>> bspline(knots, 2)
    array([0.125, 0.75, 0.125])

    >>> np.array_equal(bspline(knots, 2), quadratic(knots))
    True

    >>> np.array_equal(bspline(knots, 3), cubic(knots))
    True

    """
    ax = abs(np.asarray(x))
    res = np.zeros_like(ax)
    cond1 = np.less(ax, 1)
    if cond1.any():
        ax1 = ax[cond1]
        res[cond1] = 2.0 / 3 - 1.0 / 2 * ax1 ** 2 * (2 - ax1)
    cond2 = ~cond1 & np.less(ax, 2)
    if cond2.any():
        ax2 = ax[cond2]
        res[cond2] = 1.0 / 6 * (2 - ax2) ** 3
    return res


#pythran export quadratic(int[])
#pythran export quadratic(float[])
def quadratic(x):
    """A quadratic B-spline.

    This is a special case of `bspline`, and equivalent to ``bspline(x, 2)``.

    Parameters
    ----------
    x : array_like
        a knot vector

    Returns
    -------
    res : ndarray
        Quadratic B-spline basis function values

    See Also
    --------
    bspline : B-spline basis function of order n
    cubic : A cubic B-spline.

    Examples
    --------
    We can calculate B-Spline basis function of several orders:

    >>> from scipy.signal import bspline, cubic, quadratic
    >>> bspline(0.0, 1)
    1

    >>> knots = [-1.0, 0.0, -1.0]
    >>> bspline(knots, 2)
    array([0.125, 0.75, 0.125])

    >>> np.array_equal(bspline(knots, 2), quadratic(knots))
    True

    >>> np.array_equal(bspline(knots, 3), cubic(knots))
    True

    """
    ax = np.abs(np.asarray(x))
    res = np.zeros_like(ax)
    cond1 = np.less(ax, 0.5)
    if cond1.any():
        ax1 = ax[cond1]
        res[cond1] = 0.75 - ax1 ** 2
    cond2 = ~cond1 & np.less(ax, 1.5)
    if cond2.any():
        ax2 = ax[cond2]
        res[cond2] = (ax2 - 1.5) ** 2 / 2.0
    return res


#pythran export _coeff_smooth(float)
def _coeff_smooth(lam):
    xi = 1 - 96 * lam + 24 * lam * np.sqrt(3 + 144 * lam)
    omeg = np.arctan2(np.sqrt(144 * lam - 1), np.sqrt(xi))
    rho = (24 * lam - 1 - np.sqrt(xi)) / (24 * lam)
    rho = rho * np.sqrt((48 * lam + 24 * lam * np.sqrt(3 + 144 * lam)) / xi)
    return rho, omeg


#pythran export _hc(int, float, float, float)
#pythran export _hc(float, float, float, float)
def _hc(k, cs, rho, omega):
    return (cs / np.sin(omega) * (rho ** k) * np.sin(omega * (k + 1)) *
            np.greater(k, -1))


#pythran export _hs(int[], float, float, float)
#pythran export _hs(float[], float, float, float)
def _hs(k, cs, rho, omega):
    c0 = (cs * cs * (1 + rho * rho) / (1 - rho * rho) /
          (1 - 2 * rho * rho * np.cos(2 * omega) + rho ** 4))
    gamma = (1 - rho * rho) / (1 + rho * rho) / np.tan(omega)
    ak = abs(k)
    return c0 * rho ** ak * (np.cos(omega * ak) + gamma * np.sin(omega * ak))


#pythran export _cubic_smooth_coeff(int[], float)
#pythran export _cubic_smooth_coeff(float[], float)
def _cubic_smooth_coeff(signal, lamb):
    rho, omega = _coeff_smooth(lamb)
    cs = 1 - 2 * rho * np.cos(omega) + rho * rho
    K = len(signal)
    yp = np.zeros((K,), signal.dtype)
    k = np.arange(K)
    yp[0] = (_hc(0, cs, rho, omega) * signal[0] +
             np.add.reduce(_hc(k + 1, cs, rho, omega) * signal))

    yp[1] = (_hc(0, cs, rho, omega) * signal[0] +
             _hc(1, cs, rho, omega) * signal[1] +
             np.add.reduce(_hc(k + 2, cs, rho, omega) * signal))

    for n in range(2, K):
        yp[n] = (cs * signal[n] + 2 * rho * np.cos(omega) * yp[n - 1] -
                 rho * rho * yp[n - 2])

    y = np.zeros((K,), signal.dtype)

    y[K - 1] = np.add.reduce((_hs(k, cs, rho, omega) +
                           _hs(k + 1, cs, rho, omega)) * signal[::-1])
    y[K - 2] = np.add.reduce((_hs(k - 1, cs, rho, omega) +
                           _hs(k + 2, cs, rho, omega)) * signal[::-1])

    for n in range(K - 3, -1, -1):
        y[n] = (cs * yp[n] + 2 * rho * np.cos(omega) * y[n + 1] -
                rho * rho * y[n + 2])

    return y


#pythran export _cubic_coeff(int[])
#pythran export _cubic_coeff(float[])
def _cubic_coeff(signal):
    zi = -2 + np.sqrt(3)
    K = len(signal)
    yplus = np.zeros((K,), signal.dtype)
    powers = zi ** np.arange(K)
    yplus[0] = signal[0] + zi * np.add.reduce(powers * signal)
    for k in range(1, K):
        yplus[k] = signal[k] + zi * yplus[k - 1]
    output = np.zeros((K,), signal.dtype)
    output[K - 1] = zi / (zi - 1) * yplus[K - 1]
    for k in range(K - 2, -1, -1):
        output[k] = zi * (output[k + 1] - yplus[k])
    return output * 6.0


#pythran export _quadratic_coeff(int[])
#pythran export _quadratic_coeff(float[])
def _quadratic_coeff(signal):
    zi = -3 + 2 * np.sqrt(2.0)
    K = len(signal)
    yplus = np.zeros((K,), signal.dtype)
    powers = zi ** np.arange(K)
    yplus[0] = signal[0] + zi * np.add.reduce(powers * signal)
    for k in range(1, K):
        yplus[k] = signal[k] + zi * yplus[k - 1]
    output = np.zeros((K,), signal.dtype)
    output[K - 1] = zi / (zi - 1) * yplus[K - 1]
    for k in range(K - 2, -1, -1):
        output[k] = zi * (output[k + 1] - yplus[k])
    return output * 8.0


#pythran export cspline1d(int[])
#pythran export cspline1d(float[])
#pythran export cspline1d(int[], float)
#pythran export cspline1d(float[], float)
def cspline1d(signal, lamb=0.0):
    """
    Compute cubic spline coefficients for rank-1 array.

    Find the cubic spline coefficients for a 1-D signal assuming
    mirror-symmetric boundary conditions. To obtain the signal back from the
    spline representation mirror-symmetric-convolve these coefficients with a
    length 3 FIR window [1.0, 4.0, 1.0]/ 6.0 .

    Parameters
    ----------
    signal : ndarray
        A rank-1 array representing samples of a signal.
    lamb : float, optional
        Smoothing coefficient, default is 0.0.

    Returns
    -------
    c : ndarray
        Cubic spline coefficients.

    See Also
    --------
    cspline1d_eval : Evaluate a cubic spline at the new set of points.

    Examples
    --------
    We can filter a signal to reduce and smooth out high-frequency noise with
    a cubic spline:

    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import cspline1d, cspline1d_eval
    >>> rng = np.random.default_rng()
    >>> sig = np.repeat([0., 1., 0.], 100)
    >>> sig += rng.standard_normal(len(sig))*0.05  # add noise
    >>> time = np.linspace(0, len(sig))
    >>> filtered = cspline1d_eval(cspline1d(sig), time)
    >>> plt.plot(sig, label="signal")
    >>> plt.plot(time, filtered, label="filtered")
    >>> plt.legend()
    >>> plt.show()

    """
    if lamb != 0.0:
        return _cubic_smooth_coeff(signal, lamb)
    else:
        return _cubic_coeff(signal)


#pythran export qspline1d(int[])
#pythran export qspline1d(float[])
#pythran export qspline1d(int[], float)
#pythran export qspline1d(float[], float)
def qspline1d(signal, lamb=0.0):
    """Compute quadratic spline coefficients for rank-1 array.

    Parameters
    ----------
    signal : ndarray
        A rank-1 array representing samples of a signal.
    lamb : float, optional
        Smoothing coefficient (must be zero for now).

    Returns
    -------
    c : ndarray
        Quadratic spline coefficients.

    See Also
    --------
    qspline1d_eval : Evaluate a quadratic spline at the new set of points.

    Notes
    -----
    Find the quadratic spline coefficients for a 1-D signal assuming
    mirror-symmetric boundary conditions. To obtain the signal back from the
    spline representation mirror-symmetric-convolve these coefficients with a
    length 3 FIR window [1.0, 6.0, 1.0]/ 8.0 .

    Examples
    --------
    We can filter a signal to reduce and smooth out high-frequency noise with
    a quadratic spline:

    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import qspline1d, qspline1d_eval
    >>> rng = np.random.default_rng()
    >>> sig = np.repeat([0., 1., 0.], 100)
    >>> sig += rng.standard_normal(len(sig))*0.05  # add noise
    >>> time = np.linspace(0, len(sig))
    >>> filtered = qspline1d_eval(qspline1d(sig), time)
    >>> plt.plot(sig, label="signal")
    >>> plt.plot(time, filtered, label="filtered")
    >>> plt.legend()
    >>> plt.show()

    """
    if lamb != 0.0:
        raise ValueError("Smoothing quadratic splines not supported yet.")
    else:
        return _quadratic_coeff(signal)


#pythran export cspline1d_eval(int[], int[])
#pythran export cspline1d_eval(int[], int[], float)
#pythran export cspline1d_eval(int[], int[], float, int)
#pythran export cspline1d_eval(int[], float[])
#pythran export cspline1d_eval(int[], float[], float)
#pythran export cspline1d_eval(int[], float[], float, int)
#pythran export cspline1d_eval(float[], int[])
#pythran export cspline1d_eval(float[], int[], float)
#pythran export cspline1d_eval(float[], int[], float, int)
#pythran export cspline1d_eval(float[], float[])
#pythran export cspline1d_eval(float[], float[], float)
#pythran export cspline1d_eval(float[], float[], float, int)
def cspline1d_eval(cj, newx, dx=1.0, x0=0):
    """Evaluate a cubic spline at the new set of points.

    `dx` is the old sample-spacing while `x0` was the old origin. In
    other-words the old-sample points (knot-points) for which the `cj`
    represent spline coefficients were at equally-spaced points of:

      oldx = x0 + j*dx  j=0...N-1, with N=len(cj)

    Edges are handled using mirror-symmetric boundary conditions.

    Parameters
    ----------
    cj : ndarray
        cublic spline coefficients
    newx : ndarray
        New set of points.
    dx : float, optional
        Old sample-spacing, the default value is 1.0.
    x0 : int, optional
        Old origin, the default value is 0.

    Returns
    -------
    res : ndarray
        Evaluated a cubic spline points.

    See Also
    --------
    cspline1d : Compute cubic spline coefficients for rank-1 array.

    Examples
    --------
    We can filter a signal to reduce and smooth out high-frequency noise with
    a cubic spline:

    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import cspline1d, cspline1d_eval
    >>> rng = np.random.default_rng()
    >>> sig = np.repeat([0., 1., 0.], 100)
    >>> sig += rng.standard_normal(len(sig))*0.05  # add noise
    >>> time = np.linspace(0, len(sig))
    >>> filtered = cspline1d_eval(cspline1d(sig), time)
    >>> plt.plot(sig, label="signal")
    >>> plt.plot(time, filtered, label="filtered")
    >>> plt.legend()
    >>> plt.show()

    """
    newx = (np.asarray(newx) - x0) / float(dx)
    res = np.zeros_like(newx, dtype=cj.dtype)
    if res.size == 0:
        return res
    N = len(cj)
    cond1 = newx < 0
    cond2 = newx > (N - 1)
    cond3 = ~(cond1 | cond2)
    # handle general mirror-symmetry
    res[cond1] = cspline1d_eval(cj, -newx[cond1])
    res[cond2] = cspline1d_eval(cj, 2 * (N - 1) - newx[cond2])
    newx = newx[cond3]
    if newx.size == 0:
        return res
    result = np.zeros_like(newx, dtype=cj.dtype)
    jlower = np.floor(newx - 2).astype(int) + 1
    for i in range(4):
        thisj = jlower + i
        indj = thisj.clip(0, N - 1)  # handle edge cases
        result += cj[indj] * cubic(newx - thisj)
    res[cond3] = result
    return res


#pythran export qspline1d_eval(int[], int[])
#pythran export qspline1d_eval(int[], int[], float)
#pythran export qspline1d_eval(int[], int[], float, int)
#pythran export qspline1d_eval(int[], float[])
#pythran export qspline1d_eval(int[], float[], float)
#pythran export qspline1d_eval(int[], float[], float, int)
#pythran export qspline1d_eval(float[], int[])
#pythran export qspline1d_eval(float[], int[], float)
#pythran export qspline1d_eval(float[], int[], float, int)
#pythran export qspline1d_eval(float[], float[])
#pythran export qspline1d_eval(float[], float[], float)
#pythran export qspline1d_eval(float[], float[], float, int)
def qspline1d_eval(cj, newx, dx=1.0, x0=0):
    """Evaluate a quadratic spline at the new set of points.

    Parameters
    ----------
    cj : ndarray
        Quadratic spline coefficients
    newx : ndarray
        New set of points.
    dx : float, optional
        Old sample-spacing, the default value is 1.0.
    x0 : int, optional
        Old origin, the default value is 0.

    Returns
    -------
    res : ndarray
        Evaluated a quadratic spline points.

    See Also
    --------
    qspline1d : Compute quadratic spline coefficients for rank-1 array.

    Notes
    -----
    `dx` is the old sample-spacing while `x0` was the old origin. In
    other-words the old-sample points (knot-points) for which the `cj`
    represent spline coefficients were at equally-spaced points of::

      oldx = x0 + j*dx  j=0...N-1, with N=len(cj)

    Edges are handled using mirror-symmetric boundary conditions.

    Examples
    --------
    We can filter a signal to reduce and smooth out high-frequency noise with
    a quadratic spline:

    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import qspline1d, qspline1d_eval
    >>> rng = np.random.default_rng()
    >>> sig = np.repeat([0., 1., 0.], 100)
    >>> sig += rng.standard_normal(len(sig))*0.05  # add noise
    >>> time = np.linspace(0, len(sig))
    >>> filtered = qspline1d_eval(qspline1d(sig), time)
    >>> plt.plot(sig, label="signal")
    >>> plt.plot(time, filtered, label="filtered")
    >>> plt.legend()
    >>> plt.show()

    """
    newx = (np.asarray(newx) - x0) / dx
    res = np.zeros_like(newx)
    if res.size == 0:
        return res
    N = len(cj)
    cond1 = newx < 0
    cond2 = newx > (N - 1)
    cond3 = ~(cond1 | cond2)
    # handle general mirror-symmetry
    res[cond1] = qspline1d_eval(cj, -newx[cond1])
    res[cond2] = qspline1d_eval(cj, 2 * (N - 1) - newx[cond2])
    newx = newx[cond3]
    if newx.size == 0:
        return res
    result = np.zeros_like(newx)
    jlower = np.floor(newx - 1.5).astype(int) + 1
    for i in range(3):
        thisj = jlower + i
        indj = thisj.clip(0, N - 1)  # handle edge cases
        result += cj[indj] * quadratic(newx - thisj)
    res[cond3] = result
    return res
