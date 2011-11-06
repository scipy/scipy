
import numpy as np
from numpy import array
from scipy.interpolate import interp1d
from scipy.special import gamma


def _tukeylambda_var0(lam, order=3):
    """Taylor polynomial for the variance of the Tukey Lambda distribution.

    This function implements the Taylor polynomial (up to degree 3) at
    lambda=0 of the variance of the Tukey Lambda distribution.

    Parameters
    ----------
    lam : scalar or ndarray
        The values of lambda
    order : int (0, 1, 2, or 3)
        The order of the Taylor polynomial to compute.

    Returns
    -------
    v : scalar or ndarray
        Value of the Taylor polynomial of the variance.
    """
    if order < 0 or order > 3:
        raise ValueError("order must be at least 0 and at most 3.")

    # In the formulas in the comments, z3 = zeta(3) and z5 = zeta(5).
    # c0 = pi**2 / 3
    c0 =   3.28986813369645287
    # c1 = - (2*pi**2 + 12*z3) / 3
    c1 = -11.38796388003128288
    # c2 = (pi**2*(3*pi**2 + 80) + 480*z3) / 60
    c2 =  27.64638231176268763
    # c3 = - (pi**2*(3*pi**2 + (-20*z3 + 80)) + 480*z3 + 360*z5) / 30
    c3 = -59.82668028405662991
    c = [c0, c1, c2, c3]

    r = c[order]
    for k in range(order - 1, -1, -1):
        r = r * lam + c[k]
    return r


def tukeylambda_variance(lam):
    """Variance of the Tukey Lambda distribution.
    
    This implementation does not handle large values of lambda.  If `lam` is
    greater than 98, a NotImplementedError exception is raised.

    Parameters
    ----------
    lam : array_like
        The lambda values at which to compute the variance.

    Returns
    -------
    v : ndarray
        The variance.  For lam < -0.5, the variance is not defined, so
        np.nan is returned.  For lam = 0.5, np.inf is returned.

    Raises
    ------
    NotImplementedError if `lam` > 98.

    Notes
    -----
    In an interval around lambda=0, this function uses the degree 3 Taylor
    polynomial to compute the variance.  Otherwise it uses the standard
    formula (http://en.wikipedia.org/wiki/Tukey_lambda_distribution).  The
    Taylor polynomial is used because the standard formula has a removable
    discontinuity at lambda = 0, and does produce accurate numerical results
    near lambda = 0.  The order of the Taylor polynomial and the interval
    in which it is used were chosen to ensure an absolute error of less than
    1e-9 in a neighborhood of lambda = 0.
    """
    lam = np.asarray(lam)
    shp = lam.shape
    lam = np.atleast_1d(lam).astype(np.float64)

    if np.any(lam > 98):
        raise NotImplementedError("this function cannot compute the "
                                  "variance for lam > 98")

    # For absolute values of lam less than threshold, use the third order
    # Maclaurin series to evaluate the variance of the distribution.
    # The threshold 1.61e-3 with the order=3 Taylor polynomial gives good
    # a maximum absolute error of less than 1e-9 in the interval [-0.05, 0.05]
    # (and, in fact, for most lambda).
    threshold = 1.61e-3

    # Play games with masks to implement the conditional evaluation of
    # the distribution.
    # lambda < -0.5:  var = nan
    low_mask = lam < -0.5
    # lambda == -0.5: var = inf
    neghalf_mask = lam == -0.5
    # abs(lambda) < threshold:  var = _tlvar0(lambda) (Taylor polynomial)
    small_mask = np.abs(lam) < threshold
    # else the "regular" case:  use the explicit formula.
    reg_mask = ~(low_mask | neghalf_mask | small_mask)

    # Get the 'lam' values for the cases where they are needed.
    small = lam[small_mask]
    reg = lam[reg_mask]

    # Compute the function for each case.
    v = np.empty_like(lam)
    v[low_mask] = np.nan
    v[neghalf_mask] = np.inf
    if small.size > 0:
        v[small_mask] = _tukeylambda_var0(small, order=3)
    if reg.size > 0:
        v[reg_mask] = ( (2.0 / (reg**2))
                        * (1.0 / (1 + 2*reg)
                           - gamma(reg + 1)**2 / gamma(2*reg + 2)) )
    v.shape = shp
    return v


_interp_lambda_values = \
    array([-0.065000000000000002, -0.06459987106934198 , -0.063409336779592498,
           -0.061457712036121957, -0.058793052317185797, -0.055480970388562798,
           -0.051603020699505384, -0.047254691241535272, -0.042543052317185789,
           -0.037584120113807504, -0.032500000000000008, -0.027415879886192512,
           -0.022456947682814206, -0.01774530875846473 , -0.013396979300494632,
           -0.009519029611437205, -0.006206947682814205, -0.003542287963878045,
           -0.001590663220407504, -0.000400128930658022,  0.                  ,
            0.000710202976151315,  0.002809772626615472,  0.006206947682814208,
            0.010753255293337104,  0.016249999999999997,  0.022456947682814209,
            0.029102824943801266,  0.035897175056198737,  0.042543052317185796,
            0.048749999999999995,  0.054246744706662881,  0.058793052317185797,
            0.062190227373384532,  0.064289797023848683,  0.065000000000000002])

# These values were computed with mpmath.  These are the values of the
# kurtosis of the Tukey Lambda distribution at the values of lambda given
# in _interp_lambda_values.
_interp_kurtosis_values = \
    array([ 2.523045133845003285,  2.511668653927227268,  2.478142819737123581,
            2.424210877953232757,  2.352566761905168669,  2.266581708755009306,
            2.169995568242775175,  2.066620913287312344,  1.960096533678926667,
            1.853710289571115499,  1.750295149431033481,  1.652190218317667325,
            1.561251954554219834,  1.478898978728214075,  1.406175372864659989,
            1.343820517594302766,  1.292337050298190348,  1.252051673452250435,
            1.223165956380509289,  1.205795902428776323,  1.199999999999999956,
            1.189764208874568352,  1.159884644242310747,  1.112711345736483004,
            1.051755332560134981,  0.981198954722772632,  0.905396678790319709,
            0.828462715224070445,  0.753996213621588418,  0.684948822710563787,
            0.623608257360730001,  0.571658672469203544,  0.530279768797318307,
            0.500254854328393672,  0.482067978778517481,  0.475978600169728649])

# Create the interpolator for the Tukey Lambda kurtosis for values of
# lambda near 0.
_tukeylambda_kurtosis_interp = interp1d(_interp_lambda_values,
                                        _interp_kurtosis_values,
                                        kind=5)


def tukeylambda_kurtosis(lam):
    """Kurtosis of the Tukey Lambda distribution.

    Parameters
    ----------
    lam : array_like
        The lambda values at which to compute the variance.

    Returns
    -------
    v : ndarray
        The variance.  For lam < -0.25, the variance is not defined, so
        np.nan is returned.  For lam = 0.25, np.inf is returned.

    Raises
    ------
    NotImplementedError if `lam` > 26.

    Notes
    -----
    In the interval -0.065 < lambda < 0.065, this function uses interpolation
    (based on accurate values that were previously computed) to compute the
    kurtosis.  Otherwise it uses the standard formula, given at
        http://en.wikipedia.org/wiki/Tukey_lambda_distribution).
    The interpolation is implemented to give an absolute error of less than
    1e-10 in this interval.
    """
    # Note: interpolation was used rather than the Taylor polynomial because
    # the Taylor polynomial diverged too from the exact formula too quickly.
    lam = np.asarray(lam)
    shp = lam.shape
    lam = np.atleast_1d(lam).astype(np.float64)

    if np.any(lam) > 26:
        raise NotImplementedError("this function cannot compute the "
                                  "kurtosis for lam > 26")


    # Use masks to implement the conditional evaluation of the kurtosis.
    # lambda < -0.25:  kurtosis = nan
    low_mask = lam < -0.25
    # lambda == -0.25: kurtosis = inf
    negqrtr_mask = lam == -0.25
    # lambda near 0:  var = _tukeylambda_kurtosis_interp(lambda)
    small_mask = ((lam > _interp_lambda_values[0]) &
                  (lam < _interp_lambda_values[-1]))
    # else the "regular" case:  use the explicit formula.
    reg_mask = ~(low_mask | negqrtr_mask | small_mask)

    # Get the 'lam' values for the cases where they are needed.
    small = lam[small_mask]
    reg = lam[reg_mask]

    # Compute the function for each case.
    k = np.empty_like(lam)
    k[low_mask] = np.nan
    k[negqrtr_mask] = np.inf
    if small.size > 0:
        k[small_mask] = _tukeylambda_kurtosis_interp(small)
    if reg.size > 0:
        g1 = gamma(reg + 1)
        g2 = gamma(2*reg + 1)
        g3 = gamma(3*reg + 1)
        g4 = gamma(4*reg + 1)
        k[reg_mask] = ((2*reg + 1)**2 * g2**2 * (3*g2**2 - 4*g1*g3 + g4)
                        / 2 / (4*reg + 1) / g4 / (g1**2 - g2)**2 - 3)

    # The return value will be a numpy array; resetting the shape ensures that
    # if `lam` was a scalar, the return value is a 0-d array.
    k.shape = shp
    return k
