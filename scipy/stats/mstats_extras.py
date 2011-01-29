"""
Additional statistics functions, with support to MA.

:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_edu
:date: $Date: 2007-10-29 17:18:13 +0200 (Mon, 29 Oct 2007) $
:version: $Id: morestats.py 3473 2007-10-29 15:18:13Z jarrod.millman $
"""
__author__ = "Pierre GF Gerard-Marchant"
__docformat__ = "restructuredtext en"


__all__ = ['compare_medians_ms',
           'hdquantiles', 'hdmedian', 'hdquantiles_sd',
           'idealfourths',
           'median_cihs','mjci','mquantiles_cimj',
           'rsh',
           'trimmed_mean_ci',]

import numpy as np
from numpy import float_, int_, ndarray

import numpy.ma as ma
from numpy.ma import MaskedArray

import mstats_basic as mstats

from scipy.stats.distributions import norm, beta, t, binom


#####--------------------------------------------------------------------------
#---- --- Quantiles ---
#####--------------------------------------------------------------------------
def hdquantiles(data, prob=list([.25,.5,.75]), axis=None, var=False,):
    """Computes quantile estimates with the Harrell-Davis method, where the estimates
are calculated as a weighted linear combination of order statistics.

Parameters
----------
    data: ndarray
        Data array.
    prob: sequence
        Sequence of quantiles to compute.
    axis : int
        Axis along which to compute the quantiles. If None, use a flattened array.
    var : boolean
        Whether to return the variance of the estimate.

Returns
-------
    A (p,) array of quantiles (if ``var`` is False), or a (2,p) array of quantiles
    and variances (if ``var`` is True), where ``p`` is the number of quantiles.

Notes
-----
    The function is restricted to 2D arrays.

    """
    def _hd_1D(data,prob,var):
        "Computes the HD quantiles for a 1D array. Returns nan for invalid data."
        xsorted = np.squeeze(np.sort(data.compressed().view(ndarray)))
        # Don't use length here, in case we have a numpy scalar
        n = xsorted.size
        #.........
        hd = np.empty((2,len(prob)), float_)
        if n < 2:
            hd.flat = np.nan
            if var:
                return hd
            return hd[0]
        #.........
        v = np.arange(n+1) / float(n)
        betacdf = beta.cdf
        for (i,p) in enumerate(prob):
            _w = betacdf(v, (n+1)*p, (n+1)*(1-p))
            w = _w[1:] - _w[:-1]
            hd_mean = np.dot(w, xsorted)
            hd[0,i] = hd_mean
            #
            hd[1,i] = np.dot(w, (xsorted-hd_mean)**2)
            #
        hd[0, prob == 0] = xsorted[0]
        hd[0, prob == 1] = xsorted[-1]
        if var:
            hd[1, prob == 0] = hd[1, prob == 1] = np.nan
            return hd
        return hd[0]
    # Initialization & checks ---------
    data = ma.array(data, copy=False, dtype=float_)
    p = np.array(prob, copy=False, ndmin=1)
    # Computes quantiles along axis (or globally)
    if (axis is None) or (data.ndim == 1):
        result = _hd_1D(data, p, var)
    else:
        if data.ndim > 2:
            raise ValueError("Array 'data' must be at most two dimensional, but got data.ndim = %d" % data.ndim)
        result = ma.apply_along_axis(_hd_1D, axis, data, p, var)
    #
    return ma.fix_invalid(result, copy=False)

#..............................................................................
def hdmedian(data, axis=-1, var=False):
    """Returns the Harrell-Davis estimate of the median along the given axis.

Parameters
----------
    data: ndarray
        Data array.
    axis : int
        Axis along which to compute the quantiles. If None, use a flattened array.
    var : boolean
        Whether to return the variance of the estimate.

    """
    result = hdquantiles(data,[0.5], axis=axis, var=var)
    return result.squeeze()


#..............................................................................
def hdquantiles_sd(data, prob=list([.25,.5,.75]), axis=None):
    """Computes the standard error of the Harrell-Davis quantile estimates by jackknife.


Parameters
----------
    data: ndarray
        Data array.
    prob: sequence
        Sequence of quantiles to compute.
    axis : int
        Axis along which to compute the quantiles. If None, use a flattened array.

Notes
-----
    The function is restricted to 2D arrays.

    """
    def _hdsd_1D(data,prob):
        "Computes the std error for 1D arrays."
        xsorted = np.sort(data.compressed())
        n = len(xsorted)
        #.........
        hdsd = np.empty(len(prob), float_)
        if n < 2:
            hdsd.flat = np.nan
        #.........
        vv = np.arange(n) / float(n-1)
        betacdf = beta.cdf
        #
        for (i,p) in enumerate(prob):
            _w = betacdf(vv, (n+1)*p, (n+1)*(1-p))
            w = _w[1:] - _w[:-1]
            mx_ = np.fromiter([np.dot(w,xsorted[np.r_[range(0,k),
                                                      range(k+1,n)].astype(int_)])
                                  for k in range(n)], dtype=float_)
            mx_var = np.array(mx_.var(), copy=False, ndmin=1) * n / float(n-1)
            hdsd[i] = float(n-1) * np.sqrt(np.diag(mx_var).diagonal() / float(n))
        return hdsd
    # Initialization & checks ---------
    data = ma.array(data, copy=False, dtype=float_)
    p = np.array(prob, copy=False, ndmin=1)
    # Computes quantiles along axis (or globally)
    if (axis is None):
        result = _hdsd_1D(data, p)
    else:
        if data.ndim > 2:
            raise ValueError("Array 'data' must be at most two dimensional, but got data.ndim = %d" % data.ndim)
        result = ma.apply_along_axis(_hdsd_1D, axis, data, p)
    #
    return ma.fix_invalid(result, copy=False).ravel()


#####--------------------------------------------------------------------------
#---- --- Confidence intervals ---
#####--------------------------------------------------------------------------

def trimmed_mean_ci(data, limits=(0.2,0.2), inclusive=(True,True),
                    alpha=0.05, axis=None):
    """Returns the selected confidence interval of the trimmed mean along the
given axis.

Parameters
----------
    data : sequence
        Input data. The data is transformed to a masked array
    proportiontocut : float
        Proportion of the data to cut from each side of the data .
        As a result, (2*proportiontocut*n) values are actually trimmed.
    alpha : float
        Confidence level of the intervals.
    inclusive : tuple of boolean
        If relative==False, tuple indicating whether values exactly equal to the
        absolute limits are allowed.
        If relative==True, tuple indicating whether the number of data being masked
        on each side should be rounded (True) or truncated (False).
    axis : int
        Axis along which to cut. If None, uses a flattened version of the input.

    """
    data = ma.array(data, copy=False)
    trimmed = mstats.trimr(data, limits=limits, inclusive=inclusive, axis=axis)
    tmean = trimmed.mean(axis)
    tstde = mstats.trimmed_stde(data,limits=limits,inclusive=inclusive,axis=axis)
    df = trimmed.count(axis) - 1
    tppf = t.ppf(1-alpha/2.,df)
    return np.array((tmean - tppf*tstde, tmean+tppf*tstde))

#..............................................................................
def mjci(data, prob=[0.25,0.5,0.75], axis=None):
    """Returns the Maritz-Jarrett estimators of the standard error of selected
experimental quantiles of the data.

Parameters
-----------
    data: ndarray
        Data array.
    prob: sequence
        Sequence of quantiles to compute.
    axis : int
        Axis along which to compute the quantiles. If None, use a flattened array.

    """
    def _mjci_1D(data, p):
        data = np.sort(data.compressed())
        n = data.size
        prob = (np.array(p) * n + 0.5).astype(int_)
        betacdf = beta.cdf
        #
        mj = np.empty(len(prob), float_)
        x = np.arange(1,n+1, dtype=float_) / n
        y = x - 1./n
        for (i,m) in enumerate(prob):
            (m1,m2) = (m-1, n-m)
            W = betacdf(x,m-1,n-m) - betacdf(y,m-1,n-m)
            C1 = np.dot(W,data)
            C2 = np.dot(W,data**2)
            mj[i] = np.sqrt(C2 - C1**2)
        return mj
    #
    data = ma.array(data, copy=False)
    if data.ndim > 2:
        raise ValueError("Array 'data' must be at most two dimensional, but got data.ndim = %d" % data.ndim)
    p = np.array(prob, copy=False, ndmin=1)
    # Computes quantiles along axis (or globally)
    if (axis is None):
        return _mjci_1D(data, p)
    else:
        return ma.apply_along_axis(_mjci_1D, axis, data, p)

#..............................................................................
def mquantiles_cimj(data, prob=[0.25,0.50,0.75], alpha=0.05, axis=None):
    """
    Computes the alpha confidence interval for the selected quantiles of the
    data, with Maritz-Jarrett estimators.

    Parameters
    ----------
    data: ndarray
        Data array.
    prob: sequence
        Sequence of quantiles to compute.
    alpha : float
        Confidence level of the intervals.
    axis : integer
        Axis along which to compute the quantiles.
        If None, use a flattened array.

    """
    alpha = min(alpha, 1-alpha)
    z = norm.ppf(1-alpha/2.)
    xq = mstats.mquantiles(data, prob, alphap=0, betap=0, axis=axis)
    smj = mjci(data, prob, axis=axis)
    return (xq - z * smj, xq + z * smj)


#.............................................................................
def median_cihs(data, alpha=0.05, axis=None):
    """Computes the alpha-level confidence interval for the median of the data,
following the Hettmasperger-Sheather method.

Parameters
----------
    data : sequence
        Input data. Masked values are discarded. The input should be 1D only, or
        axis should be set to None.
    alpha : float
        Confidence level of the intervals.
    axis : integer
        Axis along which to compute the quantiles. If None, use a flattened array.
    """
    def _cihs_1D(data, alpha):
        data = np.sort(data.compressed())
        n = len(data)
        alpha = min(alpha, 1-alpha)
        k = int(binom._ppf(alpha/2., n, 0.5))
        gk = binom.cdf(n-k,n,0.5) - binom.cdf(k-1,n,0.5)
        if gk < 1-alpha:
            k -= 1
            gk = binom.cdf(n-k,n,0.5) - binom.cdf(k-1,n,0.5)
        gkk = binom.cdf(n-k-1,n,0.5) - binom.cdf(k,n,0.5)
        I = (gk - 1 + alpha)/(gk - gkk)
        lambd = (n-k) * I / float(k + (n-2*k)*I)
        lims = (lambd*data[k] + (1-lambd)*data[k-1],
                lambd*data[n-k-1] + (1-lambd)*data[n-k])
        return lims
    data = ma.rray(data, copy=False)
    # Computes quantiles along axis (or globally)
    if (axis is None):
        result = _cihs_1D(data.compressed(), alpha)
    else:
        if data.ndim > 2:
            raise ValueError("Array 'data' must be at most two dimensional, but got data.ndim = %d" % data.ndim)
        result = ma.apply_along_axis(_cihs_1D, axis, data, alpha)
    #
    return result

#..............................................................................
def compare_medians_ms(group_1, group_2, axis=None):
    """Compares the medians from two independent groups along the given axis.

The comparison is performed using the McKean-Schrader estimate of the standard
error of the medians.

Parameters
----------
    group_1 : {sequence}
        First dataset.
    group_2 : {sequence}
        Second dataset.
    axis : {integer}
        Axis along which the medians are estimated. If None, the arrays are flattened.

Returns
-------
    A (p,) array of comparison values.

    """
    (med_1, med_2) = (ma.median(group_1,axis=axis), ma.median(group_2,axis=axis))
    (std_1, std_2) = (mstats.stde_median(group_1, axis=axis),
                      mstats.stde_median(group_2, axis=axis))
    W = np.abs(med_1 - med_2) / ma.sqrt(std_1**2 + std_2**2)
    return 1 - norm.cdf(W)


def idealfourths(data, axis=None):
    """Returns an estimate of the lower and upper quartiles of the data along
    the given axis, as computed with the ideal fourths.
    """
    def _idf(data):
        x = data.compressed()
        n = len(x)
        if n < 3:
            return [np.nan,np.nan]
        (j,h) = divmod(n/4. + 5/12.,1)
        qlo = (1-h)*x[j-1] + h*x[j]
        k = n - j
        qup = (1-h)*x[k] + h*x[k-1]
        return [qlo, qup]
    data = ma.sort(data, axis=axis).view(MaskedArray)
    if (axis is None):
        return _idf(data)
    else:
        return ma.apply_along_axis(_idf, axis, data)


def rsh(data, points=None):
    """Evaluates Rosenblatt's shifted histogram estimators for each point
on the dataset 'data'.

Parameters
    data : sequence
        Input data. Masked values are ignored.
    points : sequence
        Sequence of points where to evaluate Rosenblatt shifted histogram.
        If None, use the data.
    """
    data = ma.array(data, copy=False)
    if points is None:
        points = data
    else:
        points = np.array(points, copy=False, ndmin=1)
    if data.ndim != 1:
        raise AttributeError("The input array should be 1D only !")
    n = data.count()
    r = idealfourths(data, axis=None)
    h = 1.2 * (r[-1]-r[0]) / n**(1./5)
    nhi = (data[:,None] <= points[None,:] + h).sum(0)
    nlo = (data[:,None] < points[None,:] - h).sum(0)
    return (nhi-nlo) / (2.*n*h)


###############################################################################
