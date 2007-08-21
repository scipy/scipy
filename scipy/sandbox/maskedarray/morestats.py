"""
Generic statistics functions, with support to MA.

:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_edu
:date: $Date$
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'


import numpy
from numpy import bool_, float_, int_, ndarray, \
    sqrt,\
    arange, empty,\
    r_
from numpy import array as narray
import numpy.core.numeric as numeric
from numpy.core.numeric import concatenate

import maskedarray as MA
from maskedarray.core import masked, nomask, MaskedArray, masked_array
from maskedarray.extras import apply_along_axis, dot
from maskedarray.mstats import trim_both, trimmed_stde, mquantiles, mmedian, stde_median

from scipy.stats.distributions import norm, beta, t, binom
from scipy.stats.morestats import find_repeats

__all__ = ['hdquantiles', 'hdquantiles_sd',
           'trimmed_mean_ci', 'mjci', 'rank_data']


#####--------------------------------------------------------------------------
#---- --- Quantiles ---
#####--------------------------------------------------------------------------
def hdquantiles(data, prob=list([.25,.5,.75]), axis=None, var=False,):
    """Computes quantile estimates with the Harrell-Davis method, where the estimates
    are calculated as a weighted linear combination of order statistics.
    If var=True, the variance of the estimate is also returned. 
    Depending on var, returns a (p,) array of quantiles or a (2,p) array of quantiles
    and variances.
    
:Inputs:
    data: ndarray
        Data array.    
    prob: Sequence
        List of quantiles to compute.
    axis : integer *[None]*
        Axis along which to compute the quantiles. If None, use a flattened array.
    var : boolean *[False]*
        Whether to return the variance of the estimate.
        
:Note:
    The function is restricted to 2D arrays.
    """
    def _hd_1D(data,prob,var):
        "Computes the HD quantiles for a 1D array."
        xsorted = numpy.squeeze(numpy.sort(data.compressed().view(ndarray)))
        n = len(xsorted)
        #.........
        hd = empty((2,len(prob)), float_)
        if n < 2:
            hd.flat = numpy.nan
            if var:
                return hd
            return hd[0]
        #......... 
        v = arange(n+1) / float(n)
        betacdf = beta.cdf
        for (i,p) in enumerate(prob):    
            _w = betacdf(v, (n+1)*p, (n+1)*(1-p))
            w = _w[1:] - _w[:-1]
            hd_mean = dot(w, xsorted)
            hd[0,i] = hd_mean
            #
            hd[1,i] = dot(w, (xsorted-hd_mean)**2)
            #
        hd[0, prob == 0] = xsorted[0]
        hd[0, prob == 1] = xsorted[-1]  
        if var:  
            hd[1, prob == 0] = hd[1, prob == 1] = numpy.nan
            return hd
        return hd[0]
    # Initialization & checks ---------
    data = masked_array(data, copy=False, dtype=float_)
    p = numpy.array(prob, copy=False, ndmin=1)
    # Computes quantiles along axis (or globally)
    if (axis is None): 
        result = _hd_1D(data, p, var)
    else:
        assert data.ndim <= 2, "Array should be 2D at most !"
        result = apply_along_axis(_hd_1D, axis, data, p, var)
    #
    return masked_array(result, mask=numpy.isnan(result))
    
#..............................................................................
def hdquantiles_sd(data, prob=list([.25,.5,.75]), axis=None):
    """Computes the standard error of the Harrell-Davis quantile estimates by jackknife.
    
:Inputs:
    data: ndarray
        Data array.    
    prob: Sequence
        List of quantiles to compute.
    axis : integer *[None]*
        Axis along which to compute the quantiles. If None, use a flattened array.
    var : boolean *[False]*
        Whether to return the variance of the estimate.
    stderr : boolean *[False]*
        Whether to return the standard error of the estimate.
        
:Note:
    The function is restricted to 2D arrays.
    """  
    def _hdsd_1D(data,prob):
        "Computes the std error for 1D arrays."
        xsorted = numpy.sort(data.compressed())
        n = len(xsorted)
        #.........
        hdsd = empty(len(prob), float_)
        if n < 2:
            hdsd.flat = numpy.nan
        #......... 
        vv = arange(n) / float(n-1)
        betacdf = beta.cdf
        #
        for (i,p) in enumerate(prob):    
            _w = betacdf(vv, (n+1)*p, (n+1)*(1-p)) 
            w = _w[1:] - _w[:-1]
            mx_ = numpy.fromiter([dot(w,xsorted[r_[range(0,k),
                                                   range(k+1,n)].astype(int_)])
                                  for k in range(n)], dtype=float_)
            mx_var = numpy.array(mx_.var(), copy=False, ndmin=1) * n / float(n-1)
            hdsd[i] = float(n-1) * sqrt(numpy.diag(mx_var).diagonal() / float(n))
        return hdsd
    # Initialization & checks ---------
    data = masked_array(data, copy=False, dtype=float_)
    p = numpy.array(prob, copy=False, ndmin=1)
    # Computes quantiles along axis (or globally)
    if (axis is None): 
        result = _hdsd_1D(data.compressed(), p)
    else:
        assert data.ndim <= 2, "Array should be 2D at most !"
        result = apply_along_axis(_hdsd_1D, axis, data, p)
    #
    return masked_array(result, mask=numpy.isnan(result)).ravel()


#####--------------------------------------------------------------------------
#---- --- Confidence intervals ---
#####--------------------------------------------------------------------------

def trimmed_mean_ci(data, proportiontocut=0.2, alpha=0.05, axis=None):
    """Returns the selected confidence interval of the trimmed mean along the
    given axis.
    
:Inputs:
    data : sequence
        Input data. The data is transformed to a masked array
    proportiontocut : float *[0.2]*
        Proportion of the data to cut from each side of the data . 
        As a result, (2*proportiontocut*n) values are actually trimmed.
    alpha : float *[0.05]*
        Confidence level of the intervals
    axis : integer *[None]*
        Axis along which to cut.
    """
    data = masked_array(data, copy=False)
    trimmed = trim_both(data, proportiontocut=proportiontocut, axis=axis)
    tmean = trimmed.mean(axis)
    tstde = trimmed_stde(data, proportiontocut=proportiontocut, axis=axis)
    df = trimmed.count(axis) - 1
    tppf = t.ppf(1-alpha/2.,df)
    return numpy.array((tmean - tppf*tstde, tmean+tppf*tstde))

#..............................................................................
def mjci(data, prob=[0.25,0.5,0.75], axis=None):
    """Returns the Maritz-Jarrett estimators of the standard error of selected 
    experimental quantiles of the data.
    
:Input:
    data : sequence
        Input data.
    prob : sequence *[0.25,0.5,0.75]*
        Sequence of quantiles whose standard error must be estimated.
    axis : integer *[None]*
        Axis along which to compute the standard error.
    """
    def _mjci_1D(data, p):
        data = data.compressed()
        sorted = numpy.sort(data)
        n = data.size
        prob = (numpy.array(p) * n + 0.5).astype(int_)
        betacdf = beta.cdf
        #
        mj = empty(len(prob), float_)
        x = arange(1,n+1, dtype=float_) / n
        y = x - 1./n
        for (i,m) in enumerate(prob):
            (m1,m2) = (m-1, n-m)
            W = betacdf(x,m-1,n-m) - betacdf(y,m-1,n-m)
            C1 = numpy.dot(W,sorted)
            C2 = numpy.dot(W,sorted**2)
            mj[i] = sqrt(C2 - C1**2)
        return mj
    #
    data = masked_array(data, copy=False)
    assert data.ndim <= 2, "Array should be 2D at most !"
    p = numpy.array(prob, copy=False, ndmin=1)
    # Computes quantiles along axis (or globally)
    if (axis is None): 
        return _mjci_1D(data, p)
    else:
        return apply_along_axis(_mjci_1D, axis, data, p) 

#..............................................................................
def mquantiles_cimj(data, prob=[0.25,0.50,0.75], alpha=0.05, axis=None):
    """Computes the alpha confidence interval for the selected quantiles of the
    data, with Maritz-Jarrett estimators.
    
:Input:
    data : sequence
        Input data.
    prob : sequence *[0.25,0.5,0.75]*
        Sequence of quantiles whose standard error must be estimated.
    alpha : float *[0.05]*
        Confidence degree.
    axis : integer *[None]*
        Axis along which to compute the standard error.
    """
    alpha = min(alpha, 1-alpha)
    z = norm.ppf(1-alpha/2.)
    xq = mquantiles(data, prob, alphap=0, betap=0, axis=axis)
    smj = mjci(data, prob, axis=axis)
    return (xq - z * smj, xq + z * smj)

    
#.............................................................................
def median_cihs(data, alpha=0.05, axis=None):
    """Computes the alpha-level confidence interval for the median of the data,
    following the Hettmasperger-Sheather method.
    
:Inputs:
    data : sequence
        Input data. Masked values are discarded. The input should be 1D only
    alpha : float *[0.05]*
        Confidence degree.
    """
    def _cihs_1D(data, alpha):
        data = numpy.sort(data.compressed())
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
    data = masked_array(data, copy=False)
    # Computes quantiles along axis (or globally)
    if (axis is None): 
        result = _cihs_1D(data.compressed(), p, var)
    else:
        assert data.ndim <= 2, "Array should be 2D at most !"
        result = apply_along_axis(_cihs_1D, axis, data, alpha)
    #
    return result

#..............................................................................
def compare_medians_ms(group_1, group_2, axis=None):
    """Compares the medians from two independent groups along the given axis.
    Returns an array of p values.
    The comparison is performed using the McKean-Schrader estimate of the standard
    error of the medians.    
    
:Inputs:
    group_1 : sequence
        First dataset.
    group_2 : sequence
        Second dataset.
    axis : integer *[None]*
        Axis along which the medians are estimated. If None, the arrays are flattened.
    """
    (med_1, med_2) = (mmedian(group_1, axis=axis), mmedian(group_2, axis=axis))
    (std_1, std_2) = (stde_median(group_1, axis=axis), 
                      stde_median(group_2, axis=axis)) 
    W = abs(med_1 - med_2) / sqrt(std_1**2 + std_2**2)
    return 1 - norm.cdf(W)


#####--------------------------------------------------------------------------
#---- --- Ranking ---
#####--------------------------------------------------------------------------

#..............................................................................
def rank_data(data, axis=None, use_missing=False):
    """Returns the rank (also known as order statistics) of each data point 
    along the given axis.
    If some values are tied, their rank is averaged.
    If some values are masked, their rank is set to 0 if use_missing is False, or
    set to the average rank of the unmasked values if use_missing is True.
    
:Inputs:
    data : sequence
        Input data. The data is transformed to a masked array
    axis : integer *[None]*
        Axis along which to perform the ranking. If None, the array is first
        flattened. An exception is raised if the axis is specified for arrays
        with a dimension larger than 2
    use_missing : boolean *[False]*
        Flag indicating whether the masked values have a rank of 0 (False) or
        equal to the average rank of the unmasked values (True)    
    """
    #
    def _rank1d(data, use_missing=False):
        n = data.count()
        rk = numpy.empty(data.size, dtype=float_)
        idx = data.argsort()
        rk[idx[:n]] = numpy.arange(1,n+1)
        #
        if use_missing:
            rk[idx[n:]] = (n+1)/2.
        else:
            rk[idx[n:]] = 0
        #
        repeats = find_repeats(data)
        for r in repeats[0]:
            condition = (data==r).filled(False)
            rk[condition] = rk[condition].mean()
        return rk
    #
    data = masked_array(data, copy=False)
    if axis is None:
        if data.ndim > 1:
            return _rank1d(data.ravel(), use_missing).reshape(data.shape)
        else:
            return _rank1d(data, use_missing)
    else:
        return apply_along_axis(_rank1d, axis, data, use_missing) 

###############################################################################
if __name__ == '__main__':
    data = numpy.arange(100).reshape(4,25)
#    tmp = hdquantiles(data, prob=[0.25,0.75,0.5], axis=1, var=False)
    
