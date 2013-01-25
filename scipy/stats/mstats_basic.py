"""
An extension of scipy.stats.stats to support masked arrays

:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_edu
"""
#TODO : f_value_wilks_lambda looks botched... what are dfnum & dfden for ?
#TODO : ttest_reel looks botched:  what are x1,x2,v1,v2 for ?
#TODO : reimplement ksonesamp

from __future__ import division, print_function, absolute_import

__author__ = "Pierre GF Gerard-Marchant"
__docformat__ = "restructuredtext en"

__all__ = ['argstoarray',
           'betai',
           'chisquare','count_tied_groups',
           'describe',
           'f_oneway','f_value_wilks_lambda','find_repeats','friedmanchisquare',
           'gmean',
           'hmean',
           'kendalltau','kendalltau_seasonal','kruskal','kruskalwallis',
           'ks_twosamp','ks_2samp','kurtosis','kurtosistest',
           'linregress',
           'mannwhitneyu', 'meppf','mode','moment','mquantiles','msign',
           'normaltest',
           'obrientransform',
           'pearsonr','plotting_positions','pointbiserialr',
           'rankdata',
           'scoreatpercentile','sem',
           'sen_seasonal_slopes','signaltonoise','skew','skewtest','spearmanr',
           'theilslopes','threshold','tmax','tmean','tmin','trim','trimboth',
           'trimtail','trima','trimr','trimmed_mean','trimmed_std',
           'trimmed_stde','trimmed_var','tsem','ttest_1samp','ttest_onesamp',
           'ttest_ind','ttest_rel','tvar',
           'variation',
           'winsorize',
           'zmap', 'zscore'
           ]

import numpy as np
from numpy import ndarray
import numpy.ma as ma
from numpy.ma import MaskedArray, masked, nomask

from scipy.lib.six import iteritems

import itertools
import warnings


#import scipy.stats as stats
from . import stats
import scipy.special as special
import scipy.misc as misc
#import scipy.stats.futil as futil
from . import futil

genmissingvaldoc = """
Notes
-----
    Missing values are considered pair-wise: if a value is missing in x,
    the corresponding value in y is masked.
"""
#------------------------------------------------------------------------------
def _chk_asarray(a, axis):
    if axis is None:
        a = ma.ravel(a)
        outaxis = 0
    else:
        a = ma.asanyarray(a)
        outaxis = axis
    return a, outaxis

def _chk2_asarray(a, b, axis):
    if axis is None:
        a = ma.ravel(a)
        b = ma.ravel(b)
        outaxis = 0
    else:
        a = ma.asanyarray(a)
        b = ma.asanyarray(b)
        outaxis = axis
    return a, b, outaxis

def _chk_size(a,b):
    a = ma.asanyarray(a)
    b = ma.asanyarray(b)
    (na, nb) = (a.size, b.size)
    if na != nb:
        raise ValueError("The size of the input array should match!"\
                         " (%s <> %s)" % (na,nb))
    return (a,b,na)

def argstoarray(*args):
    """
    Constructs a 2D array from a group of sequences.

    Sequences are filled with missing values to match the length of the longest
    sequence.

    Parameters
    ----------
    args : sequences
        Group of sequences.

    Returns
    -------
    argstoarray : MaskedArray
        A ( `m` x `n` ) masked array, where `m` is the number of arguments and
        `n` the length of the longest argument.

    Notes
    -----
    numpy.ma.row_stack has identical behavior, but is called with a sequence of
    sequences.

    """
    if len(args) == 1 and not isinstance(args[0], ndarray):
        output = ma.asarray(args[0])
        if output.ndim != 2:
            raise ValueError("The input should be 2D")
    else:
        n = len(args)
        m = max([len(k) for k in args])
        output = ma.array(np.empty((n,m), dtype=float), mask=True)
        for (k,v) in enumerate(args):
            output[k,:len(v)] = v
    output[np.logical_not(np.isfinite(output._data))] = masked
    return output



#####--------------------------------------------------------------------------
#---- --- Ranking ---
#####--------------------------------------------------------------------------

def find_repeats(arr):
    """Find repeats in arr and return a tuple (repeats, repeat_count).
    Masked values are discarded.

Parameters
----------
    arr : sequence
        Input array. The array is flattened if it is not 1D.

Returns
-------
    repeats : ndarray
        Array of repeated values.
    counts : ndarray
        Array of counts.

    """
    marr = ma.compressed(arr)
    if not marr.size:
        return (np.array(0), np.array(0))
    (v1, v2, n) = futil.dfreps(ma.array(ma.compressed(arr), copy=True))
    return (v1[:n], v2[:n])


def count_tied_groups(x, use_missing=False):
    """
    Counts the number of tied values.

    Parameters
    ----------
    x : sequence
        Sequence of data on which to counts the ties
    use_missing : boolean
        Whether to consider missing values as tied.

    Returns
    -------
    count_tied_groups : dict
        Returns a dictionary (nb of ties: nb of groups).

    Examples
    --------
    >>> z = [0, 0, 0, 2, 2, 2, 3, 3, 4, 5, 6]
    >>> count_tied_groups(z)
    >>> {2:1, 3:2}
    >>> # The ties were 0 (3x), 2 (3x) and 3 (2x)
    >>> z = ma.array([0, 0, 1, 2, 2, 2, 3, 3, 4, 5, 6])
    >>> count_tied_groups(z)
    >>> {2:2, 3:1}
    >>> # The ties were 0 (2x), 2 (3x) and 3 (2x)
    >>> z[[1,-1]] = masked
    >>> count_tied_groups(z, use_missing=True)
    >>> {2:2, 3:1}
    >>> # The ties were 2 (3x), 3 (2x) and masked (2x)

    """
    nmasked = ma.getmask(x).sum()
    # We need the copy as find_repeats will overwrite the initial data
    data = ma.compressed(x).copy()
    (ties, counts) = find_repeats(data)
    nties = {}
    if len(ties):
        nties = dict(zip(np.unique(counts), itertools.repeat(1)))
        nties.update(dict(zip(*find_repeats(counts))))
    if nmasked and use_missing:
        try:
            nties[nmasked] += 1
        except KeyError:
            nties[nmasked] = 1
    return nties


def rankdata(data, axis=None, use_missing=False):
    """Returns the rank (also known as order statistics) of each data point
    along the given axis.

    If some values are tied, their rank is averaged.
    If some values are masked, their rank is set to 0 if use_missing is False,
    or set to the average rank of the unmasked values if use_missing is True.

    Parameters
    ----------
        data : sequence
            Input data. The data is transformed to a masked array
        axis : {None,int}, optional
            Axis along which to perform the ranking.
            If None, the array is first flattened. An exception is raised if
            the axis is specified for arrays with a dimension larger than 2
        use_missing : {boolean}, optional
            Whether the masked values have a rank of 0 (False) or equal to the
            average rank of the unmasked values (True).
    """
    #
    def _rank1d(data, use_missing=False):
        n = data.count()
        rk = np.empty(data.size, dtype=float)
        idx = data.argsort()
        rk[idx[:n]] = np.arange(1,n+1)
        #
        if use_missing:
            rk[idx[n:]] = (n+1)/2.
        else:
            rk[idx[n:]] = 0
        #
        repeats = find_repeats(data.copy())
        for r in repeats[0]:
            condition = (data==r).filled(False)
            rk[condition] = rk[condition].mean()
        return rk
    #
    data = ma.array(data, copy=False)
    if axis is None:
        if data.ndim > 1:
            return _rank1d(data.ravel(), use_missing).reshape(data.shape)
        else:
            return _rank1d(data, use_missing)
    else:
        return ma.apply_along_axis(_rank1d,axis,data,use_missing).view(ndarray)


#####--------------------------------------------------------------------------
#---- --- Central tendency ---
#####--------------------------------------------------------------------------

def gmean(a, axis=0):
    a, axis = _chk_asarray(a, axis)
    log_a = ma.log(a)
    return ma.exp(log_a.mean(axis=axis))
gmean.__doc__ = stats.gmean.__doc__


def hmean(a, axis=0):
    a, axis = _chk_asarray(a, axis)
    if isinstance(a, MaskedArray):
        size = a.count(axis)
    else:
        size = a.shape[axis]
    return size / (1.0/a).sum(axis)
hmean.__doc__ = stats.hmean.__doc__


def mode(a, axis=0):
    def _mode1D(a):
        (rep,cnt) = find_repeats(a)
        if not cnt.ndim:
            return (0, 0)
        elif cnt.size:
            return (rep[cnt.argmax()], cnt.max())
        return (a[0], 1)
    #
    if axis is None:
        output = _mode1D(ma.ravel(a))
        output = (ma.array(output[0]), ma.array(output[1]))
    else:
        output = ma.apply_along_axis(_mode1D, axis, a)
        newshape = list(a.shape)
        newshape[axis] = 1
        slices = [slice(None)] * output.ndim
        slices[axis] = 0
        modes = output[tuple(slices)].reshape(newshape)
        slices[axis] = 1
        counts = output[tuple(slices)].reshape(newshape)
        output = (modes, counts)
    return output
mode.__doc__ = stats.mode.__doc__


#####--------------------------------------------------------------------------
#---- --- Probabilities ---
#####--------------------------------------------------------------------------

def betai(a, b, x):
    x = np.asanyarray(x)
    x = ma.where(x < 1.0, x, 1.0)  # if x > 1 then return 1.0
    return special.betainc(a, b, x)
betai.__doc__ = stats.betai.__doc__


#####--------------------------------------------------------------------------
#---- --- Correlation ---
#####--------------------------------------------------------------------------

def msign(x):
    """Returns the sign of x, or 0 if x is masked."""
    return ma.filled(np.sign(x), 0)



def pearsonr(x,y):
    """Calculates a Pearson correlation coefficient and the p-value for testing
    non-correlation.

    The Pearson correlation coefficient measures the linear relationship
    between two datasets. Strictly speaking, Pearson's correlation requires
    that each dataset be normally distributed. Like other correlation
    coefficients, this one varies between -1 and +1 with 0 implying no
    correlation. Correlations of -1 or +1 imply an exact linear
    relationship. Positive correlations imply that as x increases, so does
    y. Negative correlations imply that as x increases, y decreases.

    The p-value roughly indicates the probability of an uncorrelated system
    producing datasets that have a Pearson correlation at least as extreme
    as the one computed from these datasets. The p-values are not entirely
    reliable but are probably reasonable for datasets larger than 500 or so.

    Parameters
    ----------
    x : 1D array
    y : 1D array the same length as x

    Returns
    -------
    (Pearson's correlation coefficient,
     2-tailed p-value)

    References
    ----------
    http://www.statsoft.com/textbook/glosp.html#Pearson%20Correlation
    """
    (x, y, n) = _chk_size(x, y)
    (x, y) = (x.ravel(), y.ravel())
    # Get the common mask and the total nb of unmasked elements
    m = ma.mask_or(ma.getmask(x), ma.getmask(y))
    n -= m.sum()
    df = n-2
    if df < 0:
        return (masked, masked)
    #
    (mx, my) = (x.mean(), y.mean())
    (xm, ym) = (x-mx, y-my)
    #
    r_num = n*(ma.add.reduce(xm*ym))
    r_den = n*ma.sqrt(ma.dot(xm,xm)*ma.dot(ym,ym))
    r = (r_num / r_den)
    # Presumably, if r > 1, then it is only some small artifact of floating
    # point arithmetic.
    r = min(r, 1.0)
    r = max(r, -1.0)
    df = n-2
    #
    t = ma.sqrt(df/((1.0-r)*(1.0+r))) * r
    if t is masked:
        prob = 0.
    else:
        prob = betai(0.5*df,0.5,df/(df+t*t))
    return (r,prob)


def spearmanr(x, y, use_ties=True):
    """Calculates a Spearman rank-order correlation coefficient and the p-value
    to test for non-correlation.

    The Spearman correlation is a nonparametric measure of the linear
    relationship between two datasets. Unlike the Pearson correlation, the
    Spearman correlation does not assume that both datasets are normally
    distributed. Like other correlation coefficients, this one varies
    between -1 and +1 with 0 implying no correlation. Correlations of -1 or
    +1 imply an exact linear relationship. Positive correlations imply that
    as x increases, so does y. Negative correlations imply that as x
    increases, y decreases.

    Missing values are discarded pair-wise: if a value is missing in x, the
    corresponding value in y is masked.

    The p-value roughly indicates the probability of an uncorrelated system
    producing datasets that have a Spearman correlation at least as extreme
    as the one computed from these datasets. The p-values are not entirely
    reliable but are probably reasonable for datasets larger than 500 or so.

Parameters
----------
    x : 1D array
    y : 1D array the same length as x
        The lengths of both arrays must be > 2.
    use_ties : {True, False}, optional
        Whether the correction for ties should be computed.

Returns
-------
    (Spearman correlation coefficient,
     2-tailed p-value)

    References
    ----------
    [CRCProbStat2000] section 14.7
    """
    (x, y, n) = _chk_size(x, y)
    (x, y) = (x.ravel(), y.ravel())
    #
    m = ma.mask_or(ma.getmask(x), ma.getmask(y))
    n -= m.sum()
    if m is not nomask:
        x = ma.array(x, mask=m, copy=True)
        y = ma.array(y, mask=m, copy=True)
    df = n-2
    if df < 0:
        raise ValueError("The input must have at least 3 entries!")
    # Gets the ranks and rank differences
    rankx = rankdata(x)
    ranky = rankdata(y)
    dsq = np.add.reduce((rankx-ranky)**2)
    # Tie correction
    if use_ties:
        xties = count_tied_groups(x)
        yties = count_tied_groups(y)
        corr_x = np.sum(v*k*(k**2-1) for (k,v) in iteritems(xties))/12.
        corr_y = np.sum(v*k*(k**2-1) for (k,v) in iteritems(yties))/12.
    else:
        corr_x = corr_y = 0
    denom = n*(n**2 - 1)/6.
    if corr_x != 0 or corr_y != 0:
        rho = denom - dsq - corr_x - corr_y
        rho /= ma.sqrt((denom-2*corr_x)*(denom-2*corr_y))
    else:
        rho = 1. - dsq/denom
    #
    t = ma.sqrt(ma.divide(df,(rho+1.0)*(1.0-rho))) * rho
    if t is masked:
        prob = 0.
    else:
        prob = betai(0.5*df,0.5,df/(df+t*t))
    return rho, prob


def kendalltau(x, y, use_ties=True, use_missing=False):
    """
    Computes Kendall's rank correlation tau on two variables *x* and *y*.

    Parameters
    ----------
    xdata : sequence
        First data list (for example, time).
    ydata : sequence
        Second data list.
    use_ties : {True, False}, optional
        Whether ties correction should be performed.
    use_missing : {False, True}, optional
        Whether missing data should be allocated a rank of 0 (False) or the
        average rank (True)

    Returns
    -------
    tau : float
        Kendall tau
    prob : float
        Approximate 2-side p-value.

    """
    (x, y, n) = _chk_size(x, y)
    (x, y) = (x.flatten(), y.flatten())
    m = ma.mask_or(ma.getmask(x), ma.getmask(y))
    if m is not nomask:
        x = ma.array(x, mask=m, copy=True)
        y = ma.array(y, mask=m, copy=True)
        n -= m.sum()
    #
    if n < 2:
        return (np.nan, np.nan)
    #
    rx = ma.masked_equal(rankdata(x, use_missing=use_missing), 0)
    ry = ma.masked_equal(rankdata(y, use_missing=use_missing), 0)
    idx = rx.argsort()
    (rx, ry) = (rx[idx], ry[idx])
    C = np.sum([((ry[i+1:]>ry[i]) * (rx[i+1:]>rx[i])).filled(0).sum()
                for i in range(len(ry)-1)], dtype=float)
    D = np.sum([((ry[i+1:]<ry[i])*(rx[i+1:]>rx[i])).filled(0).sum()
                for i in range(len(ry)-1)], dtype=float)
    if use_ties:
        xties = count_tied_groups(x)
        yties = count_tied_groups(y)
        corr_x = np.sum([v*k*(k-1) for (k,v) in iteritems(xties)], dtype=float)
        corr_y = np.sum([v*k*(k-1) for (k,v) in iteritems(yties)], dtype=float)
        denom = ma.sqrt((n*(n-1)-corr_x)/2. * (n*(n-1)-corr_y)/2.)
    else:
        denom = n*(n-1)/2.
    tau = (C-D) / denom
    #
    var_s = n*(n-1)*(2*n+5)
    if use_ties:
        var_s -= np.sum(v*k*(k-1)*(2*k+5)*1. for (k,v) in iteritems(xties))
        var_s -= np.sum(v*k*(k-1)*(2*k+5)*1. for (k,v) in iteritems(yties))
        v1 = np.sum([v*k*(k-1) for (k, v) in iteritems(xties)], dtype=float) *\
             np.sum([v*k*(k-1) for (k, v) in iteritems(yties)], dtype=float)
        v1 /= 2.*n*(n-1)
        if n > 2:
            v2 = np.sum([v*k*(k-1)*(k-2) for (k,v) in iteritems(xties)],
                        dtype=float) * \
                 np.sum([v*k*(k-1)*(k-2) for (k,v) in iteritems(yties)],
                        dtype=float)
            v2 /= 9.*n*(n-1)*(n-2)
        else:
            v2 = 0
    else:
        v1 = v2 = 0
    var_s /= 18.
    var_s += (v1 + v2)
    z = (C-D)/np.sqrt(var_s)
    prob = special.erfc(abs(z)/np.sqrt(2))
    return (tau, prob)


def kendalltau_seasonal(x):
    """Computes a multivariate extension Kendall's rank correlation tau, designed
    for seasonal data.

Parameters
----------
    x: 2D array
        Array of seasonal data, with seasons in columns.
    """
    x = ma.array(x, subok=True, copy=False, ndmin=2)
    (n,m) = x.shape
    n_p = x.count(0)
    #
    S_szn = np.sum(msign(x[i:]-x[i]).sum(0) for i in range(n))
    S_tot = S_szn.sum()
    #
    n_tot = x.count()
    ties = count_tied_groups(x.compressed())
    corr_ties =  np.sum(v*k*(k-1) for (k,v) in iteritems(ties))
    denom_tot = ma.sqrt(1.*n_tot*(n_tot-1)*(n_tot*(n_tot-1)-corr_ties))/2.
    #
    R = rankdata(x, axis=0, use_missing=True)
    K = ma.empty((m,m), dtype=int)
    covmat = ma.empty((m,m), dtype=float)
#    cov_jj = ma.empty(m, dtype=float)
    denom_szn = ma.empty(m, dtype=float)
    for j in range(m):
        ties_j = count_tied_groups(x[:,j].compressed())
        corr_j = np.sum(v*k*(k-1) for (k,v) in iteritems(ties_j))
        cmb = n_p[j]*(n_p[j]-1)
        for k in range(j,m,1):
            K[j,k] = np.sum(msign((x[i:,j]-x[i,j])*(x[i:,k]-x[i,k])).sum()
                               for i in range(n))
            covmat[j,k] = (K[j,k] +4*(R[:,j]*R[:,k]).sum() - \
                           n*(n_p[j]+1)*(n_p[k]+1))/3.
            K[k,j] = K[j,k]
            covmat[k,j] = covmat[j,k]
#        cov_jj[j] = (nn_p*(2*n_p[j]+5))
#        cov_jj[j] -= np.sum(v*k*(k-1)*(2*k+5) for (k,v) in ties_j.iteritems())
#        cov_jj[j] /= 18.
        denom_szn[j] = ma.sqrt(cmb*(cmb-corr_j)) / 2.
    var_szn = covmat.diagonal()
    #
    z_szn = msign(S_szn) * (abs(S_szn)-1) / ma.sqrt(var_szn)
    z_tot_ind = msign(S_tot) * (abs(S_tot)-1) / ma.sqrt(var_szn.sum())
    z_tot_dep = msign(S_tot) * (abs(S_tot)-1) / ma.sqrt(covmat.sum())
    #
    prob_szn = special.erfc(abs(z_szn)/np.sqrt(2))
    prob_tot_ind = special.erfc(abs(z_tot_ind)/np.sqrt(2))
    prob_tot_dep = special.erfc(abs(z_tot_dep)/np.sqrt(2))
    #
    chi2_tot = (z_szn*z_szn).sum()
    chi2_trd = m * z_szn.mean()**2
    output = {'seasonal tau': S_szn/denom_szn,
              'global tau': S_tot/denom_tot,
              'global tau (alt)': S_tot/denom_szn.sum(),
              'seasonal p-value': prob_szn,
              'global p-value (indep)': prob_tot_ind,
              'global p-value (dep)': prob_tot_dep,
              'chi2 total': chi2_tot,
              'chi2 trend': chi2_trd,
              }
    return output


def pointbiserialr(x, y):
    x = ma.fix_invalid(x, copy=True).astype(bool)
    y = ma.fix_invalid(y, copy=True).astype(float)
    # Get rid of the missing data ..........
    m = ma.mask_or(ma.getmask(x), ma.getmask(y))
    if m is not nomask:
        unmask = np.logical_not(m)
        x = x[unmask]
        y = y[unmask]
    #
    n = len(x)
    # phat is the fraction of x values that are True
    phat = x.sum() / float(n)
    y0 = y[~x]  # y-values where x is False
    y1 = y[x]  # y-values where x is True
    y0m = y0.mean()
    y1m = y1.mean()
    #
    rpb = (y1m - y0m)*np.sqrt(phat * (1-phat)) / y.std()
    #
    df = n-2
    t = rpb*ma.sqrt(df/(1.0-rpb**2))
    prob = betai(0.5*df, 0.5, df/(df+t*t))
    return rpb, prob

if stats.pointbiserialr.__doc__:
    pointbiserialr.__doc__ = stats.pointbiserialr.__doc__ + genmissingvaldoc


def linregress(*args):
    if len(args) == 1:  # more than 1D array?
        args = ma.array(args[0], copy=True)
        if len(args) == 2:
            x = args[0]
            y = args[1]
        else:
            x = args[:,0]
            y = args[:,1]
    else:
        x = ma.array(args[0]).flatten()
        y = ma.array(args[1]).flatten()
    m = ma.mask_or(ma.getmask(x), ma.getmask(y))
    if m is not nomask:
        x = ma.array(x,mask=m)
        y = ma.array(y,mask=m)
    n = len(x)
    (xmean, ymean) = (x.mean(), y.mean())
    (xm, ym) = (x-xmean, y-ymean)
    (Sxx, Syy) = (ma.add.reduce(xm*xm), ma.add.reduce(ym*ym))
    Sxy = ma.add.reduce(xm*ym)
    r_den = ma.sqrt(Sxx*Syy)
    if r_den == 0.0:
        r = 0.0
    else:
        r = Sxy / r_den
        if (r > 1.0):
            r = 1.0 # from numerical error
    #z = 0.5*log((1.0+r+TINY)/(1.0-r+TINY))
    df = n-2
    t = r * ma.sqrt(df/(1.0-r*r))
    prob = betai(0.5*df,0.5,df/(df+t*t))
    slope = Sxy / Sxx
    intercept = ymean - slope*xmean
    sterrest = ma.sqrt(1.-r*r) * y.std()
    return slope, intercept, r, prob, sterrest

if stats.linregress.__doc__:
    linregress.__doc__ = stats.linregress.__doc__ + genmissingvaldoc


def theilslopes(y, x=None, alpha=0.05):
    """Computes the Theil slope over the dataset (x,y), as the median of all slopes
    between paired values.

    Parameters
    ----------
        y : sequence
            Dependent variable.
        x : {None, sequence}, optional
            Independent variable. If None, use arange(len(y)) instead.
        alpha : float
            Confidence degree.

    Returns
    -------
        medslope : float
            Theil slope
        medintercept : float
            Intercept of the Theil line, as median(y)-medslope*median(x)
        lo_slope : float
            Lower bound of the confidence interval on medslope
        up_slope : float
            Upper bound of the confidence interval on medslope

    """
    y = ma.asarray(y).flatten()
    y[-1] = masked
    n = len(y)
    if x is None:
        x = ma.arange(len(y), dtype=float)
    else:
        x = ma.asarray(x).flatten()
        if len(x) != n:
            raise ValueError("Incompatible lengths ! (%s<>%s)" % (n,len(x)))
    m = ma.mask_or(ma.getmask(x), ma.getmask(y))
    y._mask = x._mask = m
    ny = y.count()
    #
    slopes = ma.hstack([(y[i+1:]-y[i])/(x[i+1:]-x[i]) for i in range(n-1)])
    slopes.sort()
    medslope = ma.median(slopes)
    medinter = ma.median(y) - medslope*ma.median(x)
    #
    if alpha > 0.5:
        alpha = 1.-alpha
    z = stats.distributions.norm.ppf(alpha/2.)
    #
    (xties, yties) = (count_tied_groups(x), count_tied_groups(y))
    nt = ny*(ny-1)/2.
    sigsq = (ny*(ny-1)*(2*ny+5)/18.)
    sigsq -= np.sum(v*k*(k-1)*(2*k+5) for (k,v) in iteritems(xties))
    sigsq -= np.sum(v*k*(k-1)*(2*k+5) for (k,v) in iteritems(yties))
    sigma = np.sqrt(sigsq)

    Ru = min(np.round((nt - z*sigma)/2. + 1), len(slopes)-1)
    Rl = max(np.round((nt + z*sigma)/2.), 0)
    delta = slopes[[Rl,Ru]]
    return medslope, medinter, delta[0], delta[1]


def sen_seasonal_slopes(x):
    x = ma.array(x, subok=True, copy=False, ndmin=2)
    (n,_) = x.shape
    # Get list of slopes per season
    szn_slopes = ma.vstack([(x[i+1:]-x[i])/np.arange(1,n-i)[:,None]
                            for i in range(n)])
    szn_medslopes = ma.median(szn_slopes, axis=0)
    medslope = ma.median(szn_slopes, axis=None)
    return szn_medslopes, medslope


#####--------------------------------------------------------------------------
#---- --- Inferential statistics ---
#####--------------------------------------------------------------------------

def ttest_onesamp(a, popmean):
    a = ma.asarray(a)
    x = a.mean(axis=None)
    v = a.var(axis=None,ddof=1)
    n = a.count(axis=None)
    df = n-1
    svar = ((n-1)*v) / float(df)
    t = (x-popmean)/ma.sqrt(svar*(1.0/n))
    prob = betai(0.5*df,0.5,df/(df+t*t))
    return t,prob
ttest_onesamp.__doc__ = stats.ttest_1samp.__doc__
ttest_1samp = ttest_onesamp


def ttest_ind(a, b, axis=0):
    a, b, axis = _chk2_asarray(a, b, axis)
    (x1, x2) = (a.mean(axis), b.mean(axis))
    (v1, v2) = (a.var(axis=axis, ddof=1), b.var(axis=axis, ddof=1))
    (n1, n2) = (a.count(axis), b.count(axis))
    df = n1+n2-2
    svar = ((n1-1)*v1+(n2-1)*v2) / float(df)
    svar == 0
    t = (x1-x2)/ma.sqrt(svar*(1.0/n1 + 1.0/n2))  # N-D COMPUTATION HERE!!!!!!
    t = ma.filled(t, 1)           # replace NaN t-values with 1.0
    probs = betai(0.5*df,0.5,float(df)/(df+t*t)).reshape(t.shape)
    return t, probs.squeeze()
ttest_ind.__doc__ = stats.ttest_ind.__doc__


def ttest_rel(a,b,axis=None):
    a, b, axis = _chk2_asarray(a, b, axis)
    if len(a)!=len(b):
        raise ValueError('unequal length arrays')
    (x1, x2) = (a.mean(axis), b.mean(axis))
    (v1, v2) = (a.var(axis=axis, ddof=1), b.var(axis=axis, ddof=1))
    n = a.count(axis)
    df = (n-1.0)
    d = (a-b).astype('d')
    denom = ma.sqrt((n*ma.add.reduce(d*d,axis) - ma.add.reduce(d,axis)**2) /df)
    #zerodivproblem = denom == 0
    t = ma.add.reduce(d, axis) / denom
    t = ma.filled(t, 1)
    probs = betai(0.5*df,0.5,df/(df+t*t)).reshape(t.shape).squeeze()
    return t, probs
ttest_rel.__doc__ = stats.ttest_rel.__doc__


def chisquare(f_obs, f_exp=None):
    f_obs = ma.asarray(f_obs)
    if f_exp is None:
        f_exp = ma.array([f_obs.mean(axis=0)] * len(f_obs))
    f_exp = f_exp.astype(float)
    chisq = ma.add.reduce((f_obs-f_exp)**2 / f_exp)
    return chisq, stats.chisqprob(chisq, f_obs.count(0)-1)
chisquare.__doc__ = stats.chisquare.__doc__


def mannwhitneyu(x,y, use_continuity=True):
    """
    Computes the Mann-Whitney statistic

    Missing values in `x` and/or `y` are discarded.

    Parameters
    ----------
    x : sequence
        Input
    y : sequence
        Input
    use_continuity : {True, False}, optional
        Whether a continuity correction (1/2.) should be taken into account.

    Returns
    -------
    u : float
        The Mann-Whitney statistics
    prob : float
        Approximate p-value assuming a normal distribution.

    """
    x = ma.asarray(x).compressed().view(ndarray)
    y = ma.asarray(y).compressed().view(ndarray)
    ranks = rankdata(np.concatenate([x,y]))
    (nx, ny) = (len(x), len(y))
    nt = nx + ny
    U = ranks[:nx].sum() - nx*(nx+1)/2.
    U = max(U, nx*ny - U)
    u = nx*ny - U
    #
    mu = (nx*ny)/2.
    sigsq = (nt**3 - nt)/12.
    ties = count_tied_groups(ranks)
    sigsq -= np.sum(v*(k**3-k) for (k,v) in iteritems(ties))/12.
    sigsq *= nx*ny/float(nt*(nt-1))
    #
    if use_continuity:
        z = (U - 1/2. - mu) / ma.sqrt(sigsq)
    else:
        z = (U - mu) / ma.sqrt(sigsq)
    prob = special.erfc(abs(z)/np.sqrt(2))
    return (u, prob)


def kruskalwallis(*args):
    output = argstoarray(*args)
    ranks = ma.masked_equal(rankdata(output, use_missing=False), 0)
    sumrk = ranks.sum(-1)
    ngrp = ranks.count(-1)
    ntot = ranks.count()
#    ssbg = (sumrk**2/ranks.count(-1)).sum() - ranks.sum()**2/ntotal
#    H = ssbg / (ntotal*(ntotal+1)/12.)
    H = 12./(ntot*(ntot+1)) * (sumrk**2/ngrp).sum() - 3*(ntot+1)
    # Tie correction
    ties = count_tied_groups(ranks)
    T = 1. - np.sum(v*(k**3-k) for (k,v) in iteritems(ties))/float(ntot**3-ntot)
    if T == 0:
        raise ValueError('All numbers are identical in kruskal')
    H /= T
    #
    df = len(output) - 1
    prob = stats.chisqprob(H,df)
    return (H, prob)
kruskal = kruskalwallis
kruskalwallis.__doc__ = stats.kruskal.__doc__


_kolmog2 = special.kolmogorov
def _kolmog1(x,n):
    if x <= 0:
        return 0
    if x >= 1:
        return 1
    j = np.arange(np.floor(n*(1-x))+1)
    return 1 - x * np.sum(np.exp(np.log(misc.comb(n,j))
                                       + (n-j) * np.log(1-x-j/float(n))
                                       + (j-1) * np.log(x+j/float(n))))


def ks_twosamp(data1, data2, alternative="two-sided"):
    """
    Computes the Kolmogorov-Smirnov test on two samples.

    Missing values are discarded.

    Parameters
    ----------
    data1 : array_like
        First data set
    data2 : array_like
        Second data set
    alternative : {'two-sided', 'less', 'greater'}, optional
        Indicates the alternative hypothesis.  Default is 'two-sided'.

    Returns
    -------
    d : float
        Value of the Kolmogorov Smirnov test
    p : float
        Corresponding p-value.

    """
    (data1, data2) = (ma.asarray(data1), ma.asarray(data2))
    (n1, n2) = (data1.count(), data2.count())
    n = (n1*n2/float(n1+n2))
    mix = ma.concatenate((data1.compressed(), data2.compressed()))
    mixsort = mix.argsort(kind='mergesort')
    csum = np.where(mixsort<n1, 1./n1, -1./n2).cumsum()
    # Check for ties
    if len(np.unique(mix)) < (n1+n2):
        csum = csum[np.r_[np.diff(mix[mixsort]).nonzero()[0],-1]]
    #
    alternative = str(alternative).lower()[0]
    if alternative == 't':
        d = ma.abs(csum).max()
        prob = _kolmog2(np.sqrt(n)*d)
    elif alternative == 'l':
        d = -csum.min()
        prob = np.exp(-2*n*d**2)
    elif alternative == 'g':
        d = csum.max()
        prob = np.exp(-2*n*d**2)
    else:
        raise ValueError("Invalid value for the alternative hypothesis: "\
                         "should be in 'two-sided', 'less' or 'greater'")
    return (d, prob)
ks_2samp = ks_twosamp


def ks_twosamp_old(data1, data2):
    """ Computes the Kolmogorov-Smirnov statistic on 2 samples.

    Returns
    -------
    KS D-value, p-value

    """
    (data1, data2) = [ma.asarray(d).compressed() for d in (data1,data2)]
    return stats.ks_2samp(data1,data2)


#####--------------------------------------------------------------------------
#---- --- Trimming ---
#####--------------------------------------------------------------------------

def threshold(a, threshmin=None, threshmax=None, newval=0):
    """Clip array to a given value.

    Similar to numpy.clip(), except that values less than threshmin or
    greater than threshmax are replaced by newval, instead of by
    threshmin and threshmax respectively.

    Parameters
    ----------
    a : ndarray
        Input data
    threshmin : {None, float}, optional
        Lower threshold. If None, set to the minimum value.
    threshmax : {None, float}, optional
        Upper threshold. If None, set to the maximum value.
    newval : {0, float}, optional
        Value outside the thresholds.

    Returns
    -------
    a, with values less (greater) than threshmin (threshmax) replaced with newval.

    """
    a = ma.array(a, copy=True)
    mask = np.zeros(a.shape, dtype=bool)
    if threshmin is not None:
        mask |= (a < threshmin).filled(False)
    if threshmax is not None:
        mask |= (a > threshmax).filled(False)
    a[mask] = newval
    return a


def trima(a, limits=None, inclusive=(True,True)):
    """Trims an array by masking the data outside some given limits.
    Returns a masked version of the input array.

    Parameters
    ----------
    a : sequence
        Input array.
    limits : {None, tuple}, optional
        Tuple of (lower limit, upper limit) in absolute values.
        Values of the input array lower (greater) than the lower (upper) limit
        will be masked. A limit is None indicates an open interval.
    inclusive : {(True,True) tuple}, optional
        Tuple of (lower flag, upper flag), indicating whether values exactly
        equal to the lower (upper) limit are allowed.

    """
    a = ma.asarray(a)
    a.unshare_mask()
    if limits is None:
        return a
    (lower_lim, upper_lim) = limits
    (lower_in, upper_in) = inclusive
    condition = False
    if lower_lim is not None:
        if lower_in:
            condition |= (a < lower_lim)
        else:
            condition |= (a <= lower_lim)
    if upper_lim is not None:
        if upper_in:
            condition |= (a > upper_lim)
        else:
            condition |= (a >= upper_lim)
    a[condition.filled(True)] = masked
    return a


def trimr(a, limits=None, inclusive=(True, True), axis=None):
    """
    Trims an array by masking some proportion of the data on each end.
    Returns a masked version of the input array.

    Parameters
    ----------
    a : sequence
        Input array.
    limits : {None, tuple}, optional
        Tuple of the percentages to cut on each side of the array, with respect
        to the number of unmasked data, as floats between 0. and 1.
        Noting n the number of unmasked data before trimming, the
        (n*limits[0])th smallest data and the (n*limits[1])th largest data are
        masked, and the total number of unmasked data after trimming is
        n*(1.-sum(limits)).  The value of one limit can be set to None to
        indicate an open interval.
    inclusive : {(True,True) tuple}, optional
        Tuple of flags indicating whether the number of data being masked on
        the left (right) end should be truncated (True) or rounded (False) to
        integers.
    axis : {None,int}, optional
        Axis along which to trim. If None, the whole array is trimmed, but its
        shape is maintained.

    """
    def _trimr1D(a, low_limit, up_limit, low_inclusive, up_inclusive):
        n = a.count()
        idx = a.argsort()
        if low_limit:
            if low_inclusive:
                lowidx = int(low_limit*n)
            else:
                lowidx = np.round(low_limit*n)
            a[idx[:lowidx]] = masked
        if up_limit is not None:
            if up_inclusive:
                upidx = n - int(n*up_limit)
            else:
                upidx = n- np.round(n*up_limit)
            a[idx[upidx:]] = masked
        return a
    #
    a = ma.asarray(a)
    a.unshare_mask()
    if limits is None:
        return a
    # Check the limits
    (lolim, uplim) = limits
    errmsg = "The proportion to cut from the %s should be between 0. and 1."
    if lolim is not None:
        if lolim > 1. or lolim < 0:
            raise ValueError(errmsg % 'beginning' + "(got %s)" % lolim)
    if uplim is not None:
        if uplim > 1. or uplim < 0:
            raise ValueError(errmsg % 'end' + "(got %s)" % uplim)
    #
    (loinc, upinc) = inclusive
    #
    if axis is None:
        shp = a.shape
        return _trimr1D(a.ravel(),lolim,uplim,loinc,upinc).reshape(shp)
    else:
        return ma.apply_along_axis(_trimr1D, axis, a, lolim,uplim,loinc,upinc)

trimdoc = """
    Parameters
    ----------
    a : sequence
        Input array
    limits : {None, tuple}, optional
        If `relative` is False, tuple (lower limit, upper limit) in absolute values.
        Values of the input array lower (greater) than the lower (upper) limit are
        masked.

        If `relative` is True, tuple (lower percentage, upper percentage) to cut
        on each side of the  array, with respect to the number of unmasked data.

        Noting n the number of unmasked data before trimming, the (n*limits[0])th
        smallest data and the (n*limits[1])th largest data are masked, and the
        total number of unmasked data after trimming is n*(1.-sum(limits))
        In each case, the value of one limit can be set to None to indicate an
        open interval.

        If limits is None, no trimming is performed
    inclusive : {(bool, bool) tuple}, optional
        If `relative` is False, tuple indicating whether values exactly equal
        to the absolute limits are allowed.
        If `relative` is True, tuple indicating whether the number of data
        being masked on each side should be rounded (True) or truncated
        (False).
    relative : bool, optional
        Whether to consider the limits as absolute values (False) or proportions
        to cut (True).
    axis : int, optional
        Axis along which to trim.
"""


def trim(a, limits=None, inclusive=(True,True), relative=False, axis=None):
    """
    Trims an array by masking the data outside some given limits.

    Returns a masked version of the input array.

    %s

    Examples
    --------
    >>> z = [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10]
    >>> trim(z,(3,8))
    [--,--, 3, 4, 5, 6, 7, 8,--,--]
    >>> trim(z,(0.1,0.2),relative=True)
    [--, 2, 3, 4, 5, 6, 7, 8,--,--]

    """
    if relative:
        return trimr(a, limits=limits, inclusive=inclusive, axis=axis)
    else:
        return trima(a, limits=limits, inclusive=inclusive)
trim.__doc__ = trim.__doc__ % trimdoc


def trimboth(data, proportiontocut=0.2, inclusive=(True,True), axis=None):
    """
    Trims the smallest and largest data values.

    Trims the `data` by masking the ``int(proportiontocut * n)`` smallest and
    ``int(proportiontocut * n)`` largest values of data along the given axis,
    where n is the number of unmasked values before trimming.

    Parameters
    ----------
    data : ndarray
        Data to trim.
    proportiontocut : float, optional
        Percentage of trimming (as a float between 0 and 1).
        If n is the number of unmasked values before trimming, the number of
        values after trimming is ``(1 - 2*proportiontocut) * n``.
        Default is 0.2.
    inclusive : {(bool, bool) tuple}, optional
        Tuple indicating whether the number of data being masked on each side
        should be rounded (True) or truncated (False).
    axis : int, optional
        Axis along which to perform the trimming.
        If None, the input array is first flattened.

    """
    return trimr(data, limits=(proportiontocut,proportiontocut),
                 inclusive=inclusive, axis=axis)

#..............................................................................
def trimtail(data, proportiontocut=0.2, tail='left', inclusive=(True,True),
             axis=None):
    """
    Trims the data by masking values from one tail.

    Parameters
    ----------
    data : array_like
        Data to trim.
    proportiontocut : float, optional
        Percentage of trimming. If n is the number of unmasked values
        before trimming, the number of values after trimming is
        ``(1 - proportiontocut) * n``.  Default is 0.2.
    tail : {'left','right'}, optional
        If 'left' the `proportiontocut` lowest values will be masked.
        If 'right' the `proportiontocut` highest values will be masked.
        Default is 'left'.
    inclusive : {(bool, bool) tuple}, optional
        Tuple indicating whether the number of data being masked on each side
        should be rounded (True) or truncated (False).  Default is
        (True, True).
    axis : int, optional
        Axis along which to perform the trimming.
        If None, the input array is first flattened.  Default is None.

    Returns
    -------
    trimtail : ndarray
        Returned array of same shape as `data` with masked tail values.

    """
    tail = str(tail).lower()[0]
    if tail == 'l':
        limits = (proportiontocut,None)
    elif tail == 'r':
        limits = (None, proportiontocut)
    else:
        raise TypeError("The tail argument should be in ('left','right')")
    return trimr(data, limits=limits, axis=axis, inclusive=inclusive)

trim1 = trimtail

def trimmed_mean(a, limits=(0.1,0.1), inclusive=(1,1), relative=True,
                 axis=None):
    """Returns the trimmed mean of the data along the given axis.

    %s

    """ % trimdoc
    if (not isinstance(limits,tuple)) and isinstance(limits,float):
        limits = (limits, limits)
    if relative:
        return trimr(a,limits=limits,inclusive=inclusive,axis=axis).mean(axis=axis)
    else:
        return trima(a,limits=limits,inclusive=inclusive).mean(axis=axis)


def trimmed_var(a, limits=(0.1,0.1), inclusive=(1,1), relative=True,
                axis=None, ddof=0):
    """Returns the trimmed variance of the data along the given axis.

    %s
    ddof : {0,integer}, optional
        Means Delta Degrees of Freedom. The denominator used during computations
        is (n-ddof). DDOF=0 corresponds to a biased estimate, DDOF=1 to an un-
        biased estimate of the variance.

    """ % trimdoc
    if (not isinstance(limits,tuple)) and isinstance(limits,float):
        limits = (limits, limits)
    if relative:
        out = trimr(a,limits=limits, inclusive=inclusive,axis=axis)
    else:
        out = trima(a,limits=limits,inclusive=inclusive)
    return out.var(axis=axis, ddof=ddof)


def trimmed_std(a, limits=(0.1,0.1), inclusive=(1,1), relative=True,
                axis=None, ddof=0):
    """Returns the trimmed standard deviation of the data along the given axis.

    %s
    ddof : {0,integer}, optional
        Means Delta Degrees of Freedom. The denominator used during computations
        is (n-ddof). DDOF=0 corresponds to a biased estimate, DDOF=1 to an un-
        biased estimate of the variance.

    """ % trimdoc
    if (not isinstance(limits,tuple)) and isinstance(limits,float):
        limits = (limits, limits)
    if relative:
        out = trimr(a,limits=limits,inclusive=inclusive,axis=axis)
    else:
        out = trima(a,limits=limits,inclusive=inclusive)
    return out.std(axis=axis,ddof=ddof)


def trimmed_stde(a, limits=(0.1,0.1), inclusive=(1,1), axis=None):
    """
    Returns the standard error of the trimmed mean along the given axis.

    Parameters
    ----------
    a : sequence
        Input array
    limits : {(0.1,0.1), tuple of float}, optional
        tuple (lower percentage, upper percentage) to cut  on each side of the
        array, with respect to the number of unmasked data.

        If n is the number of unmasked data before trimming, the values
        smaller than ``n * limits[0]`` and the values larger than
        ``n * `limits[1]`` are masked, and the total number of unmasked
        data after trimming is ``n * (1.-sum(limits))``.  In each case,
        the value of one limit can be set to None to indicate an open interval.
        If `limits` is None, no trimming is performed.
    inclusive : {(bool, bool) tuple} optional
        Tuple indicating whether the number of data being masked on each side
        should be rounded (True) or truncated (False).
    axis : int, optional
        Axis along which to trim.

    Returns
    -------
    trimmed_stde : scalar or ndarray

    """
    #........................
    def _trimmed_stde_1D(a, low_limit, up_limit, low_inclusive, up_inclusive):
        "Returns the standard error of the trimmed mean for a 1D input data."
        n = a.count()
        idx = a.argsort()
        if low_limit:
            if low_inclusive:
                lowidx = int(low_limit*n)
            else:
                lowidx = np.round(low_limit*n)
            a[idx[:lowidx]] = masked
        if up_limit is not None:
            if up_inclusive:
                upidx = n - int(n*up_limit)
            else:
                upidx = n- np.round(n*up_limit)
            a[idx[upidx:]] = masked
        nsize = a.count()
        a[idx[:lowidx]] = a[idx[lowidx]]
        a[idx[upidx:]] = a[idx[upidx-1]]
        winstd = a.std(ddof=1)
        return winstd / ((1-low_limit-up_limit)*np.sqrt(len(a)))
    #........................
    a = ma.array(a, copy=True, subok=True)
    a.unshare_mask()
    if limits is None:
        return a.std(axis=axis,ddof=1)/ma.sqrt(a.count(axis))
    if (not isinstance(limits,tuple)) and isinstance(limits,float):
        limits = (limits, limits)
    # Check the limits
    (lolim, uplim) = limits
    errmsg = "The proportion to cut from the %s should be between 0. and 1."
    if lolim is not None:
        if lolim > 1. or lolim < 0:
            raise ValueError(errmsg % 'beginning' + "(got %s)" % lolim)
    if uplim is not None:
        if uplim > 1. or uplim < 0:
            raise ValueError(errmsg % 'end' + "(got %s)" % uplim)
    #
    (loinc, upinc) = inclusive
    if (axis is None):
        shp = a.shape
        return _trimmed_stde_1D(a.ravel(),lolim,uplim,loinc,upinc)
    else:
        if a.ndim > 2:
            raise ValueError("Array 'a' must be at most two dimensional, but got a.ndim = %d" % a.ndim)
        return ma.apply_along_axis(_trimmed_stde_1D, axis, a,
                                   lolim,uplim,loinc,upinc)


def tmean(a, limits=None, inclusive=(True,True)):
    return trima(a, limits=limits, inclusive=inclusive).mean()
tmean.__doc__ = stats.tmean.__doc__

def tvar(a, limits=None, inclusive=(True,True)):
    return trima(a, limits=limits, inclusive=inclusive).var()
tvar.__doc__ = stats.tvar.__doc__

def tmin(a, lowerlimit=None, axis=0, inclusive=True):
    a, axis = _chk_asarray(a, axis)
    am = trima(a, (lowerlimit, None), (inclusive, False))
    return ma.minimum.reduce(am, axis)
tmin.__doc__  = stats.tmin.__doc__

def tmax(a, upperlimit, axis=0, inclusive=True):
    a, axis = _chk_asarray(a, axis)
    am = trima(a, (None, upperlimit), (False, inclusive))
    return ma.maximum.reduce(am, axis)
tmax.__doc__  = stats.tmax.__doc__

def tsem(a, limits=None, inclusive=(True,True)):
    a = ma.asarray(a).ravel()
    if limits is None:
        n = float(a.count())
        return a.std()/ma.sqrt(n)
    am = trima(a.ravel(), limits, inclusive)
    sd = np.sqrt(am.var())
    return sd / am.count()
tsem.__doc__ = stats.tsem.__doc__


def winsorize(a, limits=None, inclusive=(True,True), inplace=False, axis=None):
    """
    Returns a Winsorized version of the input array.

    The (limits[0])th lowest values are set to the (limits[0])th percentile,
    and the (limits[1])th highest values are set to the (limits[1])th
    percentile.
    Masked values are skipped.


    Parameters
    ----------
    a : sequence
        Input array.
    limits : {None, tuple of float}, optional
        Tuple of the percentages to cut on each side of the array, with respect
        to the number of unmasked data, as floats between 0. and 1.
        Noting n the number of unmasked data before trimming, the
        (n*limits[0])th smallest data and the (n*limits[1])th largest data are
        masked, and the total number of unmasked data after trimming
        is n*(1.-sum(limits)) The value of one limit can be set to None to
        indicate an open interval.
    inclusive : {(True, True) tuple}, optional
        Tuple indicating whether the number of data being masked on each side
        should be rounded (True) or truncated (False).
    inplace : {False, True}, optional
        Whether to winsorize in place (True) or to use a copy (False)
    axis : {None, int}, optional
        Axis along which to trim. If None, the whole array is trimmed, but its
        shape is maintained.

    """
    def _winsorize1D(a, low_limit, up_limit, low_include, up_include):
        n = a.count()
        idx = a.argsort()
        if low_limit:
            if low_include:
                lowidx = int(low_limit*n)
            else:
                lowidx = np.round(low_limit*n)
            a[idx[:lowidx]] = a[idx[lowidx]]
        if up_limit is not None:
            if up_include:
                upidx = n - int(n*up_limit)
            else:
                upidx = n- np.round(n*up_limit)
            a[idx[upidx:]] = a[idx[upidx-1]]
        return a
    # We gonna modify a: better make a copy
    a = ma.array(a, copy=np.logical_not(inplace))
    #
    if limits is None:
        return a
    if (not isinstance(limits,tuple)) and isinstance(limits,float):
        limits = (limits, limits)
    # Check the limits
    (lolim, uplim) = limits
    errmsg = "The proportion to cut from the %s should be between 0. and 1."
    if lolim is not None:
        if lolim > 1. or lolim < 0:
            raise ValueError(errmsg % 'beginning' + "(got %s)" % lolim)
    if uplim is not None:
        if uplim > 1. or uplim < 0:
            raise ValueError(errmsg % 'end' + "(got %s)" % uplim)
    #
    (loinc, upinc) = inclusive
    #
    if axis is None:
        shp = a.shape
        return _winsorize1D(a.ravel(),lolim,uplim,loinc,upinc).reshape(shp)
    else:
        return ma.apply_along_axis(_winsorize1D, axis,a,lolim,uplim,loinc,upinc)


#####--------------------------------------------------------------------------
#---- --- Moments ---
#####--------------------------------------------------------------------------

def moment(a, moment=1, axis=0):
    a, axis = _chk_asarray(a, axis)
    if moment == 1:
        # By definition the first moment about the mean is 0.
        shape = list(a.shape)
        del shape[axis]
        if shape:
            # return an actual array of the appropriate shape
            return np.zeros(shape, dtype=float)
        else:
            # the input was 1D, so return a scalar instead of a rank-0 array
            return np.float64(0.0)
    else:
        mn = ma.expand_dims(a.mean(axis=axis), axis)
        s = ma.power((a-mn), moment)
        return s.mean(axis=axis)
moment.__doc__ = stats.moment.__doc__


def variation(a, axis=0):
    a, axis = _chk_asarray(a, axis)
    return a.std(axis)/a.mean(axis)
variation.__doc__ = stats.variation.__doc__


def skew(a, axis=0, bias=True):
    a, axis = _chk_asarray(a,axis)
    n = a.count(axis)
    m2 = moment(a, 2, axis)
    m3 = moment(a, 3, axis)
    olderr = np.seterr(all='ignore')
    try:
        vals = ma.where(m2 == 0, 0, m3 / m2**1.5)
    finally:
        np.seterr(**olderr)

    if not bias:
        can_correct = (n > 2) & (m2 > 0)
        if can_correct.any():
            m2 = np.extract(can_correct, m2)
            m3 = np.extract(can_correct, m3)
            nval = ma.sqrt((n-1.0)*n)/(n-2.0)*m3/m2**1.5
            np.place(vals, can_correct, nval)
    return vals
skew.__doc__ = stats.skew.__doc__


def kurtosis(a, axis=0, fisher=True, bias=True):
    a, axis = _chk_asarray(a, axis)
    m2 = moment(a,2,axis)
    m4 = moment(a,4,axis)
    olderr = np.seterr(all='ignore')
    try:
        vals = ma.where(m2 == 0, 0, m4 / m2**2.0)
    finally:
        np.seterr(**olderr)

    if not bias:
        n = a.count(axis)
        can_correct = (n > 3) & (m2 is not ma.masked and m2 > 0)
        if can_correct.any():
            n = np.extract(can_correct, n)
            m2 = np.extract(can_correct, m2)
            m4 = np.extract(can_correct, m4)
            nval = 1.0/(n-2)/(n-3)*((n*n-1.0)*m4/m2**2.0-3*(n-1)**2.0)
            np.place(vals, can_correct, nval+3.0)
    if fisher:
        return vals - 3
    else:
        return vals
kurtosis.__doc__ = stats.kurtosis.__doc__


def describe(a, axis=0):
    """
    Computes several descriptive statistics of the passed array.

    Parameters
    ----------
    a : array

    axis : int or None

    Returns
    -------
    n : int
        (size of the data (discarding missing values)
    mm : (int, int)
        min, max

    arithmetic mean : float

    unbiased variance : float

    biased skewness : float

    biased kurtosis : float

    Examples
    --------

    >>> ma = np.ma.array(range(6), mask=[0, 0, 0, 1, 1, 1])
    >>> describe(ma)
    (array(3),
     (0, 2),
     1.0,
     1.0,
     masked_array(data = 0.0,
                 mask = False,
           fill_value = 1e+20)
    ,
     -1.5)

    """
    a, axis = _chk_asarray(a, axis)
    n = a.count(axis)
    mm = (ma.minimum.reduce(a), ma.maximum.reduce(a))
    m = a.mean(axis)
    v = a.var(axis)
    sk = skew(a, axis)
    kurt = kurtosis(a, axis)
    return n, mm, m, v, sk, kurt

#.............................................................................
def stde_median(data, axis=None):
    """Returns the McKean-Schrader estimate of the standard error of the sample
median along the given axis. masked values are discarded.

    Parameters
    ----------
        data : ndarray
            Data to trim.
        axis : {None,int}, optional
            Axis along which to perform the trimming.
            If None, the input array is first flattened.

    """
    def _stdemed_1D(data):
        data = np.sort(data.compressed())
        n = len(data)
        z = 2.5758293035489004
        k = int(np.round((n+1)/2. - z * np.sqrt(n/4.),0))
        return ((data[n-k] - data[k-1])/(2.*z))
    #
    data = ma.array(data, copy=False, subok=True)
    if (axis is None):
        return _stdemed_1D(data)
    else:
        if data.ndim > 2:
            raise ValueError("Array 'data' must be at most two dimensional, but got data.ndim = %d" % data.ndim)
        return ma.apply_along_axis(_stdemed_1D, axis, data)

#####--------------------------------------------------------------------------
#---- --- Normality Tests ---
#####--------------------------------------------------------------------------

def skewtest(a, axis=0):
    a, axis = _chk_asarray(a, axis)
    if axis is None:
        a = a.ravel()
        axis = 0
    b2 = skew(a,axis)
    n = a.count(axis)
    if np.min(n) < 8:
        warnings.warn(
            "skewtest only valid for n>=8 ... continuing anyway, n=%i" %
            np.min(n))
    y = b2 * ma.sqrt(((n+1)*(n+3)) / (6.0*(n-2)) )
    beta2 = ( 3.0*(n*n+27*n-70)*(n+1)*(n+3) ) / ( (n-2.0)*(n+5)*(n+7)*(n+9) )
    W2 = -1 + ma.sqrt(2*(beta2-1))
    delta = 1/ma.sqrt(0.5*ma.log(W2))
    alpha = ma.sqrt(2.0/(W2-1))
    y = ma.where(y==0, 1, y)
    Z = delta*ma.log(y/alpha + ma.sqrt((y/alpha)**2+1))
    return Z, (1.0 - stats.zprob(Z))*2
skewtest.__doc__ = stats.skewtest.__doc__

def kurtosistest(a, axis=0):
    a, axis = _chk_asarray(a, axis)
    n = a.count(axis=axis).astype(float)
    if np.min(n) < 20:
        warnings.warn(
            "kurtosistest only valid for n>=20 ... continuing anyway, n=%i" %
            np.min(n))
    b2 = kurtosis(a, axis, fisher=False)
    E = 3.0*(n-1) /(n+1)
    varb2 = 24.0*n*(n-2)*(n-3) / ((n+1)*(n+1)*(n+3)*(n+5))
    x = (b2-E)/ma.sqrt(varb2)
    sqrtbeta1 = 6.0*(n*n-5*n+2)/((n+7)*(n+9)) * np.sqrt((6.0*(n+3)*(n+5))/
                                                        (n*(n-2)*(n-3)))
    A = 6.0 + 8.0/sqrtbeta1 *(2.0/sqrtbeta1 + np.sqrt(1+4.0/(sqrtbeta1**2)))
    term1 = 1 - 2./(9.0*A)
    denom = 1 + x*ma.sqrt(2/(A-4.0))
    denom[denom < 0] = masked
    term2 = ma.power((1-2.0/A)/denom,1/3.0)
    Z = ( term1 - term2 ) / np.sqrt(2/(9.0*A))
    return Z, (1.0-stats.zprob(Z))*2
kurtosistest.__doc__ = stats.kurtosistest.__doc__


def normaltest(a, axis=0):
    a, axis = _chk_asarray(a, axis)
    s,_ = skewtest(a,axis)
    k,_ = kurtosistest(a,axis)
    k2 = s*s + k*k
    return k2, stats.chisqprob(k2,2)
normaltest.__doc__ = stats.normaltest.__doc__

# Martinez-Iglewicz test
# K-S test


#####--------------------------------------------------------------------------
#---- --- Percentiles ---
#####--------------------------------------------------------------------------


def mquantiles(a, prob=list([.25,.5,.75]), alphap=.4, betap=.4, axis=None,
               limit=()):
    """
    Computes empirical quantiles for a data array.

    Samples quantile are defined by ``Q(p) = (1-gamma)*x[j] + gamma*x[j+1]``,
    where ``x[j]`` is the j-th order statistic, and gamma is a function of
    ``j = floor(n*p + m)``, ``m = alphap + p*(1 - alphap - betap)`` and
    ``g = n*p + m - j``.

    Reinterpreting the above equations to compare to **R** lead to the
    equation: ``p(k) = (k - alphap)/(n + 1 - alphap - betap)``

    Typical values of (alphap,betap) are:
        - (0,1)    : ``p(k) = k/n`` : linear interpolation of cdf
          (**R** type 4)
        - (.5,.5)  : ``p(k) = (k - 1/2.)/n`` : piecewise linear function
          (**R** type 5)
        - (0,0)    : ``p(k) = k/(n+1)`` :
          (**R** type 6)
        - (1,1)    : ``p(k) = (k-1)/(n-1)``: p(k) = mode[F(x[k])].
          (**R** type 7, **R** default)
        - (1/3,1/3): ``p(k) = (k-1/3)/(n+1/3)``: Then p(k) ~ median[F(x[k])].
          The resulting quantile estimates are approximately median-unbiased
          regardless of the distribution of x.
          (**R** type 8)
        - (3/8,3/8): ``p(k) = (k-3/8)/(n+1/4)``: Blom.
          The resulting quantile estimates are approximately unbiased
          if x is normally distributed
          (**R** type 9)
        - (.4,.4)  : approximately quantile unbiased (Cunnane)
        - (.35,.35): APL, used with PWM

    Parameters
    ----------
    a : array_like
        Input data, as a sequence or array of dimension at most 2.
    prob : array_like, optional
        List of quantiles to compute.
    alphap : float, optional
        Plotting positions parameter, default is 0.4.
    betap : float, optional
        Plotting positions parameter, default is 0.4.
    axis : int, optional
        Axis along which to perform the trimming.
        If None (default), the input array is first flattened.
    limit : tuple
        Tuple of (lower, upper) values.
        Values of `a` outside this open interval are ignored.

    Returns
    -------
    mquantiles : MaskedArray
        An array containing the calculated quantiles.

    Notes
    -----
    This formulation is very similar to **R** except the calculation of
    ``m`` from ``alphap`` and ``betap``, where in **R** ``m`` is defined
    with each type.

    References
    ----------
    .. [1] *R* statistical software at http://www.r-project.org/

    Examples
    --------
    >>> from scipy.stats.mstats import mquantiles
    >>> a = np.array([6., 47., 49., 15., 42., 41., 7., 39., 43., 40., 36.])
    >>> mquantiles(a)
    array([ 19.2,  40. ,  42.8])

    Using a 2D array, specifying axis and limit.

    >>> data = np.array([[   6.,    7.,    1.],
                         [  47.,   15.,    2.],
                         [  49.,   36.,    3.],
                         [  15.,   39.,    4.],
                         [  42.,   40., -999.],
                         [  41.,   41., -999.],
                         [   7., -999., -999.],
                         [  39., -999., -999.],
                         [  43., -999., -999.],
                         [  40., -999., -999.],
                         [  36., -999., -999.]])
    >>> mquantiles(data, axis=0, limit=(0, 50))
    array([[ 19.2 ,  14.6 ,   1.45],
           [ 40.  ,  37.5 ,   2.5 ],
           [ 42.8 ,  40.05,   3.55]])

    >>> data[:, 2] = -999.
    >>> mquantiles(data, axis=0, limit=(0, 50))
    masked_array(data =
     [[19.2 14.6 --]
     [40.0 37.5 --]
     [42.8 40.05 --]],
                 mask =
     [[False False  True]
      [False False  True]
      [False False  True]],
           fill_value = 1e+20)

    """
    def _quantiles1D(data,m,p):
        x = np.sort(data.compressed())
        n = len(x)
        if n == 0:
            return ma.array(np.empty(len(p), dtype=float), mask=True)
        elif n == 1:
            return ma.array(np.resize(x, p.shape), mask=nomask)
        aleph = (n*p + m)
        k = np.floor(aleph.clip(1, n-1)).astype(int)
        gamma = (aleph-k).clip(0,1)
        return (1.-gamma)*x[(k-1).tolist()] + gamma*x[k.tolist()]

    # Initialization & checks ---------
    data = ma.array(a, copy=False)
    if data.ndim > 2:
        raise TypeError("Array should be 2D at most !")
    #
    if limit:
        condition = (limit[0]<data) & (data<limit[1])
        data[~condition.filled(True)] = masked
    #
    p = np.array(prob, copy=False, ndmin=1)
    m = alphap + p*(1.-alphap-betap)
    # Computes quantiles along axis (or globally)
    if (axis is None):
        return _quantiles1D(data, m, p)
    return ma.apply_along_axis(_quantiles1D, axis, data, m, p)


def scoreatpercentile(data, per, limit=(), alphap=.4, betap=.4):
    """Calculate the score at the given 'per' percentile of the
    sequence a.  For example, the score at per=50 is the median.

    This function is a shortcut to mquantile

    """
    if (per < 0) or (per > 100.):
        raise ValueError("The percentile should be between 0. and 100. !"\
                         " (got %s)" % per)
    return mquantiles(data, prob=[per/100.], alphap=alphap, betap=betap,
                      limit=limit, axis=0).squeeze()


def plotting_positions(data, alpha=0.4, beta=0.4):
    """
    Returns plotting positions (or empirical percentile points) for the data.

    Plotting positions are defined as ``(i-alpha)/(n+1-alpha-beta)``, where:
        - i is the rank order statistics
        - n is the number of unmasked values along the given axis
        - alpha and beta are two parameters.

    Typical values for alpha and beta are:
        - (0,1)    : ``p(k) = k/n``, linear interpolation of cdf (R, type 4)
        - (.5,.5)  : ``p(k) = (k-1/2.)/n``, piecewise linear function
          (R, type 5)
        - (0,0)    : ``p(k) = k/(n+1)``, Weibull (R type 6)
        - (1,1)    : ``p(k) = (k-1)/(n-1)``, in this case,
          ``p(k) = mode[F(x[k])]``. That's R default (R type 7)
        - (1/3,1/3): ``p(k) = (k-1/3)/(n+1/3)``, then
          ``p(k) ~ median[F(x[k])]``.
          The resulting quantile estimates are approximately median-unbiased
          regardless of the distribution of x. (R type 8)
        - (3/8,3/8): ``p(k) = (k-3/8)/(n+1/4)``, Blom.
          The resulting quantile estimates are approximately unbiased
          if x is normally distributed (R type 9)
        - (.4,.4)  : approximately quantile unbiased (Cunnane)
        - (.35,.35): APL, used with PWM
        - (.3175, .3175): used in scipy.stats.probplot

    Parameters
    ----------
    data : array_like
        Input data, as a sequence or array of dimension at most 2.
    alpha : float, optional
        Plotting positions parameter. Default is 0.4.
    beta : float, optional
        Plotting positions parameter. Default is 0.4.

    Returns
    -------
    positions : MaskedArray
        The calculated plotting positions.

    """
    data = ma.array(data, copy=False).reshape(1,-1)
    n = data.count()
    plpos = np.empty(data.size, dtype=float)
    plpos[n:] = 0
    plpos[data.argsort()[:n]] = (np.arange(1, n+1) - alpha) / \
                                (n + 1.0 - alpha - beta)
    return ma.array(plpos, mask=data._mask)

meppf = plotting_positions

#####--------------------------------------------------------------------------
#---- --- Variability ---
#####--------------------------------------------------------------------------

def obrientransform(*args):
    """
Computes a transform on input data (any number of columns).  Used to
test for homogeneity of variance prior to running one-way stats.  Each
array in *args is one level of a factor.  If an F_oneway() run on the
transformed data and found significant, variances are unequal.   From
Maxwell and Delaney, p.112.

Returns: transformed data for use in an ANOVA
    """
    data = argstoarray(*args).T
    v = data.var(axis=0,ddof=1)
    m = data.mean(0)
    n = data.count(0).astype(float)
    # result = ((N-1.5)*N*(a-m)**2 - 0.5*v*(n-1))/((n-1)*(n-2))
    data -= m
    data **= 2
    data *= (n-1.5)*n
    data -= 0.5*v*(n-1)
    data /= (n-1.)*(n-2.)
    if not ma.allclose(v,data.mean(0)):
        raise ValueError("Lack of convergence in obrientransform.")
    return data


def signaltonoise(data, axis=0):
    """Calculates the signal-to-noise ratio, as the ratio of the mean over
    standard deviation along the given axis.

    Parameters
    ----------
        data : sequence
            Input data
        axis : {0, int}, optional
            Axis along which to compute. If None, the computation is performed
            on a flat version of the array.
"""
    data = ma.array(data, copy=False)
    m = data.mean(axis)
    sd = data.std(axis, ddof=0)
    return m/sd


def sem(a, axis=0):
    a, axis = _chk_asarray(a, axis)
    n = a.count(axis=axis)
    s = a.std(axis=axis,ddof=0) / ma.sqrt(n-1)
    return s
sem.__doc__ = stats.sem.__doc__

zmap = stats.zmap
zscore = stats.zscore


#####--------------------------------------------------------------------------
#---- --- ANOVA ---
#####--------------------------------------------------------------------------


def f_oneway(*args):
    """
Performs a 1-way ANOVA, returning an F-value and probability given
any number of groups.  From Heiman, pp.394-7.

Usage:   f_oneway (*args)    where *args is 2 or more arrays, one per
                                  treatment group
Returns: f-value, probability
"""
    # Construct a single array of arguments: each row is a group
    data = argstoarray(*args)
    ngroups = len(data)
    ntot = data.count()
    sstot = (data**2).sum() - (data.sum())**2/float(ntot)
    ssbg = (data.count(-1) * (data.mean(-1)-data.mean())**2).sum()
    sswg = sstot-ssbg
    dfbg = ngroups-1
    dfwg = ntot - ngroups
    msb = ssbg/float(dfbg)
    msw = sswg/float(dfwg)
    f = msb/msw
    prob = stats.fprob(dfbg,dfwg,f)
    return f, prob


def f_value_wilks_lambda(ER, EF, dfnum, dfden, a, b):
    """Calculation of Wilks lambda F-statistic for multivarite data, per
    Maxwell & Delaney p.657.
    """
    ER = ma.array(ER, copy=False, ndmin=2)
    EF = ma.array(EF, copy=False, ndmin=2)
    if ma.getmask(ER).any() or ma.getmask(EF).any():
        raise NotImplementedError("Not implemented when the inputs "\
                                  "have missing data")
    lmbda = np.linalg.det(EF) / np.linalg.det(ER)
    q = ma.sqrt( ((a-1)**2*(b-1)**2 - 2) / ((a-1)**2 + (b-1)**2 -5) )
    q = ma.filled(q, 1)
    n_um = (1 - lmbda**(1.0/q))*(a-1)*(b-1)
    d_en = lmbda**(1.0/q) / (n_um*q - 0.5*(a-1)*(b-1) + 1)
    return n_um / d_en



def friedmanchisquare(*args):
    """Friedman Chi-Square is a non-parametric, one-way within-subjects ANOVA.
    This function calculates the Friedman Chi-square test for repeated measures
    and returns the result, along with the associated probability value.

    Each input is considered a given group. Ideally, the number of treatments
    among each group should be equal. If this is not the case, only the first
    n treatments are taken into account, where n is the number of treatments
    of the smallest group.
    If a group has some missing values, the corresponding treatments are masked
    in the other groups.
    The test statistic is corrected for ties.

    Masked values in one group are propagated to the other groups.

    Returns: chi-square statistic, associated p-value
    """
    data = argstoarray(*args).astype(float)
    k = len(data)
    if k < 3:
        raise ValueError("Less than 3 groups (%i): " % k +\
                         "the Friedman test is NOT appropriate.")
    ranked = ma.masked_values(rankdata(data, axis=0), 0)
    if ranked._mask is not nomask:
        ranked = ma.mask_cols(ranked)
        ranked = ranked.compressed().reshape(k,-1).view(ndarray)
    else:
        ranked = ranked._data
    (k,n) = ranked.shape
    # Ties correction
    repeats = np.array([find_repeats(_) for _ in ranked.T], dtype=object)
    ties = repeats[repeats.nonzero()].reshape(-1,2)[:,-1].astype(int)
    tie_correction = 1 - (ties**3-ties).sum()/float(n*(k**3-k))
    #
    ssbg = np.sum((ranked.sum(-1) - n*(k+1)/2.)**2)
    chisq = ssbg * 12./(n*k*(k+1)) * 1./tie_correction
    return chisq, stats.chisqprob(chisq,k-1)

#-############################################################################-#
