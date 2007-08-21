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
from numpy import bool_, float_, int_, \
    sqrt
from numpy import array as narray
import numpy.core.numeric as numeric
from numpy.core.numeric import concatenate

import maskedarray as MA
from maskedarray.core import masked, nomask, MaskedArray, masked_array
from maskedarray.extras import apply_along_axis, dot

__all__ = ['cov','meppf','plotting_positions','meppf','mmedian','mquantiles',
           'stde_median','trim_tail','trim_both','trimmed_mean','trimmed_stde',
           'winsorize']

#####--------------------------------------------------------------------------
#---- -- Trimming ---
#####--------------------------------------------------------------------------

def winsorize(data, alpha=0.2):
    """Returns a Winsorized version of the input array: the (alpha/2.) lowest
    values are set to the (alpha/2.)th percentile, and the (alpha/2.) highest
    values are set to the (1-alpha/2.)th percentile 
    Masked values are skipped. The input array is first flattened.
    """
    data = masked_array(data, copy=False).ravel()
    idxsort = data.argsort()
    (nsize, ncounts) = (data.size, data.count())
    ntrim = int(alpha * ncounts)
    (xmin,xmax) = data[idxsort[[ntrim, ncounts-nsize-ntrim-1]]]
    return masked_array(numpy.clip(data, xmin, xmax), mask=data._mask) 

#..............................................................................  
def trim_both(data, proportiontocut=0.2, axis=None):
    """Trims the data by masking the int(trim*n) smallest and int(trim*n) largest 
    values of data along the given axis, where n is the number of unmasked values.
    
:Inputs: 
    data : MaskedArray
        Data to trim.
    trim : float *[0.2]*
        Percentage of trimming. If n is the number of unmasked values before trimming, 
        the number of values after trimming is (1-2*trim)*n.
    axis : integer *[None]*
        Axis along which to perform the trimming.
    """
    #...................
    def _trim_1D(data, trim):
        "Private function: return a trimmed 1D array."
        nsize = data.size
        ncounts = data.count()
        ntrim = int(trim * ncounts)
        idxsort = data.argsort()
        data[idxsort[:ntrim]] = masked
        data[idxsort[ncounts-nsize-ntrim:]] = masked
        return data
    #...................
    data = masked_array(data, copy=False, subok=True)
    data.unshare_mask()
    if (axis is None): 
        return _trim_1D(data.ravel(), proportiontocut)
    else:
        assert data.ndim <= 2, "Array should be 2D at most !"
        return apply_along_axis(_trim_1D, axis, data, proportiontocut)

#..............................................................................
def trim_tail(data, proportiontocut=0.2, tail='left', axis=None):
    """Trims the data by masking int(trim*n) values from ONE tail of the data
    along the given axis, where n is the number of unmasked values.
    
:Inputs: 
    data : MaskedArray
        Data to trim.
    trim : float *[0.2]*
        Percentage of trimming. If n is the number of unmasked values before trimming, 
        the number of values after trimming is (1-2*trim)*n.
    axis : integer *[None]*
        Axis along which to perform the trimming.
    """
    #...................
    def _trim_1D(data, trim, left):
        "Private function: return a trimmed 1D array."
        nsize = data.size
        ncounts = data.count()
        ntrim = int(trim * ncounts)
        idxsort = data.argsort()
        if left:
            data[idxsort[:ntrim]] = masked
        else:
            data[idxsort[ncounts-nsize-ntrim:]] = masked
        return data
    #...................
    data = masked_array(data, copy=False, subok=True)
    data.unshare_mask()
    #
    if not isinstance(tail, str):
        raise TypeError("The tail argument should be in ('left','right')")
    tail = tail.lower()[0]
    if tail == 'l':
        left = True
    elif tail == 'r':
        left=False
    else:
        raise ValueError("The tail argument should be in ('left','right')")
    #
    if (axis is None): 
        return _trim_1D(data.ravel(), proportiontocut, left)
    else:
        assert data.ndim <= 2, "Array should be 2D at most !"
        return apply_along_axis(_trim_1D, axis, data, proportiontocut, left)

#..............................................................................    
def trimmed_mean(data, proportiontocut=0.2, axis=None):
    """Returns the trimmed mean of the data along the given axis. Trimming is
    performed on both ends of the distribution.
    
:Inputs: 
    data : MaskedArray
        Data to trim.
    proportiontocut : float *[0.2]*
        Proportion of the data to cut from each side of the data . 
        As a result, (2*proportiontocut*n) values are actually trimmed.
    axis : integer *[None]*
        Axis along which to perform the trimming.    
    """
    return trim_both(data, proportiontocut=proportiontocut, axis=axis).mean(axis=axis)

#..............................................................................   
def trimmed_stde(data, proportiontocut=0.2, axis=None):
    """Returns the standard error of the trimmed mean for the input data, 
    along the given axis. Trimming is performed on both ends of the distribution.
    
:Inputs: 
    data : MaskedArray
        Data to trim.
    proportiontocut : float *[0.2]*
        Proportion of the data to cut from each side of the data . 
        As a result, (2*proportiontocut*n) values are actually trimmed.
    axis : integer *[None]*
        Axis along which to perform the trimming.  
    """
    #........................
    def _trimmed_stde_1D(data, trim=0.2):
        "Returns the standard error of the trimmed mean for a 1D input data."
        winsorized = winsorize(data)
        nsize = winsorized.count()
        winstd = winsorized.stdu()
        return winstd / ((1-2*trim) * numpy.sqrt(nsize))
    #........................
    data = masked_array(data, copy=False, subok=True)
    data.unshare_mask()
    if (axis is None): 
        return _trimmed_stde_1D(data.ravel(), proportiontocut)
    else:
        assert data.ndim <= 2, "Array should be 2D at most !"
        return apply_along_axis(_trimmed_stde_1D, axis, data, proportiontocut)

#.............................................................................
def stde_median(data, axis=None):
    """Returns the McKean-Schrader estimate of the standard error of the sample
    median along the given axis.
    """
    def _stdemed_1D(data):
        sorted = numpy.sort(data.compressed())
        n = len(sorted)
        z = 2.5758293035489004
        k = int(round((n+1)/2. - z * sqrt(n/4.),0))
        return ((sorted[n-k] - sorted[k-1])/(2.*z))
    #
    data = masked_array(data, copy=False, subok=True)
    if (axis is None):
        return _stdemed_1D(data)
    else:
        assert data.ndim <= 2, "Array should be 2D at most !"
        return apply_along_axis(_stdemed_1D, axis, data)


#####--------------------------------------------------------------------------
#---- --- Quantiles ---
#####--------------------------------------------------------------------------


def mquantiles(data, prob=list([.25,.5,.75]), alphap=.4, betap=.4, axis=None):
    """Computes empirical quantiles for a *1xN* data array.
Samples quantile are defined by:
*Q(p) = (1-g).x[i] +g.x[i+1]*
where *x[j]* is the jth order statistic, 
with *i = (floor(n*p+m))*, *m=alpha+p*(1-alpha-beta)* and *g = n*p + m - i)*.

Typical values of (alpha,beta) are:
    
    - (0,1)    : *p(k) = k/n* : linear interpolation of cdf (R, type 4)
    - (.5,.5)  : *p(k) = (k+1/2.)/n* : piecewise linear function (R, type 5)
    - (0,0)    : *p(k) = k/(n+1)* : (R type 6)
    - (1,1)    : *p(k) = (k-1)/(n-1)*. In this case, p(k) = mode[F(x[k])].
      That's R default (R type 7)
    - (1/3,1/3): *p(k) = (k-1/3)/(n+1/3)*. Then p(k) ~ median[F(x[k])].
      The resulting quantile estimates are approximately median-unbiased
      regardless of the distribution of x. (R type 8)
    - (3/8,3/8): *p(k) = (k-3/8)/(n+1/4)*. Blom.
      The resulting quantile estimates are approximately unbiased
      if x is normally distributed (R type 9)
    - (.4,.4)  : approximately quantile unbiased (Cunnane)
    - (.35,.35): APL, used with PWM

:Parameters:
    x : Sequence
        Input data, as a sequence or array of dimension at most 2.
    prob : Sequence *[(0.25, 0.5, 0.75)]*
        List of quantiles to compute.
    alpha : Float (*[0.4]*)
        Plotting positions parameter.
    beta : Float (*[0.4]*)
        Plotting positions parameter.
    axis : Integer *[None]*
        Axis along which to compute quantiles. If *None*, uses the whole 
        (flattened/compressed) dataset.
    """
    def _quantiles1D(data,m,p):
        x = numpy.sort(data.compressed())
        n = len(x)
        if n == 0:
            return masked_array(numpy.empty(len(p), dtype=float_), mask=True)
        elif n == 1:
            return masked_array(numpy.resize(x, p.shape), mask=nomask)
        aleph = (n*p + m)
        k = numpy.floor(aleph.clip(1, n-1)).astype(int_)
        gamma = (aleph-k).clip(0,1)
        return (1.-gamma)*x[(k-1).tolist()] + gamma*x[k.tolist()]

    # Initialization & checks ---------
    data = masked_array(data, copy=False)
    p = narray(prob, copy=False, ndmin=1)
    m = alphap + p*(1.-alphap-betap)
    # Computes quantiles along axis (or globally)
    if (axis is None): 
        return _quantiles1D(data, m, p)
    else:
        assert data.ndim <= 2, "Array should be 2D at most !"
        return apply_along_axis(_quantiles1D, axis, data, m, p)
    

def plotting_positions(data, alpha=0.4, beta=0.4):
    """Returns the plotting positions (or empirical percentile points) for the
    data.
    Plotting positions are defined as (i-alpha)/(n-alpha-beta), where:
        - i is the rank order statistics
        - n is the number of unmasked values along the given axis
        - alpha and beta are two parameters.
    
    Typical values for alpha and beta are:     
        - (0,1)    : *p(k) = k/n* : linear interpolation of cdf (R, type 4)
        - (.5,.5)  : *p(k) = (k-1/2.)/n* : piecewise linear function (R, type 5)
        - (0,0)    : *p(k) = k/(n+1)* : Weibull (R type 6)
        - (1,1)    : *p(k) = (k-1)/(n-1)*. In this case, p(k) = mode[F(x[k])].
          That's R default (R type 7)
        - (1/3,1/3): *p(k) = (k-1/3)/(n+1/3)*. Then p(k) ~ median[F(x[k])].
          The resulting quantile estimates are approximately median-unbiased
          regardless of the distribution of x. (R type 8)
        - (3/8,3/8): *p(k) = (k-3/8)/(n+1/4)*. Blom.
          The resulting quantile estimates are approximately unbiased
          if x is normally distributed (R type 9)
        - (.4,.4)  : approximately quantile unbiased (Cunnane)
        - (.35,.35): APL, used with PWM
    """
    data = masked_array(data, copy=False).reshape(1,-1)
    n = data.count()
    plpos = numpy.empty(data.size, dtype=float_)
    plpos[n:] = 0
    plpos[data.argsort()[:n]] = (numpy.arange(1,n+1) - alpha)/(n+1-alpha-beta)
    return masked_array(plpos, mask=data._mask)

meppf = plotting_positions

 
def mmedian(data, axis=None):
    """Returns the median of data along the given axis. Missing data are discarded."""
    def _median1D(data):
        x = numpy.sort(data.compressed())
        if x.size == 0:
            return masked
        return numpy.median(x)
    data = masked_array(data, subok=True, copy=True)
    if axis is None:
        return _median1D(data)
    else:
        return apply_along_axis(_median1D, axis, data)
    
   
def cov(x, y=None, rowvar=True, bias=False, strict=False):
    """
    Estimate the covariance matrix.

    If x is a vector, return the variance.  For matrices, returns the covariance 
    matrix.

    If y is given, it is treated as an additional (set of) variable(s).

    Normalization is by (N-1) where N is the number of observations (unbiased 
    estimate).  If bias is True then normalization is by N.

    If rowvar is non-zero (default), then each row is a variable with observations 
    in the columns, otherwise each column is a variable  and the observations  are 
    in the rows.
    
    If strict is True, masked values are propagated: if a masked value appears in 
    a row or column, the whole row or column is considered masked.
    """
    X = narray(x, ndmin=2, subok=True, dtype=float)
    if X.shape[0] == 1:
        rowvar = True
    if rowvar:
        axis = 0
        tup = (slice(None),None)
    else:
        axis = 1
        tup = (None, slice(None))
    #
    if y is not None:
        y = narray(y, copy=False, ndmin=2, subok=True, dtype=float)
        X = concatenate((X,y),axis)
    #
    X -= X.mean(axis=1-axis)[tup]
    n = X.count(1-axis)
    #
    if bias:
        fact = n*1.0
    else:
        fact = n-1.0
    #
    if not rowvar:
        return (dot(X.T, X.conj(), strict=False) / fact).squeeze()
    else:
        return (dot(X, X.T.conj(), strict=False) / fact).squeeze()


def idealfourths(data, axis=None):
    """Returns an estimate of the interquartile range of the data along the given
    axis, as computed with the ideal fourths.
    """
    def _idf(data):
        x = numpy.sort(data.compressed())
        n = len(x)
        (j,h) = divmod(n/4. + 5/12.,1)
        qlo = (1-h)*x[j] + h*x[j+1]
        k = n - j
        qup = (1-h)*x[k] + h*x[k-1]
        return qup - qlo
    data = masked_array(data, copy=False)
    if (axis is None): 
        return _idf(data)
    else:
        return apply_along_axis(_idf, axis, data) 
    
    
def rsh(data, points=None):
    """Evalutates Rosenblatt's shifted histogram estimators for each
    point of 'points' on the dataset 'data'.
    
:Inputs:
    data : sequence
        Input data. Masked values are discarded.
    points : 
        Sequence of points where to evaluate Rosenblatt shifted histogram. 
        If None, use the data.        
    """
    data = masked_array(data, copy=False)
    if points is None:
        points = data
    else:
        points = numpy.array(points, copy=False, ndmin=1)
    if data.ndim != 1:
        raise AttributeError("The input array should be 1D only !")
    n = data.count()
    h = 1.2 * idealfourths(data) / n**(1./5)
    nhi = (data[:,None] <= points[None,:] + h).sum(0)
    nlo = (data[:,None] < points[None,:] - h).sum(0)
    return (nhi-nlo) / (2.*n*h)

    