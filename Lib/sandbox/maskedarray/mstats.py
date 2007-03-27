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
from numpy import bool_, float_, int_
from numpy import array as narray
from numpy.core import numeric as numeric

import maskedarray as MA
from maskedarray.core import masked, nomask, MaskedArray
from maskedarray.core import masked_array as marray
from maskedarray.extras import apply_along_axis


def _quantiles_1D(data,m,p):
    """Returns quantiles for 1D compressed data. 
    Used internally by `mquantiles`.
    
:Parameters:
    data : ndarray
        Array to quantize
    m : Sequence
    p : float ndarray
        Quantiles to compute
    """
    n = data.count()
    if n == 0:
        return MA.resize(masked, len(p))
    elif n == 1:
        return MA.resize(data,len(p))
    data = data.compressed()
    aleph = (n*p + m)
    k = numpy.floor(aleph.clip(1, n-1)).astype(int_)
    gamma = (aleph-k).clip(0,1)
    y = MA.sort(data)
    return (1.-gamma)*y[(k-1).tolist()] + gamma*y[k.tolist()]

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

    # Initialization & checks ---------
    data = marray(data, copy=False)
    assert data.ndim <= 2, "Array should be 2D at most !"
    p = narray(prob, copy=False, ndmin=1)
    m = alphap + p*(1.-alphap-betap)
    # Computes quantiles along axis (or globally)
    if (axis is None): 
        return _quantiles_1D(data, m, p)
    else:
        return apply_along_axis(_quantiles_1D, axis, data, m, p)
    

def _median1d(data):
    """Returns the median of a 1D masked array. Used internally by mmedian."""
    datac = data.compressed()
    if datac.size > 0:
        return numpy.median(data.compressed())
    return masked

def _median2d_1(data):
    data = marray(data, subok=True, copy=True)
    if data._mask is nomask:
        return numpy.median(data)
    if data.ndim != 2 :
        raise ValueError("Input array should be 2D!")
    (n,p) = data.shape
    if p < n//3:
        return apply_along_axis(_median1d, 0, data)
    data.sort(axis=0)
    counts = data.count(axis=0)
    midx = (counts//2)
    midx_even = (counts%2==0)
    med = marray(numeric.empty((data.shape[-1],), dtype=data.dtype))
    med[midx_even] = (data[midx-1]+data[midx])/2.
    med[numpy.logical_not(midx_even)] = data[midx]
    if not med._mask.any():
        med._mask = nomask
    return med
             
def _median2d_2(data):
    return apply_along_axis(_median1d, 0, data)
    

def mmedian(data):
    """Returns the median of data along the first axis. Missing data are discarded."""
    data = marray(data, subok=True, copy=True)
    if data._mask is nomask:
        return numpy.median(data)
    if data.ndim == 1:
        return _median1d(data)
#    elif data.ndim == 2:
#        (n, p) = data.shape
#        if p < n//3:
#            return apply_along_axis(_median1d, 0, data)
#        data.sort(axis=0)
#        counts = data.count(axis=0)
#        midx = (counts//2)
#        midx_even = (counts%2==0)
#        med = marray(numeric.empty((p,), dtype=data.dtype))
#        med[midx_even] = (data[midx-1]+data[midx])/2.
#        med[numpy.logical_not(midx_even)] = data[midx]
#        if not med._mask.any():
#            med._mask = nomask
#        return med
    return apply_along_axis(_median1d, 0, data)
    


################################################################################
if __name__ == '__main__':
    from maskedarray.testutils import assert_almost_equal, assert_equal
    import timeit
    import maskedarray
    
    if 1:
        (n,p) = (101,30)
        x = marray(numpy.linspace(-1.,1.,n),)
        x[:10] = x[-10:] = masked
        z = marray(numpy.empty((n,p), dtype=numpy.float_))
        z[:,0] = x[:]
        idx = numpy.arange(len(x))
        for i in range(1,p):
            numpy.random.shuffle(idx)
            z[:,i] = x[idx]
    
        assert_equal(mmedian(z[:,0]), 0)
        assert_equal(mmedian(z), numpy.zeros((p,)))
        
        x = maskedarray.arange(24).reshape(3,4,2)
        x[x%3==0] = masked
        assert_equal(mmedian(x), [[12,9],[6,15],[12,9],[18,15]])
        x.shape = (4,3,2)
        assert_equal(mmedian(x),[[99,10],[11,99],[13,14]])
        x = maskedarray.arange(24).reshape(4,3,2)
        x[x%5==0] = masked
        assert_equal(mmedian(x), [[12,10],[8,9],[16,17]])
        
        
        """  [[[0 1],  [2  3],  [4 5]]
              [[6 7],  [8  9],  [10 11]]
              [[9 13]  [14 15]  [16 17]]
             [[18 19]  [20 21]  [22 23]]],

 [[[-- 1]  [2 --]  [4 5]   [-- 7]]
  [[8 --]  [10 11] [-- 13] [14 --]]
 [[16 17]  [-- 19] [20 --] [22 23]]],

        """
    
    if 0:
        print "GO!"
        med_setup1 = "from __main__ import _median2d_1,z"
        med_setup3 = "from __main__ import mmedian,z"
        med_setup2 = "from __main__ import _median2d_2,z"
        (nrep, nloop) = (3,10)
        med_r1 = timeit.Timer("_median2d_1(z)", med_setup1).repeat(nrep,nloop)
        med_r2 = timeit.Timer("_median2d_2(z)", med_setup2).repeat(nrep,nloop)
        med_r3 = timeit.Timer("mmedian(z)", med_setup3).repeat(nrep,nloop)
        med_r1 = numpy.sort(med_r1)
        med_r2 = numpy.sort(med_r2)
        med_r3 = numpy.sort(med_r3)
        print "median2d_1 : %s" % med_r1
        print "median2d_2 : %s" % med_r2
        print "median     : %s" % med_r3