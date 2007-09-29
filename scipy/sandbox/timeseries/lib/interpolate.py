"""
A collection of intperolation tools for timeseries

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import numpy.core.numeric as numeric

from scipy.interpolate import fitpack

import maskedarray as MA
from maskedarray.core import masked, nomask, getmask
from maskedarray.extras import flatnotmasked_edges
marray = MA.array

__all__ = [
    'forward_fill', 'backward_fill', 'interp_masked1d',
          ]

#####---------------------------------------------------------------------------
#---- --- Functions for filling in masked values in a masked array ---
#####---------------------------------------------------------------------------
def forward_fill(marr, maxgap=None):
    """forward_fill(marr, maxgap=None)

Forward fills masked values in a 1-d array when there are less maxgap
consecutive masked values. If maxgap is None, then forward fill all masked
values.
"""
    # Initialization ..................
    if numeric.ndim(marr) > 1:
        raise ValueError,"The input array should be 1D only!"
    marr = marray(marr, copy=True)
    if getmask(marr) is nomask or marr.size == 0:
        return marr
    #
    currGap = 0
    if maxgap is not None:
        for i in range(1, marr.size):
            if marr._mask[i]:
                currGap += 1
                if currGap <= maxgap and not marr._mask[i-1]:
                    marr._data[i] = marr._data[i-1]
                    marr._mask[i] = False
                elif currGap == maxgap + 1:
                    marr._mask[i-maxgap:i] = True
            else:
                currGap = 0               
    else:
        for i in range(1, marr.size):
            # CHECK: We should probable be able to speed things up here
            if marr._mask[i] and not marr._mask[i-1]:
                marr._data[i] = marr._data[i-1]
                marr._mask[i] = False
    return marr
#.............................................................................
def backward_fill(marr, maxgap=None):
    """backward_fill(marr, maxgap=None)

Backward fills masked values in a 1-d array when there are less than maxgap
consecutive masked values. If maxgap is None, then backward fills all
masked values.
"""
    return forward_fill(marr[::-1], maxgap=maxgap)[::-1]
#.............................................................................
def interp_masked1d(marr, kind='linear'):
    """interp_masked1d(marr, king='linear')

Interpolates masked values in marr according to method kind.
kind must be one of 'constant', 'linear', 'cubic', quintic'
"""
    if numeric.ndim(marr) > 1: 
        raise ValueError("array must be 1 dimensional!")
    #
    marr = marray(marr, copy=True)
    if getmask(marr) is nomask: 
        return marr
    #
    unmaskedIndices = (~marr._mask).nonzero()[0]
    if unmaskedIndices.size < 2: 
        return marr
    #    
    kind = kind.lower()
    if kind == 'constant': 
        return forward_fill(marr)
    try:
        k = {'linear' : 1,
             'cubic' : 3,
             'quintic' : 5}[kind.lower()]
    except KeyError:
        raise ValueError("Unsupported interpolation type.")
    
    first_unmasked, last_unmasked = flatnotmasked_edges(marr)
    
    vals = marr.data[unmaskedIndices]
    
    tck = fitpack.splrep(unmaskedIndices, vals, k=k)
    
    maskedIndices = marr._mask.nonzero()[0]
    interpIndices = maskedIndices[(maskedIndices > first_unmasked) & \
                                  (maskedIndices < last_unmasked)]
    marr[interpIndices] = fitpack.splev(interpIndices, tck).astype(marr.dtype)
    return marr
