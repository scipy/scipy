"""

A collection of filters for timeseries

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import numpy as N
from numpy import bool_, float_
narray = N.array

from scipy.signal import convolve, get_window

import maskedarray as MA
from maskedarray import MaskedArray, nomask, getmask, getmaskarray, masked
marray = MA.array

from moving_funcs import mov_average_expw, cmov_average, cmov_mean, \
                         cmov_window

__all__ = ['mov_average_expw'
           'cmov_average', 'cmov_mean', 'cmov_window']
