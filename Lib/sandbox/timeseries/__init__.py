"""TimeSeries

Support for time series in numpy/scipy

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""


__author__ = "Pierre GF Gerard-Marchant  & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import tcore
from tcore import *
import tdates
from tdates import *
import tseries
from tseries import *
import tmulti
from tmulti import *
import reportlib
from reportlib import *
import plotlib
import addons


__all__ = ['tdates', 'tseries','tmulti','reportlib']
__all__ += tdates.__all__
__all__ += tseries.__all__
__all__ += tmulti.__all__
__all__ += reportlib.__all__

