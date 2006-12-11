"""Masked arrays add-ons.

A collection of utilities for maskedarray

:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id: __init__.py 38 2006-12-09 23:01:14Z backtopop $
"""
__author__ = "Pierre GF Gerard-Marchant ($Author: backtopop $)"
__version__ = '1.0'
__revision__ = "$Revision: 38 $"
__date__     = '$Date: 2006-12-09 18:01:14 -0500 (Sat, 09 Dec 2006) $'

import core
reload(core)
from core import *

import extras
reload(extras)
from extras import *


__all__ = ['core', 'extras']
__all__ += core.__all__
__all__ += extras.__all__