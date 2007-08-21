"""Masked arrays add-ons.

A collection of utilities for maskedarray

:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import core
from core import *

import extras
from extras import *


__all__ = ['core', 'extras']
__all__ += core.__all__
__all__ += extras.__all__