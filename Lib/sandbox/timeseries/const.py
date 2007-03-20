"""
A collection of tools for timeseries

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id: tcore.py 2836 2007-03-07 16:58:14Z mattknox_ca $
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author: mattknox_ca $)"
__version__ = '1.0'
__revision__ = "$Revision: 2836 $"
__date__     = '$Date: 2007-03-07 11:58:14 -0500 (Wed, 07 Mar 2007) $'

from cseries import freq_constants

"""add constants in cseries.freq_constants dictionary to global namespace
for this module"""

__all__ = [list(freq_constants)]

_g = globals()
_g.update(freq_constants)
