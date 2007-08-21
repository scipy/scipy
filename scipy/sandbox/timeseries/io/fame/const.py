import cfame
from cfame import FAME_CONSTANTS

"""add constants in cfame.FAME_CONSTANTS dictionary to global namespace
for this module"""

_g = globals()
for var, val in FAME_CONSTANTS.iteritems():
    _g[var] = val
    
__all__ = list(FAME_CONSTANTS)