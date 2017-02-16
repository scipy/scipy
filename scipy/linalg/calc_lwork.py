"""
calc_lwork was an internal module in Scipy and has been removed.

Several functions in scipy.linalg.lapack have *_lwork variants
that perform the lwork calculation (from Scipy >= 0.15.0), or
allow passing in LWORK=-1 argument to perform the computation.

"""

from __future__ import division, print_function, absolute_import

from numpy import deprecate

from ._calc_lwork import *

@deprecate(old_name="scipy.linalg.calc_lwork", message=__doc__)
def _deprecated():
    pass
try:
    _deprecated()
except DeprecationWarning as e:
    # don't fail import if DeprecationWarnings raise error -- works around
    # the situation with Numpy's test framework
    pass
