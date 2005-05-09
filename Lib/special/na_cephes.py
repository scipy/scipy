'''This module contains GENERATED CODE which constructs numarray universal
functions from C-coded cfuncs contained in a sister extension module.

***************************** DO NOT EDIT **************************
'''

import _numarray._na_cephes
import numarray.ufunc as _uf

globals().update(_uf.make_ufuncs(_numarray._na_cephes))



_scipy_special_errprint = 0

def errprint(val=None):
    global _scipy_special_errprint
    old_val = _scipy_special_errprint
    if val is not None:
        _scipy_special_errprint = (val != 0)
    return old_val


