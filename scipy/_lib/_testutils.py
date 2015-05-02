"""
Generic test utilities and decorators.

"""

from __future__ import division, print_function, absolute_import

import os
from numpy.testing import dec


__all__ = ['knownfailure_overridable']


def knownfailure_overridable(msg=None):
    if not msg:
        msg = "Undiagnosed issues (corner cases, wrong comparison values, or otherwise)"
    msg = msg + " [Set environment variable SCIPY_XFAIL=1 to run this test nevertheless.]"

    def deco(func):
        try:
            if bool(os.environ['SCIPY_XFAIL']):
                return func
        except (ValueError, KeyError):
            pass
        return dec.knownfailureif(True, msg)(func)
    return deco
