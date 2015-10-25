"""
Generic test utilities and decorators.

"""

from __future__ import division, print_function, absolute_import

import os
import sys
from numpy.testing import dec

from nose import SkipTest

from scipy._lib.decorator import decorator


__all__ = ['knownfailure_overridable', 'suppressed_stdout', 'xslow']


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


def suppressed_stdout(f):
    import nose

    def pwrapper(*arg, **kwargs):
        oldstdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        try:
            return f(*arg, **kwargs)
        finally:
            sys.stdout.close()
            sys.stdout = oldstdout
    return nose.tools.make_decorator(f)(pwrapper)


@decorator
def xslow(func, *a, **kw):
    try:
        v = int(os.environ.get('SCIPY_XSLOW', '0'))
        if not v:
            raise ValueError()
    except ValueError:
        raise SkipTest("very slow test; set environment variable "
                       "SCIPY_XSLOW=1 to run it")
    return func(*a, **kw)
