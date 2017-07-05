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


class PytestTester(object):
    """
    Pytest test runner entry point.
    """

    def __init__(self, module_name):
        self.module_name = module_name

    def __call__(self, label=None, verbose=1, extra_argv=None, doctests=False,
                 coverage=False, tests=None):
        import pytest

        module = sys.modules[self.module_name]
        module_path = os.path.abspath(module.__path__[0])

        pytest_args = ['-l']

        if doctests:
            raise ValueError("Doctests not supported")

        if extra_argv:
            pytest_args += list(extra_argv)

        if verbose and int(verbose) > 1:
            pytest_args += ["-" + "v"*(int(verbose)-1)]

        if coverage:
            pytest_args += ["--cov=" + module_path]

        if label == "fast":
            pytest_args += ["-m", "not slow"]
        elif label != "full":
            pytest_args += ["-m", label]

        if tests is None:
            tests = [self.module_name]

        pytest_args += ['--pyargs'] + list(tests)

        try:
            code = pytest.main(pytest_args)
        except SystemExit as exc:
            code = exc.code

        return (code == 0)
