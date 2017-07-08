"""
Generic test utilities and decorators.

"""

from __future__ import division, print_function, absolute_import

import os
import sys

import warnings

from scipy._lib.decorator import decorator


__all__ = ['suppressed_stdout', 'xslow_yield']


class TestutilDeprecationWarning(DeprecationWarning):
    pass


def suppressed_stdout(f):
    import nose

    def pwrapper(*arg, **kwargs):
        warnings.warn("scipy._lib._testutils.suppressed_stdout is deprecated "
                      "-- should use pytest capture fixture instead",
                      category=TestutilDeprecationWarning, stacklevel=1)

        oldstdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        try:
            return f(*arg, **kwargs)
        finally:
            sys.stdout.close()
            sys.stdout = oldstdout
    return nose.tools.make_decorator(f)(pwrapper)


@decorator
def xslow_yield(func, *a, **kw):
    try:
        v = int(os.environ.get('SCIPY_XSLOW', '0'))
        if not v:
            raise ValueError()
    except ValueError:
        import pytest
        pytest.skip("very slow test; set environment variable "
                    "SCIPY_XSLOW=1 to run it")
    return func(*a, **kw)


def skipif_yield(condition, reason, msg=""):
    """
    Similar to pytest.mark.skipif, for use in yield tests.

    For yield tests, pytest.mark.skipif does not work as expected ---
    if any condition evaluates to true, *all* of the yielded tests are
    skipped.
    """
    @decorator
    def wrapper(func, *a, **kw):
        import pytest
        if condition:
            pytest.skip(reason)
        return func(*a, **kw)
    return wrapper


class PytestTester(object):
    """
    Pytest test runner entry point.
    """

    def __init__(self, module_name):
        self.module_name = module_name

    def __call__(self, label="fast", verbose=1, extra_argv=None, doctests=False,
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
