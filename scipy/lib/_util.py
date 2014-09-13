import sys
import warnings
import os

from nose import SkipTest
from .decorator import decorator


class DeprecatedImport(object):
    """
    Deprecated import, with redirection + warning.

    Examples
    --------
    Suppose you previously had in some module::

        from foo import spam

    If this has to be deprecated, do::

        spam = DeprecatedImport("foo.spam", "baz")

    to redirect users to use "baz" module instead.

    """

    def __init__(self, old_module_name, new_module_name):
        self._old_name = old_module_name
        self._new_name = new_module_name
        __import__(self._new_name)
        self._mod = sys.modules[self._new_name]

    def __dir__(self):
        return dir(self._mod)

    def __getattr__(self, name):
        warnings.warn("Module %s is deprecated, use %s instead"
                      % (self._old_name, self._new_name),
                      DeprecationWarning)
        return getattr(self._mod, name)


def setup_xslow():
    """Convenience function for skipping xslow tests in nose generators.

    Attaching a custom decorator to a nose test generator makes nose
    unconditionally skip tests.
    Use this one together with nose.with_setup::

        from nose import with_setup

        @with_setup(setup_xslow)
        def test_gen():
            yield check_a, a
            ...

        def check_a(a):
            ...

    """
    try:
        v = int(os.environ.get('SCIPY_XSLOW', '0'))
        if not v:
            raise ValueError()
    except ValueError:
        raise SkipTest("very slow test; set environment variable "
                       "SCIPY_XSLOW=1 to run it")

