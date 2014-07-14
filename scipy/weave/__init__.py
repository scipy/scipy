"""
C/C++ integration
=================

NOTE: this module is deprecated and will be removed from Scipy before
the 1.0 release -- use the standalone weave package
(`https://github.com/scipy/weave`_) instead.

        inline     -- a function for including C/C++ code within Python
        blitz      -- a function for compiling Numeric expressions to C++
        ext_tools  -- a module that helps construct C/C++ extension modules.
        accelerate -- a module that inline accelerates Python functions


.. note:: On Linux one needs to have the Python development headers installed
          in order to be able to compile things with the `weave` module.
          Since this is a runtime dependency these headers (typically in a
          pythonX.Y-dev package) are not always installed when installing
          scipy.

"""
from __future__ import absolute_import, print_function

import sys

from numpy import deprecate


if sys.version_info[0] >= 3:
    raise ImportError("scipy.weave only supports Python 2.x")


@deprecate(old_name="scipy.weave", new_name="weave")
def _deprecated():
    pass
try:
    _deprecated()
except DeprecationWarning as e:
    # don't fail import if DeprecationWarnings raise error -- works around
    # the situation with Numpy's test framework
    pass


from .weave_version import weave_version as __version__

try:
    from .blitz_tools import blitz, BlitzWarning
except ImportError:
    pass  # scipy (core) wasn't available

from .inline_tools import inline
from . import ext_tools
from .ext_tools import ext_module, ext_function
try:
    from .accelerate_tools import accelerate
except:
    pass

from numpy.testing import Tester
test = Tester().test
