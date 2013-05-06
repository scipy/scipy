"""
C/C++ integration
=================

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

if sys.version_info[0] >= 3:
    raise ImportError("scipy.weave only supports Python 2.x")


from .weave_version import weave_version as __version__

try:
    from .blitz_tools import blitz
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
