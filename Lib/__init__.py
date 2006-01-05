"""\
SciPy --- A scientific computing package for Python
===================================================

You can support the development of SciPy by purchasing documentation
at

  http://www.trelgol.com

It is being distributed for a fee for a limited time to try and raise
money for development.

Documentation is also available in the docstrings.

Available subpackages
---------------------
"""

import os, sys
SCIPY_IMPORT_VERBOSE = int(os.environ.get('SCIPY_IMPORT_VERBOSE','0'))

try:
    import pkg_resources # activate namespace packages (manipulates __path__)
except ImportError:
    pass

import numpy._import_tools as _ni
pkgload = _ni.PackageLoader()
del _ni

from numpy.testing import ScipyTest
test = ScipyTest('scipy').test
__all__.append('test')

from version import version as __version__
from numpy import __version__ as __numpy_version__
__all__.append('__version__')
__all__.append('__numpy_version__')

from __config__ import show as show_config
pkgload(verbose=SCIPY_IMPORT_VERBOSE,postpone=True)
