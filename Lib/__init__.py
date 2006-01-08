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

try:
    import pkg_resources as _pr # activate namespace packages (manipulates __path__)
    del _pr
except ImportError:
    pass

from numpy import show_config as show_numpy_config
if show_numpy_config is None:
    raise ImportError,"Cannot import scipy when running from numpy source directory."
from numpy import __version__ as __numpy_version__

from __config__ import show as show_config
from version import version as __version__

import numpy._import_tools as _ni
pkgload = _ni.PackageLoader()
del _ni

import os as _os
SCIPY_IMPORT_VERBOSE = int(_os.environ.get('SCIPY_IMPORT_VERBOSE','0'))
pkgload(verbose=SCIPY_IMPORT_VERBOSE,postpone=True)
del _os

from numpy.testing import ScipyTest
test = ScipyTest('scipy').test
__all__.append('test')

import numpy as _num
from numpy import *
__all__ += _num.__all__
del _num

__doc__ += pkgload.get_pkgdocs()
