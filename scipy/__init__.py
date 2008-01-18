"""\
SciPy --- A scientific computing package for Python
===================================================

You can support the development of SciPy by purchasing documentation
at

  http://www.trelgol.com

It is being distributed for a fee for a limited time to try and raise
money for development.

Documentation is also available in the docstrings.

"""

__all__ = ['pkgload','test']

from numpy import show_config as show_numpy_config
if show_numpy_config is None:
    raise ImportError,"Cannot import scipy when running from numpy source directory."
from numpy import __version__ as __numpy_version__

# Import numpy symbols to scipy name space
import numpy as _num
from numpy import oldnumeric
from numpy import *
from numpy.random import rand, randn
from numpy.fft import fft, ifft
from numpy.lib.scimath import *
_num.seterr(all='ignore')

__all__ += ['oldnumeric']+_num.__all__

__all__ += ['randn', 'rand', 'fft', 'ifft']

__doc__ += """
Contents
--------

  numpy name space
"""
del _num
# Remove the linalg imported from numpy so that the scipy.linalg package can be
# imported.
del linalg

from __config__ import show as show_config
from version import version as __version__

# Load scipy packages, their global_symbols, set up __doc__ string.
from numpy._import_tools import PackageLoader
import os as _os
SCIPY_IMPORT_VERBOSE = int(_os.environ.get('SCIPY_IMPORT_VERBOSE','-1'))
del _os
pkgload = PackageLoader()
pkgload(verbose=SCIPY_IMPORT_VERBOSE,postpone=True)
__doc__ += """

Available subpackages
---------------------
"""
__doc__ += pkgload.get_pkgdocs()

from testing.pkgtester import Tester
test = Tester().test
bench = Tester().bench
__doc__ += """

Utility tools
-------------

  test        --- Run scipy unittests
  pkgload     --- Load scipy packages
  show_config --- Show scipy build configuration
  show_numpy_config --- Show numpy build configuration
  __version__ --- Scipy version string
  __numpy_version__ --- Numpy version string

Environment variables
---------------------

  SCIPY_IMPORT_VERBOSE --- pkgload verbose flag, default is 0.
"""
