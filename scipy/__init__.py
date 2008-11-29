"""\
SciPy --- A scientific computing package for Python
===================================================

Documentation is available in the docstrings and
online at http://docs.scipy.org.
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

# Emit a warning if numpy is too old
majver, minver = [float(i) for i in _num.version.version.split('.')[:2]]
if majver < 1 or (majver == 1 and minver < 2):
    import warnings
    warnings.warn("Numpy 1.2.0 or above is recommended for this version of " \
                  "scipy (detected version %s)" % _num.version.version,
                  UserWarning)

__all__ += ['oldnumeric']+_num.__all__

__all__ += ['randn', 'rand', 'fft', 'ifft']

if __doc__:
    __doc__ += """
Contents
--------
SciPy imports all the functions from the NumPy namespace, and in
addition provides:"""

del _num
# Remove the linalg imported from numpy so that the scipy.linalg package can be
# imported.
del linalg
__all__.remove('linalg')

try:
    from __config__ import show as show_config
except ImportError, e:
    msg = """Error importing scipy: you cannot import scipy while
    being in scipy source directory; please exit the scipy source
    tree first, and relaunch your python intepreter."""
    raise ImportError(msg)
from version import version as __version__

# Load scipy packages, their global_symbols, set up __doc__ string.
from numpy._import_tools import PackageLoader
import os as _os
SCIPY_IMPORT_VERBOSE = int(_os.environ.get('SCIPY_IMPORT_VERBOSE','-1'))
del _os
pkgload = PackageLoader()
pkgload(verbose=SCIPY_IMPORT_VERBOSE,postpone=True)

if __doc__:
    __doc__ += """

Available subpackages
---------------------
"""
if __doc__:
    __doc__ += pkgload.get_pkgdocs()

from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
if __doc__:
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
