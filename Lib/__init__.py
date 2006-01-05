"""\
SciPy Core
==========

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
NUMPY_IMPORT_VERBOSE = int(os.environ.get('NUMPY_IMPORT_VERBOSE','0'))

try:
    import pkg_resources # activate namespace packages (manipulates __path__)
except ImportError:
    pass

import numpy._import_tools as ni
pkgload = ni.PackageLoader()
del ni

from numpy.testing import ScipyTest
test = ScipyTest('scipy').test
__all__.append('test')

__numpy_doc__ = """

NumPy: A scientific computing package for Python
================================================

Available subpackages
---------------------
"""

try:
    from __numpy_config__ import show as show_numpy_config
except ImportError:
    show_numpy_config = None


if show_numpy_config is not None:
    from numpy import __version__ as __numpy_version__
    __doc__ += __numpy_doc__
    pkgload(verbose=NUMPY_IMPORT_VERBOSE,postpone=True)
