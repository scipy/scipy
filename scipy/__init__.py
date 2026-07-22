"""
.. The heading is listed in the parent file `doc/reference/index.rst` to keep the
   section levels consistent.

The main ``scipy`` namespace has very few objects in it by design. Only some
generical functionality related to testing, build info, and versioning, and one
class (`LowLevelCallable`), which didn't fit into one of the submodules, is present:

.. autosummary::
   :toctree: generated/

   scipy.LowLevelCallable
   scipy.show_config
   scipy.test

The sole public attribute is:

================== ===============================================
``__version__``    SciPy version string
================== ===============================================
"""

import importlib as _importlib

from numpy import __version__ as __numpy_version__


try:
    from scipy.__config__ import show as show_config
except ImportError as e:
    msg = """Error importing SciPy: you cannot import SciPy while
    being in scipy source directory; please exit the SciPy source
    tree first and relaunch your Python interpreter."""
    raise ImportError(msg) from e


from scipy.version import version as __version__


# Allow distributors to run custom init code
from . import _distributor_init
del _distributor_init


from scipy._external.packaging_version.version import Version, parse
# In maintenance branch, change to np_maxversion N+3 if numpy is at N
np_minversion = '2.0.0'
np_maxversion = '9.9.99'
if (parse(__numpy_version__) < Version(np_minversion) or
        parse(__numpy_version__) >= Version(np_maxversion)):
    import warnings
    warnings.warn(f"A NumPy version >={np_minversion} and <{np_maxversion}"
                  f" is required for this version of SciPy (detected "
                  f"version {__numpy_version__})",
                  UserWarning, stacklevel=2)
del Version, parse


# This is the first import of an extension module within SciPy. If there's
# a general issue with the install, such that extension modules are missing
# or cannot be imported, this is where we'll get a failure - so give an
# informative error message.
try:
    from scipy._lib._ccallback import LowLevelCallable
except ImportError as e:
    msg = "The `scipy` install you are using seems to be broken, " + \
          "(extension modules cannot be imported), " + \
          "please try reinstalling."
    raise ImportError(msg) from e


from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester


submodules = [
    'cluster',
    'constants',
    'datasets',
    'differentiate',
    'fft',
    'fftpack',
    'integrate',
    'interpolate',
    'io',
    'linalg',
    'ndimage',
    'optimize',
    'signal',
    'sparse',
    'spatial',
    'special',
    'stats'
]

__all__ = submodules + [
    'LowLevelCallable',
    'test',
    'show_config',
    '__version__',
]


def __dir__():
    return __all__


def __getattr__(name):
    if name in submodules:
        return _importlib.import_module(f'scipy.{name}')
    elif name == "odr":
        raise AttributeError(
            "`scipy.odr` was deprecated in SciPy 1.17 and removed in SciPy 1.19. "
            "Please use https://pypi.org/project/odrpack/ instead."
        )
    else:
        try:
            return globals()[name]
        except KeyError:
            raise AttributeError(
                f"Module 'scipy' has no attribute '{name}'"
            )
