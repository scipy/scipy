from .gmpy2 import *
# Internal variables/functions are not imported by * above.
# These are used by some python level functions and are needed
# at the top level.
# Use try...except to for static builds were _C_API is not available.
try:
    from .gmpy2 import _C_API, _mpmath_normalize, _mpmath_create
except ImportError:
    from .gmpy2 import _mpmath_normalize, _mpmath_create
