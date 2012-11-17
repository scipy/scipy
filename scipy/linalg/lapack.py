#
# Author: Pearu Peterson, March 2002
#

__all__ = ['get_lapack_funcs']

# The following ensures that possibly missing flavor (C or Fortran) is
# replaced with the available one. If none is available, exception
# is raised at the first attempt to use the resources.
from blas import _get_funcs

from scipy.linalg import flapack
try:
    from scipy.linalg import clapack
except ImportError:
    clapack = None

_use_force_clapack = 1
if clapack is None:
    clapack = flapack
    _use_force_clapack = 0
elif hasattr(flapack,'empty_module'):
    flapack = clapack

_lapack_alias = {}

def get_lapack_funcs(names, arrays=(), dtype=None):
    """Return available LAPACK function objects from names.

    Arrays are used to determine the optimal prefix of LAPACK routines.

    Parameters
    ----------
    names : str or sequence of str
        Name(s) of LAPACK functions withouth type prefix.

    arrays : sequency of ndarrays, optional
        Arrays can be given to determine optiomal prefix of LAPACK
        routines. If not given, double-precision routines will be
        used, otherwise the most generic type in arrays will be used.

    dtype : str or dtype, optional
        Data-type specifier. Not used if `arrays` is non-empty.


    Returns
    -------
    funcs : list
        List containing the found function(s).


    Notes
    -----
    This routines automatically chooses between Fortran/C
    interfaces. Fortran code is used whenever possible for arrays with
    column major order. In all other cases, C code is preferred.

    In LAPACK, the naming convention is that all functions start with a
    type prefix, which depends on the type of the principal
    matrix. These can be one of {'s', 'd', 'c', 'z'} for the numpy
    types {float32, float64, complex64, complex128} respectevely, and
    are stored in attribute `typecode` of the returned functions.
    """
    return _get_funcs(names, arrays, dtype,
                      "LAPACK", flapack, clapack, _lapack_alias)
