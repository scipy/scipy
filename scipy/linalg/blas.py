#
# Author: Pearu Peterson, March 2002
#         refactoring by Fabian Pedregosa, March 2010
#

__all__ = ['get_blas_funcs']

import numpy as np
# The following ensures that possibly missing flavor (C or Fortran) is
# replaced with the available one. If none is available, exception
# is raised at the first attempt to use the resources.

from scipy.linalg import cblas, fblas
if hasattr(cblas,'empty_module'):
    cblas = fblas
elif hasattr(fblas,'empty_module'):
    fblas = cblas

 # 'd' will be default for 'i',..
_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z', 'G':'z'}

# some convenience alias for complex functions
_blas_alias = {'cnrm2' : 'scnrm2', 'znrm2' : 'dznrm2',
               'cdot' : 'cdotc', 'zdot' : 'zdotc',
               'cger' : 'cgerc', 'zger' : 'zgerc',
               'sdotc': 'sdot', 'sdotu': 'sdot',
               'ddotc': 'ddot', 'ddotu': 'ddot'}


def get_blas_funcs(names, arrays=(), dtype=None):
    """Return available BLAS function objects from names.

    Arrays are used to determine the optimal prefix of BLAS routines.

    Parameters
    ----------
    names : str or sequence of str
        Name(s) of BLAS functions withouth type prefix.

    arrays : sequency of ndarrays, optional
        Arrays can be given to determine optiomal prefix of BLAS
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

    In BLAS, the naming convention is that all functions start with a
    type prefix, which depends on the type of the principal
    matrix. These can be one of {'s', 'd', 'c', 'z'} for the numpy
    types {float32, float64, complex64, complex128} respectevely, and
    are stored in attribute `typecode` of the returned functions.
    """

    blas_funcs = []
    unpack = False
    dtype = np.dtype(dtype)
    module1 = (cblas, 'cblas')
    module2 = (fblas, 'fblas')

    if isinstance(names, str):
        names = (names,)
        unpack = True

    if arrays:
        # use the most generic type in arrays
        dtype, index = max(
            [(ar.dtype, i) for i, ar in enumerate(arrays)])
        if arrays[index].flags['FORTRAN']:
            # prefer Fortran for leading array with column major order
            module1, module2 = module2, module1

    prefix = _type_conv.get(dtype.char, 'd')

    for i, name in enumerate(names):
        func_name = prefix + name
        func_name = _blas_alias.get(func_name, func_name)
        func = getattr(module1[0], func_name, None)
        module_name = module1[1]
        if func is None:
            func = getattr(module2[0], func_name, None)
            module_name = module2[1]
        if func is None:
            raise ValueError(
                'BLAS function %s could not be found' % func_name)
        func.module_name, func.typecode = module_name, prefix
        blas_funcs.append(func)

    if unpack:
        return blas_funcs[0]
    else:
        return blas_funcs
