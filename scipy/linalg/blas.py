#
# Author: Pearu Peterson, March 2002
#         refactoring by Fabian Pedregosa, March 2010
#

__all__ = ['get_blas_funcs']

import numpy
# The following ensures that possibly missing flavor (C or Fortran) is
# replaced with the available one. If none is available, exception
# is raised at the first attempt to use the resources.

from scipy.linalg import cblas, fblas
if hasattr(cblas,'empty_module'):
    cblas = fblas
elif hasattr(fblas,'empty_module'):
    fblas = cblas

 # 'd' will be default for 'i',..
_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}


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
        used.

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

    n_arrays = len(arrays)
    blas_funcs = []
    unpack = False
    dtype = numpy.dtype(dtype)

    if isinstance(names, str):
        names = (names,)
        unpack = True

    for i, name in enumerate(names):
        module1 = (cblas, 'cblas')
        module2 = (fblas, 'fblas')

        if i < n_arrays:
            prefix = _type_conv.get(arrays[i].dtype.char, 'd')
            if arrays[i].flags['FORTRAN']:
                module1, module2 = module2, module1
        else:
            prefix = _type_conv.get(dtype.char, 'd')

        func_name = prefix + name
        func = getattr(module1[0], func_name, None)
        module_name = module1[1]
        if func is None:
            func = getattr(module2[0], func_name, None)
            module_name = module2[1]
        if func is None:
            raise ValueError(
                'BLAS function %s could not be found' % func_name)
        func.module_name = module_name
        func.typecode = prefix
        blas_funcs.append(func)

    if unpack:
        return blas_funcs[0]
    else:
        return blas_funcs
