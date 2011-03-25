#
# Author: Pearu Peterson, March 2002
#         refactoring by Fabian Pedregosa, March 2010
#

__all__ = ['get_blas_funcs']


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
_inv_type_conv = {'s':'f','d':'d','c':'F','z':'D'}


def get_blas_funcs(names, arrays=()):
    """Return available BLAS function objects from names.

    Arrays are used to determine the optimal prefix of BLAS routines.

    Parameters
    ----------
    names : iterable of strings
        Names of BLAS functions withouth type prefix.

    arrays : iterable of arrays, optional
        Arrays can be given to determine optiomal prefix of BLAS
        routines. If not given, double-precision routines will be
        used.

    Returns
    -------
    funcs : list
        List of functions.

    Notes
    -----
    This routines automatically chooses between Fortran/C
    interfaces. Fortran code is used whenever possible for arrays with
    column major order. In all other cases, C code is preferred.
    """

    n_arrays = len(arrays)
    blas_funcs = []

    for i, name in enumerate(names):
        module1 = (cblas, 'cblas')
        module2 = (fblas, 'fblas')
        prefix = 'd'

        if i < n_arrays:
            prefix = _type_conv.get(arrays[i].dtype.char, 'd')
            if arrays[i].flags['FORTRAN']:
                module1, module2 = module2, module1

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
        func.typecode = _inv_type_conv[prefix]
        blas_funcs.append(func)
    return blas_funcs
