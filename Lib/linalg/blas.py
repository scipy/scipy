#
# Author: Pearu Peterson, March 2002
#

__all__ = ['get_blas_funcs']

import string

# The following ensures that possibly missing flavor (C or Fortran) is
# replaced with the available one. If none is available, exception
# is raised at the first attempt to use the resources.

import fblas
import cblas
if hasattr(cblas,'empty_module'):
    cblas = fblas
elif hasattr(fblas,'empty_module'):
    fblas = cblas

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'} # 'd' will be default for 'i',..
_inv_type_conv = {'s':'f','d':'d','c':'F','z':'D'}

def get_blas_funcs(names,arrays=(),debug=0):
    """Return available BLAS function objects with names.
    arrays are used to determine the optimal prefix of
    BLAS routines."""
    ordering = []
    for i in range(len(arrays)):
        t = arrays[i].typecode()
        if not _type_conv.has_key(t): t = 'd'
        ordering.append((t,i))
    if ordering:
        ordering.sort()
        required_prefix = _type_conv[ordering[0][0]]
    else:
        required_prefix = 'd'
    typecode = _inv_type_conv[required_prefix]
    # Default lookup:
    if ordering and fblas.has_column_major_storage(arrays[ordering[0][1]]):
        # prefer Fortran code for leading array with column major order
        m1,m2 = fblas,cblas
    else:
        # in all other cases, C code is preferred
        m1,m2 = cblas,fblas
    funcs = []
    for name in names:
        if name=='ger' and typecode in 'FD':
            name = 'gerc'
        func_name = required_prefix + name
        func = getattr(m1,func_name,None)
        if func is None:
            func = getattr(m2,func_name)
            func.module_name = string.split(m2.__name__,'.')[-1]
        else:
            func.module_name = string.split(m1.__name__,'.')[-1]
        func.prefix = required_prefix
        func.typecode = typecode
        funcs.append(func)
    return tuple(funcs)
