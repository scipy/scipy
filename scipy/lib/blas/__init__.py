#
# BLAS wrappers
#

from warnings import warn
from info import __doc__

__all__ = ['fblas','cblas','get_blas_funcs']

import fblas
import cblas

from numpy import deprecate

_use_force_cblas = 1
if hasattr(cblas,'empty_module'):
    cblas = fblas
    _use_force_cblas = 0
elif hasattr(fblas,'empty_module'):
    fblas = cblas

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'} # 'd' will be default for 'i',..
_inv_type_conv = {'s':'f','d':'d','c':'F','z':'D'}

@deprecate
def get_blas_funcs(names,arrays=(),debug=0):
    """
    This function is deprecated, use scipy.linalg.get_blas_funcs instead.

    Return available BLAS function objects with names.
    arrays are used to determine the optimal prefix of
    BLAS routines.
    """

    ordering = []
    for i in range(len(arrays)):
        t = arrays[i].dtype.char
        if t not in _type_conv:
            t = 'd'
        ordering.append((t,i))
    if ordering:
        ordering.sort()
        required_prefix = _type_conv[ordering[0][0]]
    else:
        required_prefix = 'd'
    dtypechar = _inv_type_conv[required_prefix]
    # Default lookup:
    if ordering and arrays[ordering[0][1]].flags['FORTRAN']:
        # prefer Fortran code for leading array with column major order
        m1,m2 = fblas,cblas
    else:
        # in all other cases, C code is preferred
        m1,m2 = cblas,fblas
    funcs = []
    for name in names:
        if name=='ger' and dtypechar in 'FD':
            name = 'gerc'
        elif name in ('dotc', 'dotu') and dtypechar in 'fd':
            name = 'dot'
        func_name = required_prefix + name
        if name == 'nrm2' and dtypechar == 'D':
            func_name = 'dznrm2'
        elif name == 'nrm2' and dtypechar == 'F':
            func_name = 'scnrm2'
        func = getattr(m1,func_name,None)
        if func is None:
            func = getattr(m2,func_name)
            func.module_name = m2.__name__.split('.')[-1]
        else:
            func.module_name = m1.__name__.split('.')[-1]
        func.prefix = required_prefix
        func.dtypechar = dtypechar
        funcs.append(func)
    return tuple(funcs)

from numpy.testing import Tester
test = Tester().test
