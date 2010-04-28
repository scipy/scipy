#
# Author: Pearu Peterson, March 2002
#

__all__ = ['get_lapack_funcs']

# The following ensures that possibly missing flavor (C or Fortran) is
# replaced with the available one. If none is available, exception
# is raised at the first attempt to use the resources.
import types

import numpy

from scipy.linalg import flapack
from scipy.linalg import clapack
_use_force_clapack = 1
if hasattr(clapack,'empty_module'):
    clapack = flapack
    _use_force_clapack = 0
elif hasattr(flapack,'empty_module'):
    flapack = clapack

def cast_to_lapack_prefix(t):
    if issubclass(t, numpy.single):
        prefix = 's'
    elif issubclass(t, numpy.double):
        prefix = 'd'
    elif issubclass(t, numpy.longdouble):
        prefix = 'd'
    elif issubclass(t, numpy.csingle):
        prefix = 'c'
    elif issubclass(t, numpy.cdouble):
        prefix = 'z'
    elif issubclass(t, numpy.clongdouble):
        prefix = 'z'
    else:
        prefix = 'd'
    return prefix

prefix_to_order = dict(s=3, d=2, c=1, z=0)
order_to_prefix = ['s', 'd', 'c', 'z']
prefix_to_dtype = dict(s=numpy.single, d=numpy.double,
                       c=numpy.csingle, z=numpy.cdouble)

def find_best_lapack_type(arrays):
    if not arrays:
        return 'd', numpy.double, False
    ordering = []
    for i in range(len(arrays)):
        t = arrays[i].dtype.type
        prefix = cast_to_lapack_prefix(t)
        order = prefix_to_order[prefix]
        ordering.append((order, prefix, i))
    ordering.sort()
    _, required_prefix, lowest_array_index = ordering[0]
    dtype = prefix_to_dtype[required_prefix]
    isfortran = numpy.isfortran(arrays[lowest_array_index])
    return required_prefix, dtype, isfortran

def get_lapack_funcs(names, arrays=()):
    """Return available LAPACK function objects with names.
    arrays are used to determine the optimal prefix of
    LAPACK routines.
    """
    #If force_clapack is True then available Atlas routine
    #is returned for column major storaged arrays with
    #rowmajor argument set to False.
    force_clapack=False  #XXX: Don't set it true! The feature is unreliable
                         #     and may cause incorrect results.
                         #     See test_basic.test_solve.check_20Feb04_bug.

    required_prefix, dtype, isfortran = find_best_lapack_type(arrays)
    # Default lookup:
    if isfortran:
        # prefer Fortran code for leading array with column major order
        m1, m2 = flapack, clapack
    else:
        # in all other cases, C code is preferred
        m1, m2 = clapack, flapack
    if not _use_force_clapack:
        force_clapack = False
    funcs = []
    m1_name = m1.__name__.split('.')[-1]
    m2_name = m2.__name__.split('.')[-1]
    for name in names:
        func_name = required_prefix + name
        func = getattr(m1,func_name,None)
        if func is None:
            func = getattr(m2,func_name)
            func.module_name = m2_name
        else:
            func.module_name = m1_name
            if force_clapack and m1 is flapack:
                func2 = getattr(m2,func_name,None)
                if func2 is not None:
                    exec _colmajor_func_template % {'func_name':func_name}
                    func = types.FunctionType(func_code,
                                              {'clapack_func':func2},
                                              func_name)
                    func.module_name = m2_name
                    func.__doc__ = func2.__doc__
        func.prefix = required_prefix
        func.dtype = dtype
        funcs.append(func)
    return tuple(funcs)



_colmajor_func_template = '''\
def %(func_name)s(*args,**kws):
    if "rowmajor" not in kws:
        kws["rowmajor"] = 0
    return clapack_func(*args,**kws)
func_code = %(func_name)s.func_code
'''
