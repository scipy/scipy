#
# LAPACK wrappers
#

from info import __doc__

__all__ = ['get_lapack_funcs','calc_lwork','flapack','clapack']

import calc_lwork

# The following ensures that possibly missing flavor (C or Fortran) is
# replaced with the available one. If none is available, exception
# is raised at the first attempt to use the resources.

import flapack
import clapack

_use_force_clapack = 1
if hasattr(clapack,'empty_module'):
    clapack = flapack
    _use_force_clapack = 0
elif hasattr(flapack,'empty_module'):
    flapack = clapack

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'} # 'd' will be default for 'i',..
_inv_type_conv = {'s':'f','d':'d','c':'F','z':'D'}

def get_lapack_funcs(names,arrays=(),debug=0,force_clapack=1):
    """Return available LAPACK function objects with names.
    arrays are used to determine the optimal prefix of
    LAPACK routines.
    If force_clapack is True then available Atlas routine
    is returned for column major storaged arrays with
    rowmajor argument set to False.
    """
    force_clapack=0  #XXX: Don't set it true! The feature is unreliable
                     #     and may cause incorrect results.
                     #     See test_basic.test_solve.check_20Feb04_bug.

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
        m1,m2 = flapack,clapack
    else:
        # in all other cases, C code is preferred
        m1,m2 = clapack,flapack
    if not _use_force_clapack:
        force_clapack = 0
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
                    import new
                    exec _colmajor_func_template % {'func_name':func_name}
                    func = new.function(func_code,{'clapack_func':func2},func_name)
                    func.module_name = m2_name
                    func.__doc__ = func2.__doc__
        func.prefix = required_prefix
        func.dtypechar = dtypechar
        funcs.append(func)
    return tuple(funcs)

_colmajor_func_template = '''\
def %(func_name)s(*args,**kws):
    if "rowmajor" not in kws:
        kws["rowmajor"] = 0
    return clapack_func(*args,**kws)
func_code = %(func_name)s.func_code
'''

from numpy.testing import Tester
test = Tester().test
