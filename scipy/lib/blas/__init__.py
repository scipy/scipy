"""
Wrappers to BLAS library
========================

NOTE: this module is deprecated -- use scipy.linalg.blas instead!

fblas -- wrappers for Fortran [*] BLAS routines
cblas -- wrappers for ATLAS BLAS routines
get_blas_funcs -- query for wrapper functions.

[*] If ATLAS libraries are available then Fortran routines
    actually use ATLAS routines and should perform equally
    well to ATLAS routines.

Module fblas
++++++++++++

In the following all function names are shown without type prefixes.

Level 1 routines
----------------

  c,s = rotg(a,b)
  param = rotmg(d1,d2,x1,y1)
  x,y = rot(x,y,c,s,n=(len(x)-offx)/abs(incx),offx=0,incx=1,offy=0,incy=1,overwrite_x=0,overwrite_y=0)
  x,y = rotm(x,y,param,n=(len(x)-offx)/abs(incx),offx=0,incx=1,offy=0,incy=1,overwrite_x=0,overwrite_y=0)
  x,y = swap(x,y,n=(len(x)-offx)/abs(incx),offx=0,incx=1,offy=0,incy=1)
  x = scal(a,x,n=(len(x)-offx)/abs(incx),offx=0,incx=1)
  y = copy(x,y,n=(len(x)-offx)/abs(incx),offx=0,incx=1,offy=0,incy=1)
  y = axpy(x,y,n=(len(x)-offx)/abs(incx),a=1.0,offx=0,incx=1,offy=0,incy=1)
  xy = dot(x,y,n=(len(x)-offx)/abs(incx),offx=0,incx=1,offy=0,incy=1)
  xy = dotu(x,y,n=(len(x)-offx)/abs(incx),offx=0,incx=1,offy=0,incy=1)
  xy = dotc(x,y,n=(len(x)-offx)/abs(incx),offx=0,incx=1,offy=0,incy=1)
  n2 = nrm2(x,n=(len(x)-offx)/abs(incx),offx=0,incx=1)
  s = asum(x,n=(len(x)-offx)/abs(incx),offx=0,incx=1)
  k = amax(x,n=(len(x)-offx)/abs(incx),offx=0,incx=1)

  Prefixes:
    rotg,swap,copy,axpy: s,d,c,z
    amax: is,id,ic,iz
    asum,nrm2: s,d,sc,dz
    scal: s,d,c,z,sc,dz
    rotm,rotmg,dot: s,d
    dotu,dotc: c,z
    rot: s,d,cs,zd

Level 2 routines
----------------

  y = gemv(alpha,a,x,beta=0.0,y=,offx=0,incx=1,offy=0,incy=1,trans=0,overwrite_y=0)
  y = symv(alpha,a,x,beta=0.0,y=,offx=0,incx=1,offy=0,incy=1,lower=0,overwrite_y=0)
  y = hemv(alpha,a,x,beta=(0.0, 0.0),y=,offx=0,incx=1,offy=0,incy=1,lower=0,overwrite_y=0)
  x = trmv(a,x,offx=0,incx=1,lower=0,trans=0,unitdiag=0,overwrite_x=0)
  a = ger(alpha,x,y,incx=1,incy=1,a=0.0,overwrite_x=1,overwrite_y=1,overwrite_a=0)
  a = ger{u|c}(alpha,x,y,incx=1,incy=1,a=(0.0,0.0),overwrite_x=1,overwrite_y=1,overwrite_a=0)

  Prefixes:
   gemv, trmv: s,d,c,z
   symv,ger: s,d
   hemv,geru,gerc: c,z

Level 3 routines
----------------

  c = gemm(alpha,a,b,beta=0.0,c=,trans_a=0,trans_b=0,overwrite_c=0)

  Prefixes:
    gemm: s,d,c,z

Module cblas
++++++++++++

In the following all function names are shown without type prefixes.

Level 1 routines
----------------

  z = axpy(x,y,n=len(x)/abs(incx),a=1.0,incx=1,incy=incx,overwrite_y=0)

  Prefixes:
    axpy: s,d,c,z
"""
from __future__ import division, print_function, absolute_import

from warnings import warn

__all__ = ['fblas','cblas','get_blas_funcs']

from . import fblas
from . import cblas

from numpy import deprecate


@deprecate(old_name="scipy.lib.blas", new_name="scipy.linalg.blas")
def _deprecated():
    pass
try:
    _deprecated()
except DeprecationWarning as e:
    # don't fail import if DeprecationWarnings raise error -- works around
    # the situation with Numpy's test framework
    pass

_use_force_cblas = 1
if hasattr(cblas,'empty_module'):
    cblas = fblas
    _use_force_cblas = 0
elif hasattr(fblas,'empty_module'):
    fblas = cblas


_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}  # 'd' will be default for 'i',..
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
        if name == 'ger' and dtypechar in 'FD':
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
