"""
Wrappers to LAPACK library
==========================

NOTE: this module is deprecated -- use scipy.linalg.lapack instead!

  flapack -- wrappers for Fortran [*] LAPACK routines
  clapack -- wrappers for ATLAS LAPACK routines
  calc_lwork -- calculate optimal lwork parameters
  get_lapack_funcs -- query for wrapper functions.

[*] If ATLAS libraries are available then Fortran routines
    actually use ATLAS routines and should perform equally
    well to ATLAS routines.

Module flapack
++++++++++++++

In the following all function names are shown without
type prefix (s,d,c,z). Optimal values for lwork can
be computed using calc_lwork module.

Linear Equations
----------------

  Drivers::

    lu,piv,x,info = gesv(a,b,overwrite_a=0,overwrite_b=0)
    lub,piv,x,info = gbsv(kl,ku,ab,b,overwrite_ab=0,overwrite_b=0)
    c,x,info = posv(a,b,lower=0,overwrite_a=0,overwrite_b=0)

  Computational routines::

    lu,piv,info = getrf(a,overwrite_a=0)
    x,info = getrs(lu,piv,b,trans=0,overwrite_b=0)
    inv_a,info = getri(lu,piv,lwork=min_lwork,overwrite_lu=0)

    c,info = potrf(a,lower=0,clean=1,overwrite_a=0)
    x,info = potrs(c,b,lower=0,overwrite_b=0)
    inv_a,info = potri(c,lower=0,overwrite_c=0)

    inv_c,info = trtri(c,lower=0,unitdiag=0,overwrite_c=0)

Linear Least Squares (LLS) Problems
-----------------------------------

  Drivers::

    v,x,s,rank,info = gelss(a,b,cond=-1.0,lwork=min_lwork,overwrite_a=0,overwrite_b=0)

  Computational routines::

    qr,tau,info = geqrf(a,lwork=min_lwork,overwrite_a=0)
    q,info = orgqr|ungqr(qr,tau,lwork=min_lwork,overwrite_qr=0,overwrite_tau=1)

Generalized Linear Least Squares (LSE and GLM) Problems
-------------------------------------------------------

Standard Eigenvalue and Singular Value Problems
-----------------------------------------------

  Drivers::

    w,v,info = syev|heev(a,compute_v=1,lower=0,lwork=min_lwork,overwrite_a=0)
    w,v,info = syevd|heevd(a,compute_v=1,lower=0,lwork=min_lwork,overwrite_a=0)
    w,v,info = syevr|heevr(a,compute_v=1,lower=0,vrange=,irange=,atol=-1.0,lwork=min_lwork,overwrite_a=0)
    t,sdim,(wr,wi|w),vs,info = gees(select,a,compute_v=1,sort_t=0,lwork=min_lwork,select_extra_args=(),overwrite_a=0)
    wr,(wi,vl|w),vr,info = geev(a,compute_vl=1,compute_vr=1,lwork=min_lwork,overwrite_a=0)
    u,s,vt,info = gesdd(a,compute_uv=1,lwork=min_lwork,overwrite_a=0)

  Computational routines::

    ht,tau,info = gehrd(a,lo=0,hi=n-1,lwork=min_lwork,overwrite_a=0)
    ba,lo,hi,pivscale,info = gebal(a,scale=0,permute=0,overwrite_a=0)

Generalized Eigenvalue and Singular Value Problems
--------------------------------------------------

  Drivers::

    w,v,info = sygv|hegv(a,b,itype=1,compute_v=1,lower=0,lwork=min_lwork,overwrite_a=0,overwrite_b=0)
    w,v,info = sygvd|hegvd(a,b,itype=1,compute_v=1,lower=0,lwork=min_lwork,overwrite_a=0,overwrite_b=0)
    (alphar,alphai|alpha),beta,vl,vr,info = ggev(a,b,compute_vl=1,compute_vr=1,lwork=min_lwork,overwrite_a=0,overwrite_b=0)


Auxiliary routines
------------------

  a,info = lauum(c,lower=0,overwrite_c=0)
  a = laswp(a,piv,k1=0,k2=len(piv)-1,off=0,inc=1,overwrite_a=0)

Module clapack
++++++++++++++

Linear Equations
----------------

  Drivers::

    lu,piv,x,info = gesv(a,b,rowmajor=1,overwrite_a=0,overwrite_b=0)
    c,x,info = posv(a,b,lower=0,rowmajor=1,overwrite_a=0,overwrite_b=0)

  Computational routines::

    lu,piv,info = getrf(a,rowmajor=1,overwrite_a=0)
    x,info = getrs(lu,piv,b,trans=0,rowmajor=1,overwrite_b=0)
    inv_a,info = getri(lu,piv,rowmajor=1,overwrite_lu=0)

    c,info = potrf(a,lower=0,clean=1,rowmajor=1,overwrite_a=0)
    x,info = potrs(c,b,lower=0,rowmajor=1,overwrite_b=0)
    inv_a,info = potri(c,lower=0,rowmajor=1,overwrite_c=0)

    inv_c,info = trtri(c,lower=0,unitdiag=0,rowmajor=1,overwrite_c=0)

Auxiliary routines
------------------

  a,info = lauum(c,lower=0,rowmajor=1,overwrite_c=0)

Module calc_lwork
+++++++++++++++++

Optimal lwork is maxwrk. Default is minwrk.

  minwrk,maxwrk = gehrd(prefix,n,lo=0,hi=n-1)
  minwrk,maxwrk = gesdd(prefix,m,n,compute_uv=1)
  minwrk,maxwrk = gelss(prefix,m,n,nrhs)
  minwrk,maxwrk = getri(prefix,n)
  minwrk,maxwrk = geev(prefix,n,compute_vl=1,compute_vr=1)
  minwrk,maxwrk = heev(prefix,n,lower=0)
  minwrk,maxwrk = syev(prefix,n,lower=0)
  minwrk,maxwrk = gees(prefix,n,compute_v=1)
  minwrk,maxwrk = geqrf(prefix,m,n)
  minwrk,maxwrk = gqr(prefix,m,n)


"""

from __future__ import division, print_function, absolute_import

__all__ = ['get_lapack_funcs','calc_lwork','flapack','clapack']

from numpy import deprecate

from . import calc_lwork

# The following ensures that possibly missing flavor (C or Fortran) is
# replaced with the available one. If none is available, exception
# is raised at the first attempt to use the resources.


@deprecate(old_name="scipy.lib.lapack", new_name="scipy.linalg.lapack")
def _deprecated():
    pass
try:
    _deprecated()
except DeprecationWarning as e:
    # don't fail import if DeprecationWarnings raise error -- works around
    # the situation with Numpy's test framework
    pass

from . import flapack
from . import clapack

_use_force_clapack = 1
if hasattr(clapack,'empty_module'):
    clapack = flapack
    _use_force_clapack = 0
elif hasattr(flapack,'empty_module'):
    flapack = clapack


_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}  # 'd' will be default for 'i',..
_inv_type_conv = {'s':'f','d':'d','c':'F','z':'D'}


@deprecate
def get_lapack_funcs(names,arrays=(),debug=0,force_clapack=1):
    """Return available LAPACK function objects with names.
    arrays are used to determine the optimal prefix of
    LAPACK routines.
    If force_clapack is True then available Atlas routine
    is returned for column major storaged arrays with
    rowmajor argument set to False.
    """
    force_clapack = 0  # XXX: Don't set it true! The feature is unreliable
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
                    exec(_colmajor_func_template % {'func_name':func_name})
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
