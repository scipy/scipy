"""
Wrappers to BLAS library
========================

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
