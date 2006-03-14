"""
Wrappers to LAPACK library
==========================

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
postpone_import = 1
