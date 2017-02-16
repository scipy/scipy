      subroutine curev(idim,t,n,c,nc,k,u,m,x,mx,ier)
c  subroutine curev evaluates in a number of points u(i),i=1,2,...,m
c  a spline curve s(u) of degree k and dimension idim, given in its
c  b-spline representation.
c
c  calling sequence:
c     call curev(idim,t,n,c,nc,k,u,m,x,mx,ier)
c
c  input parameters:
c    idim : integer, giving the dimension of the spline curve.
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(u).
c    c    : array,length nc, which contains the b-spline coefficients.
c    nc   : integer, giving the total number of coefficients of s(u).
c    k    : integer, giving the degree of s(u).
c    u    : array,length m, which contains the points where s(u) must
c           be evaluated.
c    m    : integer, giving the number of points where s(u) must be
c           evaluated.
c    mx   : integer, giving the dimension of the array x. mx >= m*idim
c
c  output parameters:
c    x    : array,length mx,giving the value of s(u) at the different
c           points. x(idim*(i-1)+j) will contain the j-th coordinate
c           of the i-th point on the curve.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    m >= 1
c    mx >= m*idim
c    t(k+1) <= u(i) <= u(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl.
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer idim,n,nc,k,m,mx,ier
c  ..array arguments..
      real*8 t(n),c(nc),u(m),x(mx)
c  ..local scalars..
      integer i,j,jj,j1,k1,l,ll,l1,mm,nk1
      real*8 arg,sp,tb,te
c  ..local array..
      real*8 h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if (m.lt.1) go to 100
      if (m.eq.1) go to 30
      go to 10
  10  do 20 i=2,m
        if(u(i).lt.u(i-1)) go to 100
  20  continue
  30  if(mx.lt.(m*idim)) go to 100
      ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
c  main loop for the different points.
      mm = 0
      do 80 i=1,m
c  fetch a new u-value arg.
        arg = u(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
c  find the value of s(u) at u=arg.
        ll = l-k1
        do 70 j1=1,idim
          jj = ll
          sp = 0.
          do 60 j=1,k1
            jj = jj+1
            sp = sp+c(jj)*h(j)
  60      continue
          mm = mm+1
          x(mm) = sp
          ll = ll+n
  70    continue
  80  continue
 100  return
      end
