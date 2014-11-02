      subroutine splder(t,n,c,k,nu,x,y,m,e,wrk,ier)
c  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
c  the derivative of order nu of a spline s(x) of degree k,given in
c  its b-spline representation.
c
c  calling sequence:
c     call splder(t,n,c,k,nu,x,y,m,e,wrk,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    nu   : integer, specifying the order of the derivative. 0<=nu<=k
c    x    : array,length m, which contains the points where the deriv-
c           ative of s(x) must be evaluated.
c    m    : integer, giving the number of points where the derivative
c           of s(x) must be evaluated
c    e    : integer, if 0 the spline is extrapolated from the end
c           spans for points not in the support, if 1 the spline
c           evaluates to zero for those points, and if 2 ier is set to
c           1 and the subroutine returns.
c    wrk  : real array of dimension n. used as working space.
c
c  output parameters:
c    y    : array,length m, giving the value of the derivative of s(x)
c           at the different points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier = 1 : argument out of bounds and e == 2
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    0 <= nu <= k
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c++ pearu: 13 aug 20003
c++   - disabled cliping x values to interval [min(t),max(t)]
c++   - removed the restriction of the orderness of x values
c++   - fixed initialization of sp to double precision value
c
c  ..scalar arguments..
      integer n,k,nu,m,e,ier
c  ..array arguments..
      real*8 t(n),c(n),x(m),y(m),wrk(n)
c  ..local scalars..
      integer i,j,kk,k1,k2,l,ll,l1,l2,nk1,nk2,nn
      real*8 ak,arg,fac,sp,tb,te
c++..
      integer k3
c..++
c  ..local arrays ..
      real*8 h(6)
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nu.lt.0 .or. nu.gt.k) go to 200
c--      if(m-1) 200,30,10
c++..
      if(m.lt.1) go to 200
c..++
c--  10  do 20 i=2,m
c--        if(x(i).lt.x(i-1)) go to 200
c--  20  continue
      ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      k3 = k1+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
c  the derivative of order nu of a spline of degree k is a spline of
c  degree k-nu,the b-spline coefficients wrk(i) of which can be found
c  using the recurrence scheme of de boor.
      l = 1
      kk = k
      nn = n
      do 40 i=1,nk1
         wrk(i) = c(i)
  40  continue
      if(nu.eq.0) go to 100
      nk2 = nk1
      do 60 j=1,nu
         ak = kk
         nk2 = nk2-1
         l1 = l
         do 50 i=1,nk2
            l1 = l1+1
            l2 = l1+kk
            fac = t(l2)-t(l1)
            if(fac.le.0.) go to 50
            wrk(i) = ak*(wrk(i+1)-wrk(i))/fac
  50     continue
         l = l+1
         kk = kk-1
  60  continue
      if(kk.ne.0) go to 100
c  if nu=k the derivative is a piecewise constant function
      j = 1
      do 90 i=1,m
        arg = x(i)
c++..
c  check if arg is in the support
        if (arg .lt. tb .or. arg .gt. te) then
            if (e .eq. 0) then
                goto 65
            else if (e .eq. 1) then
                y(i) = 0
                goto 90
            else if (e .eq. 2) then
                ier = 1
                goto 200
            endif
        endif
c  search for knot interval t(l) <= arg < t(l+1)
 65     if(arg.ge.t(l) .or. l+1.eq.k2) go to 70
        l1 = l
        l = l-1
        j = j-1
        go to 65
c..++
  70    if(arg.lt.t(l+1) .or. l.eq.nk1) go to 80
        l = l+1
        j = j+1
        go to 70
  80    y(i) = wrk(j)
  90  continue
      go to 200

 100  l = k1
      l1 = l+1
      k2 = k1-nu
c  main loop for the different points.
      do 180 i=1,m
c  fetch a new x-value arg.
        arg = x(i)
c  check if arg is in the support
        if (arg .lt. tb .or. arg .gt. te) then
            if (e .eq. 0) then
                goto 135
            else if (e .eq. 1) then
                y(i) = 0
                goto 180
            else if (e .eq. 2) then
                ier = 1
                goto 200
            endif
        endif
c  search for knot interval t(l) <= arg < t(l+1)
 135    if(arg.ge.t(l) .or. l1.eq.k3) go to 140
        l1 = l
        l = l-1
        go to 135
c..++
 140    if(arg.lt.t(l1) .or. l.eq.nk1) go to 150
        l = l1
        l1 = l+1
        go to 140
c  evaluate the non-zero b-splines of degree k-nu at arg.
 150    call fpbspl(t,n,kk,arg,l,h)
c  find the value of the derivative at x=arg.
        sp = 0.0d0
        ll = l-k1
        do 160 j=1,k2
          ll = ll+1
          sp = sp+wrk(ll)*h(j)
 160    continue
        y(i) = sp
 180  continue
 200  return
      end
