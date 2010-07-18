      subroutine splev(t,n,c,k,x,y,m,e,ier)
c  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
c  a spline s(x) of degree k, given in its b-spline representation.
c
c  calling sequence:
c     call splev(t,n,c,k,x,y,m,e,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    x    : array,length m, which contains the points where s(x) must
c           be evaluated.
c    m    : integer, giving the number of points where s(x) must be
c           evaluated.
c    e    : integer, if 0 the spline is extrapolated from the end
c           spans for points not in the support, if 1 the spline
c           evaluates to zero for those points, and if 2 ier is set to
c           1 and the subroutine returns.
c
c  output parameter:
c    y    : array,length m, giving the value of s(x) at the different
c           points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier = 1 : argument out of bounds and e == 2
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    m >= 1
c--    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl.
c
c  references :
c    de boor c  : on calculating with b-splines, j. approximation theory
c                 6 (1972) 50-62.
c    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
c                 applics 10 (1972) 134-149.
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
c++ pearu: 11 aug 2003
c++   - disabled cliping x values to interval [min(t),max(t)]
c++   - removed the restriction of the orderness of x values
c++   - fixed initialization of sp to double precision value
c
c  ..scalar arguments..
      integer n, k, m, e, ier
c  ..array arguments..
      real*8 t(n), c(n), x(m), y(m)
c  ..local scalars..
      integer i, j, k1, l, ll, l1, nk1
c++..
      integer k2
c..++
      real*8 arg, sp, tb, te
c  ..local array..
      real*8 h(20)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
c--      if(m-1) 100,30,10
c++..
      if (m .lt. 1) go to 100
c..++
c--  10  do 20 i=2,m
c--        if(x(i).lt.x(i-1)) go to 100
c--  20  continue
  30  ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k + 1
c++..
      k2 = k1 + 1
c..++
      nk1 = n - k1
      tb = t(k1)
      te = t(nk1 + 1)
      l = k1
      l1 = l + 1
c  main loop for the different points.
      do 80 i = 1, m
c  fetch a new x-value arg.
        arg = x(i)
c  check if arg is in the support
        if (arg .lt. tb .or. arg .gt. te) then
            if (e .eq. 0) then
                goto 35
            else if (e .eq. 1) then
                y(i) = 0
                goto 80
            else if (e .eq. 2) then
                ier = 1
                goto 100
            endif
        endif
c  search for knot interval t(l) <= arg < t(l+1)
c++..
 35     if (arg .ge. t(l) .or. l1 .eq. k2) go to 40
        l1 = l
        l = l - 1
        go to 35
c..++
  40    if(arg .lt. t(l1) .or. l .eq. nk1) go to 50
        l = l1
        l1 = l + 1
        go to 40
c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t, n, k, arg, l, h)
c  find the value of s(x) at x=arg.
        sp = 0.0d0
        ll = l - k1
        do 60 j = 1, k1
          ll = ll + 1
          sp = sp + c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end
