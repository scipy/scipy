      subroutine spalde(t,n,c,k1,x,d,ier)
c  subroutine spalde evaluates at a point x all the derivatives
c              (j-1)
c      d(j) = s     (x) , j=1,2,...,k1
c  of a spline s(x) of order k1 (degree k=k1-1), given in its b-spline
c  representation.
c
c  calling sequence:
c     call spalde(t,n,c,k1,x,d,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k1   : integer, giving the order of s(x) (order=degree+1)
c    x    : real, which contains the point where the derivatives must
c           be evaluated.
c
c  output parameters:
c    d    : array,length k1, containing the derivative values of s(x).
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    t(k1) <= x <= t(n-k1+1)
c
c  further comments:
c    if x coincides with a knot, right derivatives are computed
c    ( left derivatives if x = t(n-k1+1) ).
c
c  other subroutines required: fpader.
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
c  ..scalar arguments..
      integer n,k1,ier
      real*8 x
c  ..array arguments..
      real*8 t(n),c(n),d(k1)
c  ..local scalars..
      integer l,nk1
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      nk1 = n-k1
      if(x.lt.t(k1) .or. x.gt.t(nk1+1)) go to 300
c  search for knot interval t(l) <= x < t(l+1)
      l = k1
 100  if(x.lt.t(l+1) .or. l.eq.nk1) go to 200
      l = l+1
      go to 100
 200  if(t(l).ge.t(l+1)) go to 300
      ier = 0
c  calculate the derivatives.
      call fpader(t,n,c,k1,x,l,d)
 300  return
      end
