      real*8 function splint(t,n,c,k,a,b,wrk)
      implicit none
c  function splint calculates the integral of a spline function s(x)
c  of degree k, which is given in its normalized b-spline representation
c
c  calling sequence:
c     aint = splint(t,n,c,k,a,b,wrk)
c
c  input parameters:
c    t    : array,length n,which contains the position of the knots
c           of s(x).
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, containing the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    a,b  : real values, containing the end points of the integration
c           interval. s(x) is considered to be identically zero outside
c           the interval (t(k+1),t(n-k)).
c
c  output parameter:
c    aint : real, containing the integral of s(x) between a and b.
c    wrk  : real array, length n.  used as working space
c           on output, wrk will contain the integrals of the normalized
c           b-splines defined on the set of knots.
c
c  other subroutines required: fpintb.
c
c  references :
c    gaffney p.w. : the calculation of indefinite integrals of b-splines
c                   j. inst. maths applics 17 (1976) 37-41.
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
      real*8 a,b
      integer n,k
c  ..array arguments..
      real*8 t(n),c(n),wrk(n)
c  ..local scalars..
      integer i,nk1
c  ..
      nk1 = n-k-1
c  calculate the integrals wrk(i) of the normalized b-splines
c  ni,k+1(x), i=1,2,...nk1.
      call fpintb(t,n,wrk,nk1,a,b)
c  calculate the integral of s(x).
      splint = 0.0d0
      do 10 i=1,nk1
        splint = splint+c(i)*wrk(i)
  10  continue
      return
      end
