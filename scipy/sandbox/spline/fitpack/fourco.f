      subroutine fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
c  subroutine fourco calculates the integrals
c                    /t(n-3)
c    ress(i) =      !        s(x)*sin(alfa(i)*x) dx    and
c              t(4)/
c                    /t(n-3)
c    resc(i) =      !        s(x)*cos(alfa(i)*x) dx, i=1,...,m,
c              t(4)/
c  where s(x) denotes a cubic spline which is given in its
c  b-spline representation.
c
c  calling sequence:
c     call fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
c
c  input parameters:
c    t    : real array,length n, containing the knots of s(x).
c    n    : integer, containing the total number of knots. n>=10.
c    c    : real array,length n, containing the b-spline coefficients.
c    alfa : real array,length m, containing the parameters alfa(i).
c    m    : integer, specifying the number of integrals to be computed.
c    wrk1 : real array,length n. used as working space
c    wrk2 : real array,length n. used as working space
c
c  output parameters:
c    ress : real array,length m, containing the integrals ress(i).
c    resc : real array,length m, containing the integrals resc(i).
c    ier  : error flag:
c      ier=0 : normal return.
c      ier=10: invalid input data (see restrictions).
c
c  restrictions:
c    n >= 10
c    t(4) < t(5) < ... < t(n-4) < t(n-3).
c    t(1) <= t(2) <= t(3) <= t(4).
c    t(n-3) <= t(n-2) <= t(n-1) <= t(n).
c
c  other subroutines required: fpbfou,fpcsin
c
c  references :
c    dierckx p. : calculation of fouriercoefficients of discrete
c                 functions using cubic splines. j. computational
c                 and applied mathematics 3 (1977) 207-209.
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
      integer n,m,ier
c  ..array arguments..
      real*8 t(n),c(n),wrk1(n),wrk2(n),alfa(m),ress(m),resc(m)
c  ..local scalars..
      integer i,j,n4
      real*8 rs,rc
c  ..
      n4 = n-4
c  before starting computations a data check is made. in the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(n.lt.10) go to 50
      j = n
      do 10 i=1,3
        if(t(i).gt.t(i+1)) go to 50
        if(t(j).lt.t(j-1)) go to 50
        j = j-1
  10  continue
      do 20 i=4,n4
        if(t(i).ge.t(i+1)) go to 50
  20  continue
      ier = 0
c  main loop for the different alfa(i).
      do 40 i=1,m
c  calculate the integrals
c    wrk1(j) = integral(nj,4(x)*sin(alfa*x))    and
c    wrk2(j) = integral(nj,4(x)*cos(alfa*x)),  j=1,2,...,n-4,
c  where nj,4(x) denotes the normalised cubic b-spline defined on the
c  knots t(j),t(j+1),...,t(j+4).
         call fpbfou(t,n,alfa(i),wrk1,wrk2)
c  calculate the integrals ress(i) and resc(i).
         rs = 0.
         rc = 0.
         do 30 j=1,n4
            rs = rs+c(j)*wrk1(j)
            rc = rc+c(j)*wrk2(j)
  30     continue
         ress(i) = rs
         resc(i) = rc
  40  continue
  50  return
      end
