      subroutine insert(iopt,t,n,c,k,x,tt,nn,cc,nest,ier)
c  subroutine insert inserts a new knot x into a spline function s(x)
c  of degree k and calculates the b-spline representation of s(x) with
c  respect to the new set of knots. in addition, if iopt.ne.0, s(x)
c  will be considered as a periodic spline with period per=t(n-k)-t(k+1)
c  satisfying the boundary constraints
c       t(i+n-2*k-1) = t(i)+per  ,i=1,2,...,2*k+1
c       c(i+n-2*k-1) = c(i)      ,i=1,2,...,k
c  in that case, the knots and b-spline coefficients returned will also
c  satisfy these boundary constraints, i.e.
c       tt(i+nn-2*k-1) = tt(i)+per  ,i=1,2,...,2*k+1
c       cc(i+nn-2*k-1) = cc(i)      ,i=1,2,...,k
c
c  calling sequence:
c     call insert(iopt,t,n,c,k,x,tt,nn,cc,nest,ier)
c
c  input parameters:
c    iopt : integer flag, specifying whether (iopt.ne.0) or not (iopt=0)
c           the given spline must be considered as being periodic.
c    t    : array,length nest, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length nest, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    x    : real, which gives the location of the knot to be inserted.
c    nest : integer specifying the dimension of the arrays t,c,tt and cc
c           nest > n.
c
c  output parameters:
c    tt   : array,length nest, which contains the position of the knots
c           after insertion.
c    nn   : integer, giving the total number of knots after insertion
c    cc   : array,length nest, which contains the b-spline coefficients
c           of s(x) with respect to the new set of knots.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    nest > n
c    t(k+1) <= x <= t(n-k)
c    in case of a periodic spline (iopt.ne.0) there must be
c       either at least k interior knots t(j) satisfying t(k+1)<t(j)<=x
c       or at least k interior knots t(j) satisfying x<=t(j)<t(n-k)
c
c  other subroutines required: fpinst.
c
c  further comments:
c   subroutine insert may be called as follows
c        call insert(iopt,t,n,c,k,x,t,n,c,nest,ier)
c   in which case the new representation will simply replace the old one
c
c  references :
c    boehm w : inserting new knots into b-spline curves. computer aided
c              design 12 (1980) 199-201.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : february 2007 (second interval search added)
c
c  ..scalar arguments..
      integer iopt,n,k,nn,nest,ier
      real*8 x
c  ..array arguments..
      real*8 t(nest),c(nest),tt(nest),cc(nest)
c  ..local scalars..
      integer kk,k1,l,nk
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nest.le.n) go to 40
      k1 = k+1
      nk = n-k
      if(x.lt.t(k1) .or. x.gt.t(nk)) go to 40
c  search for knot interval t(l) <= x < t(l+1).
      l = k1
  10  if(x.lt.t(l+1)) go to 20
      l = l+1
      if(l.eq.nk) go to 14
      go to 10
c  if no interval found above, then reverse the search and 
c  look for knot interval t(l) < x <= t(l+1).
  14  l = nk-1
  16  if(x.gt.t(l)) go to 20
      l = l-1
      if(l.eq.k) go to 40
      go to 16
  20  if(t(l).ge.t(l+1)) go to 40
      if(iopt.eq.0) go to 30
      kk = 2*k
      if(l.le.kk .and. l.ge.(n-kk)) go to 40
  30  ier = 0
c  insert the new knot.
      call fpinst(iopt,t,n,c,k,x,l,tt,nn,cc,nest)
  40  return
      end
