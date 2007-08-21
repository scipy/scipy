      subroutine cualde(idim,t,n,c,nc,k1,u,d,nd,ier)
c  subroutine cualde evaluates at the point u all the derivatives
c                     (l)
c     d(idim*l+j) = sj   (u) ,l=0,1,...,k, j=1,2,...,idim
c  of a spline curve s(u) of order k1 (degree k=k1-1) and dimension idim
c  given in its b-spline representation.
c
c  calling sequence:
c     call cualde(idim,t,n,c,nc,k1,u,d,nd,ier)
c
c  input parameters:
c    idim : integer, giving the dimension of the spline curve.
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(u).
c    c    : array,length nc, which contains the b-spline coefficients.
c    nc   : integer, giving the total number of coefficients of s(u).
c    k1   : integer, giving the order of s(u) (order=degree+1).
c    u    : real, which contains the point where the derivatives must
c           be evaluated.
c    nd   : integer, giving the dimension of the array d. nd >= k1*idim
c
c  output parameters:
c    d    : array,length nd,giving the different curve derivatives.
c           d(idim*l+j) will contain the j-th coordinate of the l-th
c           derivative of the curve at the point u.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    nd >= k1*idim
c    t(k1) <= u <= t(n-k1+1)
c
c  further comments:
c    if u coincides with a knot, right derivatives are computed
c    ( left derivatives if u = t(n-k1+1) ).
c
c  other subroutines required: fpader.
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
      integer idim,n,nc,k1,nd,ier
      real*8 u
c  ..array arguments..
      real*8 t(n),c(nc),d(nd)
c  ..local scalars..
      integer i,j,kk,l,m,nk1
c  ..local array..
      real*8 h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nd.lt.(k1*idim)) go to 500
      nk1 = n-k1
      if(u.lt.t(k1) .or. u.gt.t(nk1+1)) go to 500
c  search for knot interval t(l) <= u < t(l+1)
      l = k1
 100  if(u.lt.t(l+1) .or. l.eq.nk1) go to 200
      l = l+1
      go to 100
 200  if(t(l).ge.t(l+1)) go to 500
      ier = 0
c  calculate the derivatives.
      j = 1
      do 400 i=1,idim
        call fpader(t,n,c(j),k1,u,l,h)
        m = i
        do 300 kk=1,k1
          d(m) = h(kk)
          m = m+idim
 300    continue
        j = j+n
 400  continue
 500  return
      end
