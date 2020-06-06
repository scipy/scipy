      subroutine pardeu(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,y,z,m,
     * wrk,lwrk,iwrk,kwrk,ier)
c  subroutine pardeu evaluates on a set of points (x(i),y(i)),i=1,...,m
c  the partial derivative ( order nux,nuy) of a bivariate spline
c  s(x,y) of degrees kx and ky, given in the b-spline representation.
c
c  calling sequence:
c     call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,
c    * iwrk,kwrk,ier)
c
c  input parameters:
c   tx    : real array, length nx, which contains the position of the
c           knots in the x-direction.
c   nx    : integer, giving the total number of knots in the x-direction
c   ty    : real array, length ny, which contains the position of the
c           knots in the y-direction.
c   ny    : integer, giving the total number of knots in the y-direction
c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients.
c   kx,ky : integer values, giving the degrees of the spline.
c   nux   : integer values, specifying the order of the partial
c   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
c   kx,ky : integer values, giving the degrees of the spline.
c   x     : real array of dimension (mx).
c   y     : real array of dimension (my).
c   m     : on entry m must specify the number points. m >= 1.
c   wrk   : real array of dimension lwrk. used as workspace.
c   lwrk  : integer, specifying the dimension of wrk.
c           lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1)
c   iwrk  : integer array of dimension kwrk. used as workspace.
c   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
c
c  output parameters:
c   z     : real array of dimension (m).
c           on successful exit z(i) contains the value of the
c           specified partial derivative of s(x,y) at the point
c           (x(i),y(i)),i=1,...,m.
c   ier   : integer error flag
c   ier=0 : normal return
c   ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   lwrk>=m*(kx+1-nux)+m*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1),
c
c  other subroutines required:
c    fpbisp,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1989
c
c  ..scalar arguments..
      integer nx,ny,kx,ky,m,lwrk,kwrk,ier,nux,nuy
c  ..array arguments..
      integer iwrk(kwrk)
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(m),y(m),z(m),
     * wrk(lwrk)
c  ..local scalars..
      integer i,iwx,iwy,j,kkx,kky,kx1,ky1,lx,ly,lwest,l1,l2,mm,m0,m1,
     * nc,nkx1,nky1,nxx,nyy
      real*8 ak,fac
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      nc = nkx1*nky1
      if(nux.lt.0 .or. nux.ge.kx) go to 400
      if(nuy.lt.0 .or. nuy.ge.ky) go to 400
      lwest = nc +(kx1-nux)*m+(ky1-nuy)*m
      if(lwrk.lt.lwest) go to 400
      if(kwrk.lt.(m+m)) go to 400
      if (m.lt.1) go to 400
      ier = 0
      nxx = nkx1
      nyy = nky1
      kkx = kx
      kky = ky
c  the partial derivative of order (nux,nuy) of a bivariate spline of
c  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
c  we calculate the b-spline coefficients of this spline
      do 70 i=1,nc
        wrk(i) = c(i)
  70  continue
      if(nux.eq.0) go to 200
      lx = 1
      do 100 j=1,nux
        ak = kkx
        nxx = nxx-1
        l1 = lx
        m0 = 1
        do 90 i=1,nxx
          l1 = l1+1
          l2 = l1+kkx
          fac = tx(l2)-tx(l1)
          if(fac.le.0.) go to 90
          do 80 mm=1,nyy
            m1 = m0+nyy
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+1
  80      continue
  90    continue
        lx = lx+1
        kkx = kkx-1
 100  continue
 200  if(nuy.eq.0) go to 300
      ly = 1
      do 230 j=1,nuy
        ak = kky
        nyy = nyy-1
        l1 = ly
        do 220 i=1,nyy
          l1 = l1+1
          l2 = l1+kky
          fac = ty(l2)-ty(l1)
          if(fac.le.0.) go to 220
          m0 = i
          do 210 mm=1,nxx
            m1 = m0+1
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+nky1
 210      continue
 220    continue
        ly = ly+1
        kky = kky-1
 230  continue
      m0 = nyy
      m1 = nky1
      do 250 mm=2,nxx
        do 240 i=1,nyy
          m0 = m0+1
          m1 = m1+1
          wrk(m0) = wrk(m1)
 240    continue
        m1 = m1+nuy
 250  continue
c  we partition the working space and evaluate the partial derivative
 300  iwx = 1+nxx*nyy
      iwy = iwx+m*(kx1-nux)
      do 390 i=1,m
        call fpbisp(tx(nux+1),nx-2*nux,ty(nuy+1),ny-2*nuy,wrk,kkx,kky,
     *   x(i),1,y(i),1,z(i),wrk(iwx),wrk(iwy),iwrk(1),iwrk(2))
 390  continue
 400  return
      end
