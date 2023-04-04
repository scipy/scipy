      recursive subroutine pardtc(tx,nx,ty,ny,c,kx,ky,nux,nuy,
     *   newc,ier)
      implicit none
c  subroutine pardtc takes the knots and coefficients of a bivariate
c  spline, and returns the coefficients for a new bivariate spline that
c  evaluates the partial derivative (order nux, nuy) of the original
c  spline.
c
c  calling sequence:
c     call pardtc(tx,nx,ty,ny,c,kx,ky,nux,nuy,newc,ier)
c
c  input parameters:
c   tx    : real array, length nx, which contains the position of the
c           knots in the x-direction.
c   nx    : integer, giving the total number of knots in the x-direction
c           (hidden)
c   ty    : real array, length ny, which contains the position of the
c           knots in the y-direction.
c   ny    : integer, giving the total number of knots in the y-direction
c           (hidden)
c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients.
c   kx,ky : integer values, giving the degrees of the spline.
c   nux   : integer values, specifying the order of the partial
c   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
c
c  output parameters:
c   newc  : real array containing the coefficients of the derivative.
c           the dimension is (nx-nux-kx-1)*(ny-nuy-ky-1).
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   0 <= nux < kx, 0 <= nuy < kyc
c
c  other subroutines required:
c    none
c
c  references :
c   de boor c  : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  based on the subroutine "parder" by Paul Dierckx.
c
c  author :
c    Cong Ma
c    Department of Mathematics and Applied Mathematics, U. of Cape Town
c    Cross Campus Road, Rondebosch 7700, Cape Town, South Africa.
c    e-mail : cong.ma@uct.ac.za
c
c  latest update : may 2019
c
c  ..scalar arguments..
      integer nx,ny,kx,ky,nux,nuy,ier, nc
c  ..array arguments..
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),
     * newc((nx-kx-1)*(ny-ky-1))
c  ..local scalars..
      integer i,j,kx1,ky1,lx,ly,l1,l2,m,m0,m1,
     * nkx1,nky1,nxx,nyy,newkx,newky
      real*8 ak,fac
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nux.lt.0 .or. nux.ge.kx) go to 400
      if(nuy.lt.0 .or. nuy.ge.ky) go to 400
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      nc = nkx1*nky1
      ier = 0
      nxx = nkx1
      nyy = nky1
      newkx = kx
      newky = ky
c  the partial derivative of order (nux,nuy) of a bivariate spline of
c  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
c  we calculate the b-spline coefficients of this spline
c  that is to say newkx = kx - nux, newky = ky - nuy
      do 70 i=1,nc
        newc(i) = c(i)
  70  continue
      if(nux.eq.0) go to 200
      lx = 1
      do 100 j=1,nux
        ak = newkx
        nxx = nxx-1
        l1 = lx
        m0 = 1
        do 90 i=1,nxx
          l1 = l1+1
          l2 = l1+newkx
          fac = tx(l2)-tx(l1)
          if(fac.le.0.) go to 90
          do 80 m=1,nyy
            m1 = m0+nyy
            newc(m0) = (newc(m1)-newc(m0))*ak/fac
            m0  = m0+1
  80      continue
  90    continue
        lx = lx+1
        newkx = newkx-1
 100  continue
 200  if(nuy.eq.0) go to 400
c orig: if(nuy.eq.0) go to 300
      ly = 1
      do 230 j=1,nuy
        ak = newky
        nyy = nyy-1
        l1 = ly
        do 220 i=1,nyy
          l1 = l1+1
          l2 = l1+newky
          fac = ty(l2)-ty(l1)
          if(fac.le.0.) go to 220
          m0 = i
          do 210 m=1,nxx
            m1 = m0+1
            newc(m0) = (newc(m1)-newc(m0))*ak/fac
            m0  = m0+nky1
 210      continue
 220    continue
        ly = ly+1
        newky = newky-1
 230  continue
      m0 = nyy
      m1 = nky1
      do 250 m=2,nxx
        do 240 i=1,nyy
          m0 = m0+1
          m1 = m1+1
          newc(m0) = newc(m1)
 240    continue
        m1 = m1+nuy
 250  continue
c300  iwx = 1+nxx*nyy
c     iwy = iwx+mx*(kx1-nux)
c
c from parder.f:
c     call fpbisp(tx(nux+1),nx-2*nux,ty(nuy+1),ny-2*nuy,newc,newkx,newky,
c    * x,mx,y,my,z,newc(iwx),newc(iwy),iwrk(1),iwrk(mx+1))
c
c from bispev.f:
c     call fpbisp(tx,       nx,      ty,       ny,      c,   kx,   ky,
c    * x,mx,y,my,z,wrk(1),   wrk(iw),  iwrk(1),iwrk(mx+1))
c
c from fpbisp.f:
c          fpbisp(tx,       nx,      ty,       ny,      c,   kx,   ky,
c    * x,mx,y,my,z,wx,       wy,       lx,     ly)
 400  return
      end

