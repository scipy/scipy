      real*8 function dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk)
c  function dblint calculates the double integral
c         / xe  / ye
c        |     |      s(x,y) dx dy
c    xb /  yb /
c  with s(x,y) a bivariate spline of degrees kx and ky, given in the
c  b-spline representation.
c
c  calling sequence:
c     aint = dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk)
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
c   xb,xe : real values, containing the boundaries of the integration
c   yb,ye   domain. s(x,y) is considered to be identically zero out-
c           side the rectangle (tx(kx+1),tx(nx-kx))*(ty(ky+1),ty(ny-ky))
c
c  output parameters:
c   aint  : real , containing the double integral of s(x,y).
c   wrk   : real array of dimension at least (nx+ny-kx-ky-2).
c           used as working space.
c           on exit, wrk(i) will contain the integral
c                / xe
c               | ni,kx+1(x) dx , i=1,2,...,nx-kx-1
c           xb /
c           with ni,kx+1(x) the normalized b-spline defined on
c           the knots tx(i),...,tx(i+kx+1)
c           wrk(j+nx-kx-1) will contain the integral
c                / ye
c               | nj,ky+1(y) dy , j=1,2,...,ny-ky-1
c           yb /
c           with nj,ky+1(y) the normalized b-spline defined on
c           the knots ty(j),...,ty(j+ky+1)
c
c  other subroutines required: fpintb
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
c  latest update : march 1989
c
c  ..scalar arguments..
      integer nx,ny,kx,ky
      real*8 xb,xe,yb,ye
c  ..array arguments..
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),wrk(nx+ny-kx-ky-2)
c  ..local scalars..
      integer i,j,l,m,nkx1,nky1
      real*8 res
c  ..
      nkx1 = nx-kx-1
      nky1 = ny-ky-1
c  we calculate the integrals of the normalized b-splines ni,kx+1(x)
      call fpintb(tx,nx,wrk,nkx1,xb,xe)
c  we calculate the integrals of the normalized b-splines nj,ky+1(y)
      call fpintb(ty,ny,wrk(nkx1+1),nky1,yb,ye)
c  calculate the integral of s(x,y)
      dblint = 0.
      do 200 i=1,nkx1
        res = wrk(i)
        if(res.eq.0.) go to 200
        m = (i-1)*nky1
        l = nkx1
        do 100 j=1,nky1
          m = m+1
          l = l+1
          dblint = dblint+res*wrk(l)*c(m)
 100    continue
 200  continue
      return
      end
