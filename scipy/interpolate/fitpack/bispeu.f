      subroutine bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk, ier)
c  subroutine bispeu evaluates on a set of points (x(i),y(i)),i=1,...,m
c  a bivariate spline s(x,y) of degrees kx and ky, given in the
c  b-spline representation.
c
c  calling sequence:
c     call bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk,
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
c   x     : real array of dimension (mx).
c   y     : real array of dimension (my).
c   m     : on entry m must specify the number points. m >= 1.
c   wrk   : real array of dimension lwrk. used as workspace.
c   lwrk  : integer, specifying the dimension of wrk.
c           lwrk >= kx+ky+2
c
c  output parameters:
c   z     : real array of dimension m.
c           on succesful exit z(i) contains the value of s(x,y)
c           at the point (x(i),y(i)), i=1,...,m.
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   m >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
c   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
c   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
c
c  other subroutines required:
c    fpbisp,fpbspl
c
c  ..scalar arguments..
      integer nx,ny,kx,ky,m,lwrk,kwrk,ier
c  ..array arguments..
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(m),y(m),z(m),
     *     wrk(lwrk)
c  ..local scalars..
      integer iwrk(2)
      integer i,iw,lwest
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      lwest = kx+ky+2
      if (lwrk.lt.lwest) go to 100
      if (m.lt.1) go to 100
      ier = 0
      do 10 i=1,m
         call fpbisp(tx,nx,ty,ny,c,kx,ky,x(i),1,y(i),1,z(i),wrk(1),
     *        wrk(kx+2),iwrk(1),iwrk(2))
 10   continue
 100  return
      end
