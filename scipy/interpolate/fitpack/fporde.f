      recursive subroutine fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,
     *   index,nreg)
c  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
c  according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
c  to. for each panel a stack is constructed  containing the numbers
c  of data points lying inside; index(j),j=1,2,...,nreg points to the
c  first data point in the jth panel while nummer(i),i=1,2,...,m gives
c  the number of the next data point in the panel.
c  ..
c  ..scalar arguments..
      integer m,kx,ky,nx,ny,nreg
c  ..array arguments..
      real*8 x(m),y(m),tx(nx),ty(ny)
      integer nummer(m),index(nreg)
c  ..local scalars..
      real*8 xi,yi
      integer i,im,k,kx1,ky1,k1,l,l1,nk1x,nk1y,num,nyy
c  ..
      kx1 = kx+1
      ky1 = ky+1
      nk1x = nx-kx1
      nk1y = ny-ky1
      nyy = nk1y-ky
      do 10 i=1,nreg
        index(i) = 0
  10  continue
      do 60 im=1,m
        xi = x(im)
        yi = y(im)
        l = kx1
        l1 = l+1
  20    if(xi.lt.tx(l1) .or. l.eq.nk1x) go to 30
        l = l1
        l1 = l+1
        go to 20
  30    k = ky1
        k1 = k+1
  40    if(yi.lt.ty(k1) .or. k.eq.nk1y) go to 50
        k = k1
        k1 = k+1
        go to 40
  50    num = (l-kx1)*nyy+k-ky
        nummer(im) = index(num)
        index(num) = im
  60  continue
      return
      end

