      recursive subroutine fprpsp(nt,np,co,si,c,f,ncoff)
      implicit none
c  given the coefficients of a spherical spline function, subroutine
c  fprpsp calculates the coefficients in the standard b-spline re-
c  presentation of this bicubic spline.
c  ..
c  ..scalar arguments
      integer nt,np,ncoff
c  ..array arguments
      real*8 co(np),si(np),c(ncoff),f(ncoff)
c  ..local scalars
      real*8 cn,c1,c2,c3
      integer i,ii,j,k,l,ncof,npp,np4,nt4
c  ..
      nt4 = nt-4
      np4 = np-4
      npp = np4-3
      ncof = 6+npp*(nt4-4)
      c1 = c(1)
      cn = c(ncof)
      j = ncoff
      do 10 i=1,np4
         f(i) = c1
         f(j) = cn
         j = j-1
  10  continue
      i = np4
      j=1
      do 70 l=3,nt4
         ii = i
         if(l.eq.3 .or. l.eq.nt4) go to 30
         do 20 k=1,npp
            i = i+1
            j = j+1
            f(i) = c(j)
  20     continue
         go to 50
  30     if(l.eq.nt4) c1 = cn
         c2 = c(j+1)
         c3 = c(j+2)
         j = j+2
         do 40 k=1,npp
            i = i+1
            f(i) = c1+c2*co(k)+c3*si(k)
  40     continue
  50     do 60 k=1,3
            ii = ii+1
            i = i+1
            f(i) = f(ii)
  60     continue
  70  continue
      do 80 i=1,ncoff
         c(i) = f(i)
  80  continue
      return
      end
