      recursive subroutine fprppo(nu,nv,if1,if2,cosi,ratio,c,f,ncoff)
      implicit none
c  given the coefficients of a constrained bicubic spline, as determined
c  in subroutine fppola, subroutine fprppo calculates the coefficients
c  in the standard b-spline representation of bicubic splines.
c  ..
c  ..scalar arguments..
      real*8 ratio
      integer nu,nv,if1,if2,ncoff
c  ..array arguments
      real*8 c(ncoff),f(ncoff),cosi(5,nv)
c  ..local scalars..
      integer i,iopt,ii,j,k,l,nu4,nvv
c  ..
      nu4 = nu-4
      nvv = nv-7
      iopt = if1+1
      do 10 i=1,ncoff
         f(i) = 0.
  10  continue
      i = 0
      do 120 l=1,nu4
         ii = i
         if(l.gt.iopt) go to 80
         go to (20,40,60),l
  20     do 30 k=1,nvv
            i = i+1
            f(i) = c(1)
  30     continue
         j = 1
         go to 100
  40     do 50 k=1,nvv
            i = i+1
            f(i) = c(1)+c(2)*cosi(1,k)+c(3)*cosi(2,k)
  50     continue
         j = 3
         go to 100
  60     do 70 k=1,nvv
            i = i+1
            f(i) = c(1)+ratio*(c(2)*cosi(1,k)+c(3)*cosi(2,k))+
     *             c(4)*cosi(3,k)+c(5)*cosi(4,k)+c(6)*cosi(5,k)
  70     continue
         j = 6
         go to 100
  80     if(l.eq.nu4 .and. if2.ne.0) go to 120
         do 90 k=1,nvv
            i = i+1
            j = j+1
            f(i) = c(j)
  90     continue
 100     do 110 k=1,3
            ii = ii+1
            i = i+1
            f(i) = f(ii)
 110     continue
 120  continue
      do 130 i=1,ncoff
         c(i) = f(i)
 130  continue
      return
      end

