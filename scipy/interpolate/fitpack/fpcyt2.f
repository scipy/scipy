      subroutine fpcyt2(a,n,b,c,nn)
c subroutine fpcyt2 solves a linear n x n system
c         a * c = b
c where matrix a is a cyclic tridiagonal matrix, decomposed
c using subroutine fpsyt1.
c  ..
c  ..scalar arguments..
      integer n,nn
c  ..array arguments..
      real*8 a(nn,6),b(n),c(n)
c  ..local scalars..
      real*8 cc,sum
      integer i,j,j1,n1
c  ..
      c(1) = b(1)*a(1,4)
      sum = c(1)*a(1,5)
      n1 = n-1
      do 10 i=2,n1
         c(i) = (b(i)-a(i,1)*c(i-1))*a(i,4)
         sum = sum+c(i)*a(i,5)
  10  continue
      cc = (b(n)-sum)*a(n,4)
      c(n) = cc
      c(n1) = c(n1)-cc*a(n1,6)
      j = n1
      do 20 i=3,n
         j1 = j-1
         c(j1) = c(j1)-c(j)*a(j1,3)*a(j1,4)-cc*a(j1,6)
         j = j1
  20  continue
      return
      end
