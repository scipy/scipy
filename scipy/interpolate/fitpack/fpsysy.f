      subroutine fpsysy(a,n,g)
c subroutine fpsysy solves a linear n x n symmetric system
c    (a) * (b) = (g)
c on input, vector g contains the right hand side ; on output it will
c contain the solution (b).
c  ..
c  ..scalar arguments..
      integer n
c  ..array arguments..
      real*8 a(6,6),g(6)
c  ..local scalars..
      real*8 fac
      integer i,i1,j,k
c  ..
      g(1) = g(1)/a(1,1)
      if(n.eq.1) return
c  decomposition of the symmetric matrix (a) = (l) * (d) *(l)'
c  with (l) a unit lower triangular matrix and (d) a diagonal
c  matrix
      do 10 k=2,n
         a(k,1) = a(k,1)/a(1,1)
  10  continue
      do 40 i=2,n
         i1 = i-1
         do 30 k=i,n
            fac = a(k,i)
            do 20 j=1,i1
               fac = fac-a(j,j)*a(k,j)*a(i,j)
  20        continue
            a(k,i) = fac
            if(k.gt.i) a(k,i) = fac/a(i,i)
  30     continue
  40  continue
c  solve the system (l)*(d)*(l)'*(b) = (g).
c  first step : solve (l)*(d)*(c) = (g).
      do 60 i=2,n
         i1 = i-1
         fac = g(i)
         do 50 j=1,i1
            fac = fac-g(j)*a(j,j)*a(i,j)
  50     continue
         g(i) = fac/a(i,i)
  60  continue
c  second step : solve (l)'*(b) = (c)
      i = n
      do 80 j=2,n
         i1 = i
         i = i-1
         fac = g(i)
         do 70 k=i1,n
            fac = fac-g(k)*a(k,i)
  70     continue
         g(i) = fac
  80  continue
      return
      end
