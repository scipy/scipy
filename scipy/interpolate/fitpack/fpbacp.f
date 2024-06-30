      recursive subroutine fpbacp(a,b,z,n,k,c,k1,nest)
      implicit none
c  subroutine fpbacp calculates the solution of the system of equations
c  g * c = z  with g  a n x n upper triangular matrix of the form
c            ! a '   !
c        g = !   ' b !
c            ! 0 '   !
c  with b a n x k matrix and a a (n-k) x (n-k) upper triangular
c  matrix of bandwidth k1.
c  ..
c  ..scalar arguments..
      integer n,k,k1,nest
c  ..array arguments..
      real*8 a(nest,k1),b(nest,k),z(n),c(n)
c  ..local scalars..
      integer i,i1,j,l,l0,l1,n2
      real*8 store
c  ..
      n2 = n-k
      l = n
      do 30 i=1,k
        store = z(l)
        j = k+2-i
        if(i.eq.1) go to 20
        l0 = l
        do 10 l1=j,k
          l0 = l0+1
          store = store-c(l0)*b(l,l1)
  10    continue
  20    c(l) = store/b(l,j-1)
        l = l-1
        if(l.eq.0) go to 80
  30  continue
      do 50 i=1,n2
        store = z(i)
        l = n2
        do 40 j=1,k
          l = l+1
          store = store-c(l)*b(i,j)
  40    continue
        c(i) = store
  50  continue
      i = n2
      c(i) = c(i)/a(i,1)
      if(i.eq.1) go to 80
      do 70 j=2,n2
        i = i-1
        store = c(i)
        i1 = k
        if(j.le.k) i1=j-1
        l = i
        do 60 l0=1,i1
          l = l+1
          store = store-c(l)*a(i,l0+1)
  60    continue
        c(i) = store/a(i,1)
  70  continue
  80  return
      end
