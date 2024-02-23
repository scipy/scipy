      recursive subroutine fpchec(x,m,t,n,k,ier)
      implicit none
c  subroutine fpchec verifies the number and the position of the knots
c  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
c  and the position of the data points x(i),i=1,2,...,m. if all of the
c  following conditions are fulfilled, the error parameter ier is set
c  to zero. if one of the conditions is violated ier is set to ten.
c      1) k+1 <= n-k-1 <= m
c      2) t(1) <= t(2) <= ... <= t(k+1)
c         t(n-k) <= t(n-k+1) <= ... <= t(n)
c      3) t(k+1) < t(k+2) < ... < t(n-k)
c      4) t(k+1) <= x(i) <= t(n-k)
c      5) the conditions specified by schoenberg and whitney must hold
c         for at least one subset of data points, i.e. there must be a
c         subset of data points y(j) such that
c             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
c  ..
c  ..scalar arguments..
      integer m,n,k,ier
c  ..array arguments..
      real*8 x(m),t(n)
c  ..local scalars..
      integer i,j,k1,k2,l,nk1,nk2,nk3
      real*8 tj,tl
c  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      ier = 10
c  check condition no 1
      if(nk1.lt.k1 .or. nk1.gt.m)then; ier=10; go to 80; endif
c  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1))then; ier=20; go to 80; endif
        if(t(j).lt.t(j-1))then; ier=20; go to 80; endif
        j = j-1
  20  continue
c  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1))then; ier=30; go to 80; endif
  30  continue
c  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2))then; ier=40; go to 80;
      endif
c  check condition no 5
      if(x(1).ge.t(k2) .or. x(m).le.t(nk1))then; ier=50; go to 80;
      endif
      i = 1
      l = k2
      nk3 = nk1-1
      if(nk3.lt.2) go to 70
      do 60 j=2,nk3
        tj = t(j)
        l = l+1
        tl = t(l)
  40    i = i+1
        if(i.ge.m)then; ier=50; go to 80; endif
        if(x(i).le.tj) go to 40
        if(x(i).ge.tl)then; ier=50; go to 80; endif
  60  continue
  70  ier = 0
  80  return
      end
