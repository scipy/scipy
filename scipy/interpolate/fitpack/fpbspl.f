      recursive subroutine fpbspl(t,n,k,x,l,h)
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  Travis Oliphant  2007
c    changed so that weighting of 0 is used when knots with
c      multiplicity are present.
c    Also, notice that l+k <= n and 1 <= l+1-k
c      or else the routine will be accessing memory outside t
c      Thus it is imperative that that k <= l <= n-k but this
c      is not checked.
c  ..
c  ..scalar arguments..
      real*8 x
      integer n,k,l
c  ..array arguments..
      real*8 t(n),h(20)
c  ..local scalars..
      real*8 f,one
      integer i,j,li,lj
c  ..local arrays..
      real*8 hh(19)
c  ..
      one = 0.1d+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.0d0
        do 20 i=1,j
          li = l+i
          lj = li-j
          if (t(li).ne.t(lj)) goto 15
          h(i+1) = 0.0d0 
          goto 20
  15      f = hh(i)/(t(li)-t(lj)) 
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end
