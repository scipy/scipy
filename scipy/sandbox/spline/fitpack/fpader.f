      subroutine fpader(t,n,c,k1,x,l,d)
c  subroutine fpader calculates the derivatives
c             (j-1)
c     d(j) = s     (x) , j=1,2,...,k1
c  of a spline of order k1 at the point t(l)<=x<t(l+1), using the
c  stable recurrence scheme of de boor
c  ..
c  ..scalar arguments..
      real*8 x
      integer n,k1,l
c  ..array arguments..
      real*8 t(n),c(n),d(k1)
c  ..local scalars..
      integer i,ik,j,jj,j1,j2,ki,kj,li,lj,lk
      real*8 ak,fac,one
c  ..local array..
      real*8 h(6)
c  ..
      one = 0.1d+01
      lk = l-k1
      do 100 i=1,k1
        ik = i+lk
        h(i) = c(ik)
 100  continue
      kj = k1
      fac = one
      do 700 j=1,k1
        ki = kj
        j1 = j+1
        if(j.eq.1) go to 300
        i = k1
        do 200 jj=j,k1
          li = i+lk
          lj = li+kj
          h(i) = (h(i)-h(i-1))/(t(lj)-t(li))
          i = i-1
 200    continue
 300    do 400 i=j,k1
          d(i) = h(i)
 400    continue
        if(j.eq.k1) go to 600
        do 500 jj=j1,k1
          ki = ki-1
          i = k1
          do 500 j2=jj,k1
            li = i+lk
            lj = li+ki
            d(i) = ((x-t(li))*d(i)+(t(lj)-x)*d(i-1))/(t(lj)-t(li))
            i = i-1
 500    continue
 600    d(j) = d(k1)*fac
        ak = k1-j
        fac = fac*ak
        kj = kj-1
 700  continue
      return
      end
