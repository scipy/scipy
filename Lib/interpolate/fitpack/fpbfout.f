      subroutine fpbfou(t,n,par,ress,resc)
c  subroutine fpbfou calculates the integrals
c                    /t(n-3)
c    ress(j) =      !        nj,4(x)*sin(par*x) dx    and
c              t(4)/
c                    /t(n-3)
c    resc(j) =      !        nj,4(x)*cos(par*x) dx ,  j=1,2,...n-4
c              t(4)/
c  where nj,4(x) denotes the cubic b-spline defined on the knots
c  t(j),t(j+1),...,t(j+4).
c
c  calling sequence:
c     call fpbfou(t,n,par,ress,resc)
c
c  input parameters:
c    t    : real array,length n, containing the knots.
c    n    : integer, containing the number of knots.
c    par  : real, containing the value of the parameter par.
c
c  output parameters:
c    ress  : real array,length n, containing the integrals ress(j).
c    resc  : real array,length n, containing the integrals resc(j).
c
c  restrictions:
c    n >= 10, t(4) < t(5) < ... < t(n-4) < t(n-3).
c  ..
c  ..scalar arguments..
      integer n
      real*8 par
c  ..array arguments..
      real*8 t(n),ress(n),resc(n)
c  ..local scalars..
      integer i,ic,ipj,is,j,jj,jp1,jp4,k,li,lj,ll,nmj,nm3,nm7
      real*8 ak,beta,con1,con2,c1,c2,delta,eps,fac,f1,f2,f3,one,quart,
     * sign,six,s1,s2,term
c  ..local arrays..
      real*8 co(5),si(5),hs(5),hc(5),rs(3),rc(3)
c  ..function references..
      real*8 cos,sin,abs
c  ..
c  initialization.
      one = 0.1e+01
      six = 0.6e+01
      eps = 0.1e-07
      quart = 0.25e0
      con1 = 0.5e-01
      con2 = 0.12e+03
      nm3 = n-3
      nm7 = n-7
      if(par.ne.0.) term = six/par
      beta = par*t(4)
      co(1) = cos(beta)
      si(1) = sin(beta)
c  calculate the integrals ress(j) and resc(j), j=1,2,3 by setting up
c  a divided difference table.
      do 30 j=1,3
        jp1 = j+1
        jp4 = j+4
        beta = par*t(jp4)
        co(jp1) = cos(beta)
        si(jp1) = sin(beta)
        call fpcsin(t(4),t(jp4),par,si(1),co(1),si(jp1),co(jp1),
     *  rs(j),rc(j))
        i = 5-j
        hs(i) = 0.
        hc(i) = 0.
        do 10 jj=1,j
          ipj = i+jj
          hs(ipj) = rs(jj)
          hc(ipj) = rc(jj)
  10    continue
        do 20 jj=1,3
          if(i.lt.jj) i = jj
          k = 5
          li = jp4
          do 20 ll=i,4
            lj = li-jj
            fac = t(li)-t(lj)
            hs(k) = (hs(k)-hs(k-1))/fac
            hc(k) = (hc(k)-hc(k-1))/fac
            k = k-1
            li = li-1
  20    continue
        ress(j) = hs(5)-hs(4)
        resc(j) = hc(5)-hc(4)
  30  continue
      if(nm7.lt.4) go to 160
c  calculate the integrals ress(j) and resc(j),j=4,5,...,n-7.
      do 150 j=4,nm7
        jp4 = j+4
        beta = par*t(jp4)
        co(5) = cos(beta)
        si(5) = sin(beta)
        delta = t(jp4)-t(j)
c  the way of computing ress(j) and resc(j) depends on the value of
c  beta = par*(t(j+4)-t(j)).
        beta = delta*par
        if(abs(beta).le.one) go to 60
c  if !beta! > 1 the integrals are calculated by setting up a divided
c  difference table.
        do 40 k=1,5
          hs(k) = si(k)
          hc(k) = co(k)
  40    continue
        do 50 jj=1,3
          k = 5
          li = jp4
          do 50 ll=jj,4
            lj = li-jj
            fac = par*(t(li)-t(lj))
            hs(k) = (hs(k)-hs(k-1))/fac
            hc(k) = (hc(k)-hc(k-1))/fac
            k = k-1
            li = li-1
  50    continue
        s2 = (hs(5)-hs(4))*term
        c2 = (hc(5)-hc(4))*term
        go to 130
c  if !beta! <= 1 the integrals are calculated by evaluating a series
c  expansion.
  60    f3 = 0.
        do 70 i=1,4
          ipj = i+j
          hs(i) = par*(t(ipj)-t(j))
          hc(i) = hs(i)
          f3 = f3+hs(i)
  70    continue
        f3 = f3*con1
        c1 = quart
        s1 = f3
        if(abs(f3).le.eps) go to 120
        sign = one
        fac = con2
        k = 5
        is = 0
        do 110 ic=1,20
          k = k+1
          ak = k
          fac = fac*ak
          f1 = 0.
          f3 = 0.
          do 80 i=1,4
            f1 = f1+hc(i)
            f2 = f1*hs(i)
            hc(i) = f2
            f3 = f3+f2
  80      continue
          f3 = f3*six/fac
          if(is.eq.0) go to 90
          is = 0
          s1 = s1+f3*sign
          go to 100
  90      sign = -sign
          is = 1
          c1 = c1+f3*sign
 100      if(abs(f3).le.eps) go to 120
 110    continue
 120    s2 = delta*(co(1)*s1+si(1)*c1)
        c2 = delta*(co(1)*c1-si(1)*s1)
 130    ress(j) = s2
        resc(j) = c2
        do 140 i=1,4
          co(i) = co(i+1)
          si(i) = si(i+1)
 140    continue
 150  continue
c  calculate the integrals ress(j) and resc(j),j=n-6,n-5,n-4 by setting
c  up a divided difference table.
 160  do 190 j=1,3
        nmj = nm3-j
        i = 5-j
        call fpcsin(t(nm3),t(nmj),par,si(4),co(4),si(i-1),co(i-1),
     *  rs(j),rc(j))
        hs(i) = 0.
        hc(i) = 0.
        do 170 jj=1,j
          ipj = i+jj
          hc(ipj) = rc(jj)
          hs(ipj) = rs(jj)
 170    continue
        do 180 jj=1,3
          if(i.lt.jj) i = jj
          k = 5
          li = nmj
          do 180 ll=i,4
            lj = li+jj
            fac = t(lj)-t(li)
            hs(k) = (hs(k-1)-hs(k))/fac
            hc(k) = (hc(k-1)-hc(k))/fac
            k = k-1
            li = li+1
 180    continue
        ress(nmj) = hs(4)-hs(5)
        resc(nmj) = hc(4)-hc(5)
 190  continue
      return
      end
