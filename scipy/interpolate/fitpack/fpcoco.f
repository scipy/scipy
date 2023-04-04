      recursive subroutine fpcoco(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,
     *   n,t,c,sq,sx,bind,e,wrk,lwrk,iwrk,kwrk,ier)
      implicit none
c  ..scalar arguments..
      real*8 s,sq
      integer iopt,m,nest,maxtr,maxbin,n,lwrk,kwrk,ier
c  ..array arguments..
      integer iwrk(kwrk)
      real*8 x(m),y(m),w(m),v(m),t(nest),c(nest),sx(m),e(nest),wrk(lwrk)
     *
      logical bind(nest)
c  ..local scalars..
      integer i,ia,ib,ic,iq,iu,iz,izz,i1,j,k,l,l1,m1,nmax,nr,n4,n6,n8,
     * ji,jib,jjb,jl,jr,ju,mb,nm
      real*8 sql,sqmax,term,tj,xi,half
c  ..subroutine references..
c    fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno
c  ..
c  set constant
      half = 0.5e0
c  determine the maximal admissible number of knots.
      nmax = m+4
c  the initial choice of knots depends on the value of iopt.
c    if iopt=0 the program starts with the minimal number of knots
c    so that can be guarantied that the concavity/convexity constraints
c    will be satisfied.
c    if iopt = 1 the program will continue from the point on where she
c    left at the foregoing call.
      if(iopt.gt.0) go to 80
c  find the minimal number of knots.
c  a knot is located at the data point x(i), i=2,3,...m-1 if
c    1) v(i) ^= 0    and
c    2) v(i)*v(i-1) <= 0  or  v(i)*v(i+1) <= 0.
      m1 = m-1
      n = 4
      do 20 i=2,m1
        if(v(i).eq.0. .or. (v(i)*v(i-1).gt.0. .and.
     *  v(i)*v(i+1).gt.0.)) go to 20
        n = n+1
c  test whether the required storage space exceeds the available one.
        if(n+4.gt.nest) go to 200
        t(n) = x(i)
  20  continue
c  find the position of the knots t(1),...t(4) and t(n-3),...t(n) which
c  are needed for the b-spline representation of s(x).
      do 30 i=1,4
        t(i) = x(1)
        n = n+1
        t(n) = x(m)
  30  continue
c  test whether the minimum number of knots exceeds the maximum number.
      if(n.gt.nmax) go to 210
c  main loop for the different sets of knots.
c  find corresponding values e(j) to the knots t(j+3),j=1,2,...n-6
c    e(j) will take the value -1,1, or 0 according to the requirement
c    that s(x) must be locally convex or concave at t(j+3) or that the
c    sign of s''(x) is unrestricted at that point.
  40  i= 1
      xi = x(1)
      j = 4
      tj = t(4)
      n6 = n-6
      do 70 l=1,n6
  50    if(xi.eq.tj) go to 60
        i = i+1
        xi = x(i)
        go to 50
  60    e(l) = v(i)
        j = j+1
        tj = t(j)
  70  continue
c  we partition the working space
      nm = n+maxbin
      mb = maxbin+1
      ia = 1
      ib = ia+4*n
      ic = ib+nm*maxbin
      iz = ic+n
      izz = iz+n
      iu = izz+n
      iq = iu+maxbin
      ji = 1
      ju = ji+maxtr
      jl = ju+maxtr
      jr = jl+maxtr
      jjb = jr+maxtr
      jib = jjb+mb
c  given the set of knots t(j),j=1,2,...n, find the least-squares cubic
c  spline which satisfies the imposed concavity/convexity constraints.
      call fpcosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,nm,mb,wrk(ia),
     *
     * wrk(ib),wrk(ic),wrk(iz),wrk(izz),wrk(iu),wrk(iq),iwrk(ji),
     * iwrk(ju),iwrk(jl),iwrk(jr),iwrk(jjb),iwrk(jib),ier)
c  if sq <= s or in case of abnormal exit from fpcosp, control is
c  repassed to the driver program.
      if(sq.le.s .or. ier.gt.0) go to 300
c  calculate for each knot interval t(l-1) <= xi <= t(l) the
c  sum((wi*(yi-s(xi)))**2).
c  find the interval t(k-1) <= x <= t(k) for which this sum is maximal
c  on the condition that this interval contains at least one interior
c  data point x(nr) and that s(x) is not given there by a straight line.
  80  sqmax = 0.
      sql = 0.
      l = 5
      nr = 0
      i1 = 1
      n4 = n-4
      do 110 i=1,m
        term = (w(i)*(sx(i)-y(i)))**2
        if(x(i).lt.t(l) .or. l.gt.n4) go to 100
        term = term*half
        sql = sql+term
        if(i-i1.le.1 .or. (bind(l-4).and.bind(l-3))) go to 90
        if(sql.le.sqmax) go to 90
        k = l
        sqmax = sql
        nr = i1+(i-i1)/2
  90    l = l+1
        i1 = i
        sql = 0.
 100    sql = sql+term
 110  continue
      if(m-i1.le.1 .or. (bind(l-4).and.bind(l-3))) go to 120
      if(sql.le.sqmax) go to 120
      k = l
      nr = i1+(m-i1)/2
c  if no such interval is found, control is repassed to the driver
c  program (ier = -1).
 120  if(nr.eq.0) go to 190
c  if s(x) is given by the same straight line in two succeeding knot
c  intervals t(l-1) <= x <= t(l) and t(l) <= x <= t(l+1),delete t(l)
      n8 = n-8
      l1 = 0
      if(n8.le.0) go to 150
      do 140 i=1,n8
        if(.not. (bind(i).and.bind(i+1).and.bind(i+2))) go to 140
        l = i+4-l1
        if(k.gt.l) k = k-1
        n = n-1
        l1 = l1+1
        do 130 j=l,n
          t(j) = t(j+1)
 130    continue
 140  continue
c  test whether we cannot further increase the number of knots.
 150  if(n.eq.nmax) go to 180
      if(n.eq.nest) go to 170
c  locate an additional knot at the point x(nr).
      j = n
      do 160 i=k,n
        t(j+1) = t(j)
        j = j-1
 160  continue
      t(k) = x(nr)
      n = n+1
c  restart the computations with the new set of knots.
      go to 40
c  error codes and messages.
 170  ier = -3
      go to 300
 180  ier = -2
      go to 300
 190  ier = -1
      go to 300
 200  ier = 4
      go to 300
 210  ier = 5
 300  return
      end
