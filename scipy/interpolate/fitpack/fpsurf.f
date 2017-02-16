      subroutine fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kxx,kyy,s,nxest,
     * nyest,eta,tol,maxit,nmax,km1,km2,ib1,ib3,nc,intest,nrest,
     * nx0,tx,ny0,ty,c,fp,fp0,fpint,coord,f,ff,a,q,bx,by,spx,spy,h,
     * index,nummer,wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      real*8 xb,xe,yb,ye,s,eta,tol,fp,fp0
      integer iopt,m,kxx,kyy,nxest,nyest,maxit,nmax,km1,km2,ib1,ib3,
     * nc,intest,nrest,nx0,ny0,lwrk,ier
c  ..array arguments..
      real*8 x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),c(nc),fpint(intest),
     * coord(intest),f(nc),ff(nc),a(nc,ib1),q(nc,ib3),bx(nmax,km2),
     * by(nmax,km2),spx(m,km1),spy(m,km1),h(ib3),wrk(lwrk)
      integer index(nrest),nummer(m)
c  ..local scalars..
      real*8 acc,arg,cos,dmax,fac1,fac2,fpmax,fpms,f1,f2,f3,hxi,p,pinv,
     * piv,p1,p2,p3,sigma,sin,sq,store,wi,x0,x1,y0,y1,zi,eps,
     * rn,one,con1,con9,con4,half,ten
      integer i,iband,iband1,iband3,iband4,ibb,ichang,ich1,ich3,ii,
     * in,irot,iter,i1,i2,i3,j,jrot,jxy,j1,kx,kx1,kx2,ky,ky1,ky2,l,
     * la,lf,lh,lwest,lx,ly,l1,l2,n,ncof,nk1x,nk1y,nminx,nminy,nreg,
     * nrint,num,num1,nx,nxe,nxx,ny,nye,nyy,n1,rank
c  ..local arrays..
      real*8 hx(6),hy(6)
c  ..function references..
      real*8 abs,fprati,sqrt
      integer min0
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
      ten = 0.1e+02
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c given a set of knots we compute the least-squares spline sinf(x,y),  c
c and the corresponding weighted sum of squared residuals fp=f(p=inf). c
c if iopt=-1  sinf(x,y) is the requested approximation.                c
c if iopt=0 or iopt=1 we check whether we can accept the knots:        c
c   if fp <=s we will continue with the current set of knots.          c
c   if fp > s we will increase the number of knots and compute the     c
c      corresponding least-squares spline until finally  fp<=s.        c
c the initial choice of knots depends on the value of s and iopt.      c
c   if iopt=0 we first compute the least-squares polynomial of degree  c
c     kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        c
c     fp0=f(0) denotes the corresponding weighted sum of squared       c
c     residuals                                                        c
c   if iopt=1 we start with the knots found at the last call of the    c
c     routine, except for the case that s>=fp0; then we can compute    c
c     the least-squares polynomial directly.                           c
c eventually the independent variables x and y (and the corresponding  c
c parameters) will be switched if this can reduce the bandwidth of the c
c system to be solved.                                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ichang denotes whether(1) or not(-1) the directions have been inter-
c  changed.
      ichang = -1
      x0 = xb
      x1 = xe
      y0 = yb
      y1 = ye
      kx = kxx
      ky = kyy
      kx1 = kx+1
      ky1 = ky+1
      nxe = nxest
      nye = nyest
      eps = sqrt(eta)
      if(iopt.lt.0) go to 20
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      if(iopt.eq.0) go to 10
      if(fp0.gt.s) go to 20
c  initialization for the least-squares polynomial.
  10  nminx = 2*kx1
      nminy = 2*ky1
      nx = nminx
      ny = nminy
      ier = -2
      go to 30
  20  nx = nx0
      ny = ny0
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  30  do 420 iter=1,m
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(x,y).
        l = nx
        do 40 i=1,kx1
          tx(i) = x0
          tx(l) = x1
          l = l-1
  40    continue
        l = ny
        do 50 i=1,ky1
          ty(i) = y0
          ty(l) = y1
          l = l-1
  50    continue
c  find nrint, the total number of knot intervals and nreg, the number
c  of panels in which the approximation domain is subdivided by the
c  intersection of knots.
        nxx = nx-2*kx1+1
        nyy = ny-2*ky1+1
        nrint = nxx+nyy
        nreg = nxx*nyy
c  find the bandwidth of the observation matrix a.
c  if necessary, interchange the variables x and y, in order to obtain
c  a minimal bandwidth.
        iband1 = kx*(ny-ky1)+ky
        l = ky*(nx-kx1)+kx
        if(iband1.le.l) go to 130
        iband1 = l
        ichang = -ichang
        do 60 i=1,m
          store = x(i)
          x(i) = y(i)
          y(i) = store
  60    continue
        store = x0
        x0 = y0
        y0 = store
        store = x1
        x1 = y1
        y1 = store
        n = min0(nx,ny)
        do 70 i=1,n
          store = tx(i)
          tx(i) = ty(i)
          ty(i) = store
  70    continue
        n1 = n+1
        if (nx.lt.ny) go to 80
        if (nx.eq.ny) go to 120
        go to 100
  80    do 90 i=n1,ny
          tx(i) = ty(i)
  90    continue
        go to 120
 100    do 110 i=n1,nx
          ty(i) = tx(i)
 110    continue
 120    l = nx
        nx = ny
        ny = l
        l = nxe
        nxe = nye
        nye = l
        l = nxx
        nxx = nyy
        nyy = l
        l = kx
        kx = ky
        ky = l
        kx1 = kx+1
        ky1 = ky+1
 130    iband = iband1+1
c  arrange the data points according to the panel they belong to.
        call fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
c  find ncof, the number of b-spline coefficients.
        nk1x = nx-kx1
        nk1y = ny-ky1
        ncof = nk1x*nk1y
c  initialize the observation matrix a.
        do 140 i=1,ncof
          f(i) = 0.
          do 140 j=1,iband
            a(i,j) = 0.
 140    continue
c  initialize the sum of squared residuals.
        fp = 0.
c  fetch the data points in the new order. main loop for the
c  different panels.
        do 250 num=1,nreg
c  fix certain constants for the current panel; jrot records the column
c  number of the first non-zero element in a row of the observation
c  matrix according to a data point of the panel.
          num1 = num-1
          lx = num1/nyy
          l1 = lx+kx1
          ly = num1-lx*nyy
          l2 = ly+ky1
          jrot = lx*nk1y+ly
c  test whether there are still data points in the panel.
          in = index(num)
 150      if(in.eq.0) go to 250
c  fetch a new data point.
          wi = w(in)
          zi = z(in)*wi
c  evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in).
          call fpbspl(tx,nx,kx,x(in),l1,hx)
c  evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in).
          call fpbspl(ty,ny,ky,y(in),l2,hy)
c  store the value of these b-splines in spx and spy respectively.
          do 160 i=1,kx1
            spx(in,i) = hx(i)
 160      continue
          do 170 i=1,ky1
            spy(in,i) = hy(i)
 170      continue
c  initialize the new row of observation matrix.
          do 180 i=1,iband
            h(i) = 0.
 180      continue
c  calculate the non-zero elements of the new row by making the cross
c  products of the non-zero b-splines in x- and y-direction.
          i1 = 0
          do 200 i=1,kx1
            hxi = hx(i)
            j1 = i1
            do 190 j=1,ky1
              j1 = j1+1
              h(j1) = hxi*hy(j)*wi
 190        continue
            i1 = i1+nk1y
 200      continue
c  rotate the row into triangle by givens transformations .
          irot = jrot
          do 220 i=1,iband
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 220
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),cos,sin)
c  apply that transformation to the right hand side.
            call fprota(cos,sin,zi,f(irot))
            if(i.eq.iband) go to 230
c  apply that transformation to the left hand side.
            i2 = 1
            i3 = i+1
            do 210 j=i3,iband
              i2 = i2+1
              call fprota(cos,sin,h(j),a(irot,i2))
 210        continue
 220      continue
c  add the contribution of the row to the sum of squares of residual
c  right hand sides.
 230      fp = fp+zi**2
c  find the number of the next data point in the panel.
          in = nummer(in)
          go to 150
 250    continue
c  find dmax, the maximum value for the diagonal elements in the reduced
c  triangle.
        dmax = 0.
        do 260 i=1,ncof
          if(a(i,1).le.dmax) go to 260
          dmax = a(i,1)
 260    continue
c  check whether the observation matrix is rank deficient.
        sigma = eps*dmax
        do 270 i=1,ncof
          if(a(i,1).le.sigma) go to 280
 270    continue
c  backward substitution in case of full rank.
        call fpback(a,f,ncof,iband,c,nc)
        rank = ncof
        do 275 i=1,ncof
          q(i,1) = a(i,1)/dmax
 275    continue
        go to 300
c  in case of rank deficiency, find the minimum norm solution.
c  check whether there is sufficient working space
 280    lwest = ncof*iband+ncof+iband
        if(lwrk.lt.lwest) go to 780
        do 290 i=1,ncof
          ff(i) = f(i)
          do 290 j=1,iband
            q(i,j) = a(i,j)
 290    continue
        lf =1
        lh = lf+ncof
        la = lh+iband
        call fprank(q,ff,ncof,iband,nc,sigma,c,sq,rank,wrk(la),
     *    wrk(lf),wrk(lh))
        do 295 i=1,ncof
          q(i,1) = q(i,1)/dmax
 295    continue
c  add to the sum of squared residuals, the contribution of reducing
c  the rank.
        fp = fp+sq
 300    if(ier.eq.(-2)) fp0 = fp
c  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) go to 820
        fpms = fp-s
        if(abs(fpms).le.acc) then
          if (fp.le.0) go to 815
          go to 820
        endif
c  test whether we can accept the choice of knots.
        if(fpms.lt.0.) go to 430
c  test whether we cannot further increase the number of knots.
        if(ncof.gt.m) go to 790
        ier = 0
c  search where to add a new knot.
c  find for each interval the sum of squared residuals fpint for the
c  data points having the coordinate belonging to that knot interval.
c  calculate also coord which is the same sum, weighted by the position
c  of the data points considered.
        do 320 i=1,nrint
          fpint(i) = 0.
          coord(i) = 0.
 320    continue
        do 360 num=1,nreg
          num1 = num-1
          lx = num1/nyy
          l1 = lx+1
          ly = num1-lx*nyy
          l2 = ly+1+nxx
          jrot = lx*nk1y+ly
          in = index(num)
 330      if(in.eq.0) go to 360
          store = 0.
          i1 = jrot
          do 350 i=1,kx1
            hxi = spx(in,i)
            j1 = i1
            do 340 j=1,ky1
              j1 = j1+1
              store = store+hxi*spy(in,j)*c(j1)
 340        continue
            i1 = i1+nk1y
 350      continue
          store = (w(in)*(z(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*x(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*y(in)
          in = nummer(in)
          go to 330
 360    continue
c  find the interval for which fpint is maximal on the condition that
c  there still can be added a knot.
 370    l = 0
        fpmax = 0.
        l1 = 1
        l2 = nrint
        if(nx.eq.nxe) l1 = nxx+1
        if(ny.eq.nye) l2 = nxx
        if(l1.gt.l2) go to 810
        do 380 i=l1,l2
          if(fpmax.ge.fpint(i)) go to 380
          l = i
          fpmax = fpint(i)
 380    continue
c  test whether we cannot further increase the number of knots.
        if(l.eq.0) go to 785
c  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
c  test in what direction the new knot is going to be added.
        if(l.gt.nxx) go to 400
c  addition in the x-direction.
        jxy = l+kx1
        fpint(l) = 0.
        fac1 = tx(jxy)-arg
        fac2 = arg-tx(jxy-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 370
        j = nx
        do 390 i=jxy,nx
          tx(j+1) = tx(j)
          j = j-1
 390    continue
        tx(jxy) = arg
        nx = nx+1
        go to 420
c  addition in the y-direction.
 400    jxy = l+ky1-nxx
        fpint(l) = 0.
        fac1 = ty(jxy)-arg
        fac2 = arg-ty(jxy-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 370
        j = ny
        do 410 i=jxy,ny
          ty(j+1) = ty(j)
          j = j-1
 410    continue
        ty(jxy) = arg
        ny = ny+1
c  restart the computations with the new set of knots.
 420  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 430  if(ier.eq.(-2)) go to 830
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(x,y)                c
c *****************************************************                c
c we have determined the number of knots and their position. we now    c
c compute the b-spline coefficients of the smoothing spline sp(x,y).   c
c the observation matrix a is extended by the rows of a matrix,        c
c expressing that sp(x,y) must be a polynomial of degree kx in x and   c
c ky in y. the corresponding weights of these additional rows are set  c
c to 1./p.  iteratively we than have to determine the value of p       c
c such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.           c
c we already know that the least-squares polynomial corresponds to     c
c p=0  and that the least-squares spline corresponds to p=infinity.    c
c the iteration process which is proposed here makes use of rational   c
c interpolation. since f(p) is a convex and strictly decreasing        c
c function of p, it can be approximated by a rational function r(p)=   c
c (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values c
c of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the c
c new value of p such that r(p)=s. convergence is guaranteed by taking c
c f1 > 0 and f3 < 0.                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      kx2 = kx1+1
c  test whether there are interior knots in the x-direction.
      if(nk1x.eq.kx1) go to 440
c  evaluate the discotinuity jumps of the kx-th order derivative of
c  the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1.
      call fpdisc(tx,nx,kx2,bx,nmax)
 440  ky2 = ky1 + 1
c  test whether there are interior knots in the y-direction.
      if(nk1y.eq.ky1) go to 450
c  evaluate the discontinuity jumps of the ky-th order derivative of
c  the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1.
      call fpdisc(ty,ny,ky2,by,nmax)
c  initial value for p.
 450  p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 460 i=1,ncof
        p = p+a(i,1)
 460  continue
      rn = ncof
      p = rn/p
c  find the bandwidth of the extended observation matrix.
      iband3 = kx1*nk1y
      iband4 = iband3 +1
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 770 iter=1,maxit
        pinv = one/p
c  store the triangularized observation matrix into q.
        do 480 i=1,ncof
          ff(i) = f(i)
          do 470 j=1,iband
            q(i,j) = a(i,j)
 470      continue
          ibb = iband+1
          do 480 j=ibb,iband4
            q(i,j) = 0.
 480    continue
        if(nk1y.eq.ky1) go to 560
c  extend the observation matrix with the rows of a matrix, expressing
c  that for x=cst. sp(x,y) must be a polynomial in y of degree ky.
        do 550 i=ky2,nk1y
          ii = i-ky1
          do 550 j=1,nk1x
c  initialize the new row.
            do 490 l=1,iband
              h(l) = 0.
 490        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            do 500 l=1,ky2
              h(l) = by(ii,l)*pinv
 500        continue
            zi = 0.
            jrot = (j-1)*nk1y+ii
c  rotate the new row into triangle by givens transformations without
c  square roots.
            do 540 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if(piv.eq.0.) then
                if (i2.le.0) go to 550
                go to 520
              endif
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),cos,sin)
c  apply that givens transformation to the right hand side.
              call fprota(cos,sin,zi,ff(irot))
              if(i2.eq.0) go to 550
c  apply that givens transformation to the left hand side.
              do 510 l=1,i2
                l1 = l+1
                call fprota(cos,sin,h(l1),q(irot,l1))
 510          continue
 520          do 530 l=1,i2
                h(l) = h(l+1)
 530          continue
              h(i2+1) = 0.
 540        continue
 550    continue
 560    if(nk1x.eq.kx1) go to 640
c  extend the observation matrix with the rows of a matrix expressing
c  that for y=cst. sp(x,y) must be a polynomial in x of degree kx.
        do 630 i=kx2,nk1x
          ii = i-kx1
          do 630 j=1,nk1y
c  initialize the new row
            do 570 l=1,iband4
              h(l) = 0.
 570        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            j1 = 1
            do 580 l=1,kx2
              h(j1) = bx(ii,l)*pinv
              j1 = j1+nk1y
 580        continue
            zi = 0.
            jrot = (i-kx2)*nk1y+j
c  rotate the new row into triangle by givens transformations .
            do 620 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if(piv.eq.0.) then
                if (i2.le.0) go to 630
                go to 600
              endif
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),cos,sin)
c  apply that givens transformation to the right hand side.
              call fprota(cos,sin,zi,ff(irot))
              if(i2.eq.0) go to 630
c  apply that givens transformation to the left hand side.
              do 590 l=1,i2
                l1 = l+1
                call fprota(cos,sin,h(l1),q(irot,l1))
 590          continue
 600          do 610 l=1,i2
                h(l) = h(l+1)
 610          continue
              h(i2+1) = 0.
 620        continue
 630    continue
c  find dmax, the maximum value for the diagonal elements in the
c  reduced triangle.
 640    dmax = 0.
        do 650 i=1,ncof
          if(q(i,1).le.dmax) go to 650
          dmax = q(i,1)
 650    continue
c  check whether the matrix is rank deficient.
        sigma = eps*dmax
        do 660 i=1,ncof
          if(q(i,1).le.sigma) go to 670
 660    continue
c  backward substitution in case of full rank.
        call fpback(q,ff,ncof,iband4,c,nc)
        rank = ncof
        go to 675
c  in case of rank deficiency, find the minimum norm solution.
 670    lwest = ncof*iband4+ncof+iband4
        if(lwrk.lt.lwest) go to 780
        lf = 1
        lh = lf+ncof
        la = lh+iband4
        call fprank(q,ff,ncof,iband4,nc,sigma,c,sq,rank,wrk(la),
     *   wrk(lf),wrk(lh))
 675    do 680 i=1,ncof
          q(i,1) = q(i,1)/dmax
 680    continue
c  compute f(p).
        fp = 0.
        do 720 num = 1,nreg
          num1 = num-1
          lx = num1/nyy
          ly = num1-lx*nyy
          jrot = lx*nk1y+ly
          in = index(num)
 690      if(in.eq.0) go to 720
          store = 0.
          i1 = jrot
          do 710 i=1,kx1
            hxi = spx(in,i)
            j1 = i1
            do 700 j=1,ky1
              j1 = j1+1
              store = store+hxi*spy(in,j)*c(j1)
 700        continue
            i1 = i1+nk1y
 710      continue
          fp = fp+(w(in)*(z(in)-store))**2
          in = nummer(in)
          go to 690
 720    continue
c  test whether the approximation sp(x,y) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).le.acc) go to 820
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 795
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 740
        if((f2-f3).gt.acc) go to 730
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 770
 730    if(f2.lt.0.) ich3 = 1
 740    if(ich1.ne.0) go to 760
        if((f1-f2).gt.acc) go to 750
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 770
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 770
 750    if(f2.gt.0.) ich1 = 1
c  test whether the iteration process proceeds as theoretically
c  expected.
 760    if(f2.ge.f1 .or. f2.le.f3) go to 800
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 770  continue
c  error codes and messages.
 780  ier = lwest
      go to 830
 785  ier = 5
      go to 830
 790  ier = 4
      go to 830
 795  ier = 3
      go to 830
 800  ier = 2
      go to 830
 810  ier = 1
      go to 830
 815  ier = -1
      fp = 0.
 820  if(ncof.ne.rank) ier = -rank
c  test whether x and y are in the original order.
 830  if(ichang.lt.0) go to 930
c  if not, interchange x and y once more.
      l1 = 1
      do 840 i=1,nk1x
        l2 = i
        do 840 j=1,nk1y
          f(l2) = c(l1)
          l1 = l1+1
          l2 = l2+nk1x
 840  continue
      do 850 i=1,ncof
        c(i) = f(i)
 850  continue
      do 860 i=1,m
        store = x(i)
        x(i) = y(i)
        y(i) = store
 860  continue
      n = min0(nx,ny)
      do 870 i=1,n
        store = tx(i)
        tx(i) = ty(i)
        ty(i) = store
 870  continue
      n1 = n+1
      if (nx.lt.ny) go to 880
      if (nx.eq.ny) go to 920
      go to 900
 880  do 890 i=n1,ny
        tx(i) = ty(i)
 890  continue
      go to 920
 900  do 910 i=n1,nx
        ty(i) = tx(i)
 910  continue
 920  l = nx
      nx = ny
      ny = l
 930  if(iopt.lt.0) go to 940
      nx0 = nx
      ny0 = ny
 940  return
      end

