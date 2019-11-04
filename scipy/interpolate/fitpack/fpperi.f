      subroutine fpperi(iopt,x,y,w,m,k,s,nest,tol,maxit,k1,k2,n,t,c,
     * fp,fpint,z,a1,a2,b,g1,g2,q,nrdata,ier)
c  ..
c  ..scalar arguments..
      real*8 s,tol,fp
      integer iopt,m,k,nest,maxit,k1,k2,n,ier
c  ..array arguments..
      real*8 x(m),y(m),w(m),t(nest),c(nest),fpint(nest),z(nest),
     * a1(nest,k1),a2(nest,k),b(nest,k2),g1(nest,k2),g2(nest,k1),
     * q(m,k1)
      integer nrdata(nest)
c  ..local scalars..
      real*8 acc,cos,c1,d1,fpart,fpms,fpold,fp0,f1,f2,f3,p,per,pinv,piv,
     *
     * p1,p2,p3,sin,store,term,wi,xi,yi,rn,one,con1,con4,con9,half
      integer i,ich1,ich3,ij,ik,it,iter,i1,i2,i3,j,jk,jper,j1,j2,kk,
     * kk1,k3,l,l0,l1,l5,mm,m1,new,nk1,nk2,nmax,nmin,nplus,npl1,
     * nrint,n10,n11,n7,n8
c  ..local arrays..
      real*8 h(6),h1(7),h2(6)
c  ..function references..
      real*8 abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares periodic spline   c
c  sinf(x). if the sum f(p=inf) <= s we accept the choice of knots.    c
c  the initial choice of knots depends on the value of s and iopt.     c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+2*k.                                        c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial of         c
c      degree k; n = nmin = 2*k+2. since s(x) must be periodic we      c
c      find that s(x) is a constant function.                          c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the least-squares periodic polynomial.      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      m1 = m-1
      kk = k
      kk1 = k1
      k3 = 3*k+1
      nmin = 2*k1
c  determine the length of the period of s(x).
      per = x(m)-x(1)
      if(iopt.lt.0) go to 50
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  determine nmax, the number of knots for periodic spline interpolation
      nmax = m+2*k
      if(s.gt.0. .or. nmax.eq.nmin) go to 30
c  if s=0, s(x) is an interpolating spline.
      n = nmax
c  test whether the required storage space exceeds the available one.
      if(n.gt.nest) go to 620
c  find the position of the interior knots in case of interpolation.
   5  if((k/2)*2 .eq. k) go to 20
      do 10 i=2,m1
        j = i+k
        t(j) = x(i)
  10  continue
      if(s.gt.0.) go to 50
      kk = k-1
      kk1 = k
      if(kk.gt.0) go to 50
      t(1) = t(m)-per
      t(2) = x(1)
      t(m+1) = x(m)
      t(m+2) = t(3)+per
      do 15 i=1,m1
        c(i) = y(i)
  15  continue
      c(m) = c(1)
      fp = 0.
      fpint(n) = fp0
      fpint(n-1) = 0.
      nrdata(n) = 0
      go to 630
  20  do 25 i=2,m1
        j = i+k
        t(j) = (x(i)+x(i-1))*half
  25  continue
      go to 50
c  if s > 0 our initial choice depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  periodic polynomial. (i.e. a constant function).
c  if iopt=1 and fp0>s we start computing the least-squares periodic
c  spline according the set of knots found at the last call of the
c  routine.
  30  if(iopt.eq.0) go to 35
      if(n.eq.nmin) go to 35
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 50
c  the case that s(x) is a constant function is treated separetely.
c  find the least-squares constant c1 and compute fp0 at the same time.
  35  fp0 = 0.
      d1 = 0.
      c1 = 0.
      do 40 it=1,m1
        wi = w(it)
        yi = y(it)*wi
        call fpgivs(wi,d1,cos,sin)
        call fprota(cos,sin,yi,c1)
        fp0 = fp0+yi**2
  40  continue
      c1 = c1/d1
c  test whether that constant function is a solution of our problem.
      fpms = fp0-s
      if(fpms.lt.acc .or. nmax.eq.nmin) go to 640
      fpold = fp0
c  test whether the required storage space exceeds the available one.
      if(nmin.ge.nest) go to 620
c  start computing the least-squares periodic spline with one
c  interior knot.
      nplus = 1
      n = nmin+1
      mm = (m+1)/2
      t(k2) = x(mm)
      nrdata(1) = mm-2
      nrdata(2) = m1-mm
c  main loop for the different sets of knots. m is a save upper
c  bound for the number of trials.
  50  do 340 iter=1,m
c  find nrint, the number of knot intervals.
        nrint = n-nmin+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(x). if we take
c      t(k+1) = x(1), t(n-k) = x(m)
c      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
c      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
c  then s(x) is a periodic spline with period per if the b-spline
c  coefficients satisfy the following conditions
c      c(n7+j) = c(j), j=1,...k   (**)   with n7=n-2*k-1.
        t(k1) = x(1)
        nk1 = n-k1
        nk2 = nk1+1
        t(nk2) = x(m)
        do 60 j=1,k
          i1 = nk2+j
          i2 = nk2-j
          j1 = k1+j
          j2 = k1-j
          t(i1) = t(j1)+per
          t(j2) = t(i2)-per
  60    continue
c  compute the b-spline coefficients c(j),j=1,...n7 of the least-squares
c  periodic spline sinf(x). the observation matrix a is built up row
c  by row while taking into account condition (**) and is reduced to
c  triangular form by givens transformations .
c  at the same time fp=f(p=inf) is computed.
c  the n7 x n7 triangularised upper matrix a has the form
c            ! a1 '    !
c        a = !    ' a2 !
c            ! 0  '    !
c  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
c  matrix of bandwidth k+1 ( n10 = n7-k).
c  initialization.
        do 70 i=1,nk1
          z(i) = 0.
          do 70 j=1,kk1
            a1(i,j) = 0.
  70    continue
        n7 = nk1-k
        n10 = n7-kk
        jper = 0
        fp = 0.
        l = k1
        do 290 it=1,m1
c  fetch the current data point x(it),y(it)
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
c  search for knot interval t(l) <= xi < t(l+1).
  80      if(xi.lt.t(l+1)) go to 85
          l = l+1
          go to 80
c  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  85      call fpbspl(t,n,k,xi,l,h)
          do 90 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  90      continue
          l5 = l-k1
c  test whether the b-splines nj,k+1(x),j=1+n7,...nk1 are all zero at xi
          if(l5.lt.n10) go to 285
          if(jper.ne.0) go to 160
c  initialize the matrix a2.
          do 95 i=1,n7
          do 95 j=1,kk
              a2(i,j) = 0.
  95      continue
          jk = n10+1
          do 110 i=1,kk
            ik = jk
            do 100 j=1,kk1
              if(ik.le.0) go to 105
              a2(ik,i) = a1(ik,j)
              ik = ik-1
 100        continue
 105        jk = jk+1
 110      continue
          jper = 1
c  if one of the b-splines nj,k+1(x),j=n7+1,...nk1 is not zero at xi
c  we take account of condition (**) for setting up the new row
c  of the observation matrix a. this row is stored in the arrays h1
c  (the part with respect to a1) and h2 (the part with
c  respect to a2).
 160      do 170 i=1,kk
            h1(i) = 0.
            h2(i) = 0.
 170      continue
          h1(kk1) = 0.
          j = l5-n10
          do 210 i=1,kk1
            j = j+1
            l0 = j
 180        l1 = l0-kk
            if(l1.le.0) go to 200
            if(l1.le.n10) go to 190
            l0 = l1-n10
            go to 180
 190        h1(l1) = h(i)
            go to 210
 200        h2(l0) = h2(l0)+h(i)
 210      continue
c  rotate the new row of the observation matrix into triangle
c  by givens transformations.
          if(n10.le.0) go to 250
c  rotation with the rows 1,2,...n10 of matrix a.
          do 240 j=1,n10
            piv = h1(1)
            if(piv.ne.0.) go to 214
            do 212 i=1,kk
              h1(i) = h1(i+1)
 212        continue
            h1(kk1) = 0.
            go to 240
c  calculate the parameters of the givens transformation.
 214        call fpgivs(piv,a1(j,1),cos,sin)
c  transformation to the right hand side.
            call fprota(cos,sin,yi,z(j))
c  transformations to the left hand side with respect to a2.
            do 220 i=1,kk
              call fprota(cos,sin,h2(i),a2(j,i))
 220        continue
            if(j.eq.n10) go to 250
            i2 = min0(n10-j,kk)
c  transformations to the left hand side with respect to a1.
            do 230 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),a1(j,i1))
              h1(i) = h1(i1)
 230        continue
            h1(i1) = 0.
 240      continue
c  rotation with the rows n10+1,...n7 of matrix a.
 250      do 270 j=1,kk
            ij = n10+j
            if(ij.le.0) go to 270
            piv = h2(j)
            if(piv.eq.0.) go to 270
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a2(ij,j),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,z(ij))
            if(j.eq.kk) go to 280
            j1 = j+1
c  transformations to left hand side.
            do 260 i=j1,kk
              call fprota(cos,sin,h2(i),a2(ij,i))
 260        continue
 270      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 280      fp = fp+yi**2
          go to 290
c  rotation of the new row of the observation matrix into
c  triangle in case the b-splines nj,k+1(x),j=n7+1,...n-k-1 are all zero
c  at xi.
 285      j = l5
          do 140 i=1,kk1
            j = j+1
            piv = h(i)
            if(piv.eq.0.) go to 140
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a1(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,z(j))
            if(i.eq.kk1) go to 150
            i2 = 1
            i3 = i+1
c  transformations to left hand side.
            do 130 i1=i3,kk1
              i2 = i2+1
              call fprota(cos,sin,h(i1),a1(j,i2))
 130        continue
 140      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 150      fp = fp+yi**2
 290    continue
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
c  backward substitution to obtain the b-spline coefficients c(j),j=1,.n
        call fpbacp(a1,a2,z,n7,kk,c,kk1,nest)
c  calculate from condition (**) the coefficients c(j+n7),j=1,2,...k.
        do 295 i=1,k
          j = i+n7
          c(j) = c(i)
 295    continue
        if(iopt.lt.0) go to 660
c  test whether the approximation sinf(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 660
c  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.) go to 350
c  if n=nmax, sinf(x) is an interpolating spline.
        if(n.eq.nmax) go to 630
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of the
c  storage capacity limitation.
        if(n.eq.nest) go to 620
c  determine the number of knots nplus we are going to add.
        npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
        fpold = fp
c  compute the sum(wi*(yi-s(xi))**2) for each knot interval
c  t(j+k) <= xi <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.
        i = 1
        l = k1
        do 320 it=1,m1
          if(x(it).lt.t(l)) go to 300
          new = 1
          l = l+1
 300      term = 0.
          l0 = l-k2
          do 310 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 310      continue
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new.eq.0) go to 320
          if(l.gt.k2) go to 315
          fpint(nrint) = term
          new = 0
          go to 320
 315      store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 320    continue
        fpint(nrint) = fpint(nrint)+fpart
        do 330 l=1,nplus
c  add a new knot
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
c  if n=nmax we locate the knots as for interpolation.
          if(n.eq.nmax) go to 5
c  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 340
 330    continue
c  restart the computations with the new set of knots.
 340  continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 2: determination of the smoothing periodic spline sp(x).       c
c  *************************************************************       c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing spline    c
c  sp(x). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(x) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/sqrt(p).      c
c  iteratively we then have to determine the value of p such that      c
c  f(p)=sum(w(i)*(y(i)-sp(x(i)))**2) be = s. we already know that      c
c  the least-squares constant function corresponds to p=0, and that    c
c  the least-squares periodic spline corresponds to p=infinity. the    c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
 350  call fpdisc(t,n,k2,b,nest)
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      n11 = n10-1
      n8 = n7-1
      p = 0.
      l = n7
      do 352 i=1,k
         j = k+1-i
         p = p+a2(l,j)
         l = l-1
         if(l.eq.0) go to 356
 352  continue
      do 354 i=1,n10
         p = p+a1(i,1)
 354  continue
 356  rn = n7
      p = rn/p
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p) = s.
      do 595 iter=1,maxit
c  form the matrix g  as the matrix a extended by the rows of matrix b.
c  the rows of matrix b with weight 1/p are rotated into
c  the triangularised observation matrix a.
c  after triangularisation our n7 x n7 matrix g takes the form
c            ! g1 '    !
c        g = !    ' g2 !
c            ! 0  '    !
c  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
c  matrix of bandwidth k+2. ( n11 = n7-k-1)
        pinv = one/p
c  store matrix a into g
        do 360 i=1,n7
          c(i) = z(i)
          g1(i,k1) = a1(i,k1)
          g1(i,k2) = 0.
          g2(i,1) = 0.
          do 360 j=1,k
            g1(i,j) = a1(i,j)
            g2(i,j+1) = a2(i,j)
 360    continue
        l = n10
        do 370 j=1,k1
          if(l.le.0) go to 375
          g2(l,1) = a1(l,j)
          l = l-1
 370    continue
 375    do 540 it=1,n8
c  fetch a new row of matrix b and store it in the arrays h1 (the part
c  with respect to g1) and h2 (the part with respect to g2).
          yi = 0.
          do 380 i=1,k1
            h1(i) = 0.
            h2(i) = 0.
 380      continue
          h1(k2) = 0.
          if(it.gt.n11) go to 420
          l = it
          l0 = it
          do 390 j=1,k2
            if(l0.eq.n10) go to 400
            h1(j) = b(it,j)*pinv
            l0 = l0+1
 390      continue
          go to 470
 400      l0 = 1
          do 410 l1=j,k2
            h2(l0) = b(it,l1)*pinv
            l0 = l0+1
 410      continue
          go to 470
 420      l = 1
          i = it-n10
          do 460 j=1,k2
            i = i+1
            l0 = i
 430        l1 = l0-k1
            if(l1.le.0) go to 450
            if(l1.le.n11) go to 440
            l0 = l1-n11
            go to 430
 440        h1(l1) = b(it,j)*pinv
            go to 460
 450        h2(l0) = h2(l0)+b(it,j)*pinv
 460      continue
          if(n11.le.0) go to 510
c  rotate this row into triangle by givens transformations without
c  square roots.
c  rotation with the rows l,l+1,...n11.
 470      do 500 j=l,n11
            piv = h1(1)
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,g1(j,1),cos,sin)
c  transformation to right hand side.
            call fprota(cos,sin,yi,c(j))
c  transformation to the left hand side with respect to g2.
            do 480 i=1,k1
              call fprota(cos,sin,h2(i),g2(j,i))
 480        continue
            if(j.eq.n11) go to 510
            i2 = min0(n11-j,k1)
c  transformation to the left hand side with respect to g1.
            do 490 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),g1(j,i1))
              h1(i) = h1(i1)
 490        continue
            h1(i1) = 0.
 500      continue
c  rotation with the rows n11+1,...n7
 510      do 530 j=1,k1
            ij = n11+j
            if(ij.le.0) go to 530
            piv = h2(j)
c  calculate the parameters of the givens transformation
            call fpgivs(piv,g2(ij,j),cos,sin)
c  transformation to the right hand side.
            call fprota(cos,sin,yi,c(ij))
            if(j.eq.k1) go to 540
            j1 = j+1
c  transformation to the left hand side.
            do 520 i=j1,k1
              call fprota(cos,sin,h2(i),g2(ij,i))
 520        continue
 530      continue
 540    continue
c  backward substitution to obtain the b-spline coefficients
c  c(j),j=1,2,...n7 of sp(x).
        call fpbacp(g1,g2,c,n7,k1,c,k2,nest)
c  calculate from condition (**) the b-spline coefficients c(n7+j),j=1,.
        do 545 i=1,k
          j = i+n7
          c(j) = c(i)
 545    continue
c  computation of f(p).
        fp = 0.
        l = k1
        do 570 it=1,m1
          if(x(it).lt.t(l)) go to 550
          l = l+1
 550      l0 = l-k2
          term = 0.
          do 560 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 560      continue
          fp = fp+(w(it)*(term-y(it)))**2
 570    continue
c  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 660
c  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 600
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 580
        if((f2-f3) .gt. acc) go to 575
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 +p2*con1
        go to 595
 575    if(f2.lt.0.) ich3 = 1
 580    if(ich1.ne.0) go to 590
        if((f1-f2) .gt. acc) go to 585
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 595
        if(p.ge.p3) p = p2*con1 +p3*con9
        go to 595
 585    if(f2.gt.0.) ich1 = 1
c  test whether the iteration process proceeds as theoretically
c  expected.
 590    if(f2.ge.f1 .or. f2.le.f3) go to 610
c  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 595  continue
c  error codes and messages.
 600  ier = 3
      go to 660
 610  ier = 2
      go to 660
 620  ier = 1
      go to 660
 630  ier = -1
      go to 660
 640  ier = -2
c  the least-squares constant function c1 is a solution of our problem.
c  a constant function is a spline of degree k with all b-spline
c  coefficients equal to that constant c1.
      do 650 i=1,k1
        rn = k1-i
        t(i) = x(1)-rn*per
        c(i) = c1
        j = i+k1
        rn = i-1
        t(j) = x(m)+rn*per
 650  continue
      n = nmin
      fp = fp0
      fpint(n) = fp0
      fpint(n-1) = 0.
      nrdata(n) = 0
 660  return
      end
