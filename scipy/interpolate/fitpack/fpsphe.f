      subroutine fpsphe(iopt,m,teta,phi,r,w,s,ntest,npest,eta,tol,maxit,
     *
     * ib1,ib3,nc,ncc,intest,nrest,nt,tt,np,tp,c,fp,sup,fpint,coord,f,
     * ff,row,coco,cosi,a,q,bt,bp,spt,spp,h,index,nummer,wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      integer iopt,m,ntest,npest,maxit,ib1,ib3,nc,ncc,intest,nrest,
     * nt,np,lwrk,ier
      real*8 s,eta,tol,fp,sup
c  ..array arguments..
      real*8 teta(m),phi(m),r(m),w(m),tt(ntest),tp(npest),c(nc),
     * fpint(intest),coord(intest),f(ncc),ff(nc),row(npest),coco(npest),
     *
     * cosi(npest),a(ncc,ib1),q(ncc,ib3),bt(ntest,5),bp(npest,5),
     * spt(m,4),spp(m,4),h(ib3),wrk(lwrk)
      integer index(nrest),nummer(m)
c  ..local scalars..
      real*8 aa,acc,arg,cn,co,c1,dmax,d1,d2,eps,facc,facs,fac1,fac2,fn,
     * fpmax,fpms,f1,f2,f3,hti,htj,p,pi,pinv,piv,pi2,p1,p2,p3,ri,si,
     * sigma,sq,store,wi,rn,one,con1,con9,con4,half,ten
      integer i,iband,iband1,iband3,iband4,ich1,ich3,ii,ij,il,in,irot,
     * iter,i1,i2,i3,j,jlt,jrot,j1,j2,l,la,lf,lh,ll,lp,lt,lwest,l1,l2,
     * l3,l4,ncof,ncoff,npp,np4,nreg,nrint,nrr,nr1,ntt,nt4,nt6,num,
     * num1,rank
c  ..local arrays..
      real*8 ht(4),hp(4)
c  ..function references..
      real*8 abs,atan,fprati,sqrt,cos,sin
      integer min0
c  ..subroutine references..
c   fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota,fprpsp
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
      ten = 0.1e+02
      pi = atan(one)*4
      pi2 = pi+pi
      eps = sqrt(eta)
      if(iopt.lt.0) go to 70
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      if(iopt.eq.0) go to 10
      if(s.lt.sup) then
        if (np.lt.11) go to 60
        go to 70
      endif
c  if iopt=0 we begin by computing the weighted least-squares polynomial
c  of the form
c     s(teta,phi) = c1*f1(teta) + cn*fn(teta)
c  where f1(teta) and fn(teta) are the cubic polynomials satisfying
c     f1(0) = 1, f1(pi) = f1'(0) = f1'(pi) = 0 ; fn(teta) = 1-f1(teta).
c  the corresponding weighted sum of squared residuals gives the upper
c  bound sup for the smoothing factor s.
  10  sup = 0.
      d1 = 0.
      d2 = 0.
      c1 = 0.
      cn = 0.
      fac1 = pi*(one + half)
      fac2 = (one + one)/pi**3
      aa = 0.
      do 40 i=1,m
         wi = w(i)
         ri = r(i)*wi
         arg = teta(i)
         fn = fac2*arg*arg*(fac1-arg)
         f1 = (one-fn)*wi
         fn = fn*wi
         if(fn.eq.0.) go to 20
         call fpgivs(fn,d1,co,si)
         call fprota(co,si,f1,aa)
         call fprota(co,si,ri,cn)
 20      if(f1.eq.0.) go to 30
         call fpgivs(f1,d2,co,si)
         call fprota(co,si,ri,c1)
 30      sup = sup+ri*ri
 40   continue
      if(d2.ne.0.) c1 = c1/d2
      if(d1.ne.0.) cn = (cn-aa*c1)/d1
c  find the b-spline representation of this least-squares polynomial
      nt = 8
      np = 8
      do 50 i=1,4
         c(i) = c1
         c(i+4) = c1
         c(i+8) = cn
         c(i+12) = cn
         tt(i) = 0.
         tt(i+4) = pi
         tp(i) = 0.
         tp(i+4) = pi2
  50  continue
      fp = sup
c  test whether the least-squares polynomial is an acceptable solution
      fpms = sup-s
      if(fpms.lt.acc) go to 960
c  test whether we cannot further increase the number of knots.
  60  if(npest.lt.11 .or. ntest.lt.9) go to 950
c  find the initial set of interior knots of the spherical spline in
c  case iopt = 0.
      np = 11
      tp(5) = pi*half
      tp(6) = pi
      tp(7) = tp(5)+pi
      nt = 9
      tt(5) = tp(5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1 : computation of least-squares spherical splines.            c
c  ********************************************************            c
c  if iopt < 0 we compute the least-squares spherical spline according c
c  to the given set of knots.                                          c
c  if iopt >=0 we compute least-squares spherical splines with increas-c
c  ing numbers of knots until the corresponding sum f(p=inf)<=s.       c
c  the initial set of knots then depends on the value of iopt:         c
c    if iopt=0 we start with one interior knot in the teta-direction   c
c              (pi/2) and three in the phi-direction (pi/2,pi,3*pi/2). c
c    if iopt>0 we start with the set of knots found at the last call   c
c              of the routine.                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  70  do 570 iter=1,m
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(teta,phi).
         l1 = 4
         l2 = l1
         l3 = np-3
         l4 = l3
         tp(l2) = 0.
         tp(l3) = pi2
         do 80 i=1,3
            l1 = l1+1
            l2 = l2-1
            l3 = l3+1
            l4 = l4-1
            tp(l2) = tp(l4)-pi2
            tp(l3) = tp(l1)+pi2
  80     continue
        l = nt
        do 90 i=1,4
          tt(i) = 0.
          tt(l) = pi
          l = l-1
  90    continue
c  find nrint, the total number of knot intervals and nreg, the number
c  of panels in which the approximation domain is subdivided by the
c  intersection of knots.
        ntt = nt-7
        npp = np-7
        nrr = npp/2
        nr1 = nrr+1
        nrint = ntt+npp
        nreg = ntt*npp
c  arrange the data points according to the panel they belong to.
        call fporde(teta,phi,m,3,3,tt,nt,tp,np,nummer,index,nreg)
c  find the b-spline coefficients coco and cosi of the cubic spline
c  approximations sc(phi) and ss(phi) for cos(phi) and sin(phi).
        do 100 i=1,npp
           coco(i) = 0.
           cosi(i) = 0.
           do 100 j=1,npp
              a(i,j) = 0.
 100    continue
c  the coefficients coco and cosi are obtained from the conditions
c  sc(tp(i))=cos(tp(i)),resp. ss(tp(i))=sin(tp(i)),i=4,5,...np-4.
        do 150 i=1,npp
           l2 = i+3
           arg = tp(l2)
           call fpbspl(tp,np,3,arg,l2,hp)
           do 110 j=1,npp
              row(j) = 0.
 110       continue
           ll = i
           do 120 j=1,3
              if(ll.gt.npp) ll= 1
              row(ll) = row(ll)+hp(j)
              ll = ll+1
 120       continue
           facc = cos(arg)
           facs = sin(arg)
           do 140 j=1,npp
              piv = row(j)
              if(piv.eq.0.) go to 140
              call fpgivs(piv,a(j,1),co,si)
              call fprota(co,si,facc,coco(j))
              call fprota(co,si,facs,cosi(j))
              if(j.eq.npp) go to 150
              j1 = j+1
              i2 = 1
              do 130 l=j1,npp
                 i2 = i2+1
                 call fprota(co,si,row(l),a(j,i2))
 130          continue
 140       continue
 150    continue
        call fpback(a,coco,npp,npp,coco,ncc)
        call fpback(a,cosi,npp,npp,cosi,ncc)
c  find ncof, the dimension of the spherical spline and ncoff, the
c  number of coefficients in the standard b-spline representation.
        nt4 = nt-4
        np4 = np-4
        ncoff = nt4*np4
        ncof = 6+npp*(ntt-1)
c  find the bandwidth of the observation matrix a.
        iband = 4*npp
        if(ntt.eq.4) iband = 3*(npp+1)
        if(ntt.lt.4) iband = ncof
        iband1 = iband-1
c  initialize the observation matrix a.
        do 160 i=1,ncof
          f(i) = 0.
          do 160 j=1,iband
            a(i,j) = 0.
 160    continue
c  initialize the sum of squared residuals.
        fp = 0.
c  fetch the data points in the new order. main loop for the
c  different panels.
        do 340 num=1,nreg
c  fix certain constants for the current panel; jrot records the column
c  number of the first non-zero element in a row of the observation
c  matrix according to a data point of the panel.
          num1 = num-1
          lt = num1/npp
          l1 = lt+4
          lp = num1-lt*npp+1
          l2 = lp+3
          lt = lt+1
          jrot = 0
          if(lt.gt.2) jrot = 3+(lt-3)*npp
c  test whether there are still data points in the current panel.
          in = index(num)
 170      if(in.eq.0) go to 340
c  fetch a new data point.
          wi = w(in)
          ri = r(in)*wi
c  evaluate for the teta-direction, the 4 non-zero b-splines at teta(in)
          call fpbspl(tt,nt,3,teta(in),l1,ht)
c  evaluate for the phi-direction, the 4 non-zero b-splines at phi(in)
          call fpbspl(tp,np,3,phi(in),l2,hp)
c  store the value of these b-splines in spt and spp resp.
          do 180 i=1,4
            spp(in,i) = hp(i)
            spt(in,i) = ht(i)
 180      continue
c  initialize the new row of observation matrix.
          do 190 i=1,iband
            h(i) = 0.
 190      continue
c  calculate the non-zero elements of the new row by making the cross
c  products of the non-zero b-splines in teta- and phi-direction and
c  by taking into account the conditions of the spherical splines.
          do 200 i=1,npp
             row(i) = 0.
 200      continue
c  take into account the condition (3) of the spherical splines.
          ll = lp
          do 210 i=1,4
             if(ll.gt.npp) ll=1
             row(ll) = row(ll)+hp(i)
             ll = ll+1
 210      continue
c  take into account the other conditions of the spherical splines.
          if(lt.gt.2 .and. lt.lt.(ntt-1)) go to 230
          facc = 0.
          facs = 0.
          do 220 i=1,npp
             facc = facc+row(i)*coco(i)
             facs = facs+row(i)*cosi(i)
 220     continue
c  fill in the non-zero elements of the new row.
 230     j1 = 0
         do 280 j =1,4
            jlt = j+lt
            htj = ht(j)
            if(jlt.gt.2 .and. jlt.le.nt4) go to 240
            j1 = j1+1
            h(j1) = h(j1)+htj
            go to 280
 240        if(jlt.eq.3 .or. jlt.eq.nt4) go to 260
            do 250 i=1,npp
               j1 = j1+1
               h(j1) = row(i)*htj
 250        continue
            go to 280
 260        if(jlt.eq.3) go to 270
            h(j1+1) = facc*htj
            h(j1+2) = facs*htj
            h(j1+3) = htj
            j1 = j1+2
            go to 280
 270        h(1) = h(1)+htj
            h(2) = facc*htj
            h(3) = facs*htj
            j1 = 3
 280      continue
          do 290 i=1,iband
            h(i) = h(i)*wi
 290      continue
c  rotate the row into triangle by givens transformations.
          irot = jrot
          do 310 i=1,iband
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 310
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
c  apply that transformation to the right hand side.
            call fprota(co,si,ri,f(irot))
            if(i.eq.iband) go to 320
c  apply that transformation to the left hand side.
            i2 = 1
            i3 = i+1
            do 300 j=i3,iband
              i2 = i2+1
              call fprota(co,si,h(j),a(irot,i2))
 300        continue
 310      continue
c  add the contribution of the row to the sum of squares of residual
c  right hand sides.
 320      fp = fp+ri**2
c  find the number of the next data point in the panel.
          in = nummer(in)
          go to 170
 340    continue
c  find dmax, the maximum value for the diagonal elements in the reduced
c  triangle.
        dmax = 0.
        do 350 i=1,ncof
          if(a(i,1).le.dmax) go to 350
          dmax = a(i,1)
 350    continue
c  check whether the observation matrix is rank deficient.
        sigma = eps*dmax
        do 360 i=1,ncof
          if(a(i,1).le.sigma) go to 370
 360    continue
c  backward substitution in case of full rank.
        call fpback(a,f,ncof,iband,c,ncc)
        rank = ncof
        do 365 i=1,ncof
          q(i,1) = a(i,1)/dmax
 365    continue
        go to 390
c  in case of rank deficiency, find the minimum norm solution.
 370    lwest = ncof*iband+ncof+iband
        if(lwrk.lt.lwest) go to 925
        lf = 1
        lh = lf+ncof
        la = lh+iband
        do 380 i=1,ncof
          ff(i) = f(i)
          do 380 j=1,iband
            q(i,j) = a(i,j)
 380    continue
        call fprank(q,ff,ncof,iband,ncc,sigma,c,sq,rank,wrk(la),
     *   wrk(lf),wrk(lh))
        do 385 i=1,ncof
          q(i,1) = q(i,1)/dmax
 385    continue
c  add to the sum of squared residuals, the contribution of reducing
c  the rank.
        fp = fp+sq
c  find the coefficients in the standard b-spline representation of
c  the spherical spline.
 390    call fprpsp(nt,np,coco,cosi,c,ff,ncoff)
c  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) then
          if (fp.le.0) go to 970
          go to 980
        endif
        fpms = fp-s
        if(abs(fpms).le.acc) then
          if (fp.le.0) go to 970
          go to 980
        endif
c  if f(p=inf) < s, accept the choice of knots.
        if(fpms.lt.0.) go to 580
c  test whether we cannot further increase the number of knots.
        if(ncof.gt.m) go to 935
c  search where to add a new knot.
c  find for each interval the sum of squared residuals fpint for the
c  data points having the coordinate belonging to that knot interval.
c  calculate also coord which is the same sum, weighted by the position
c  of the data points considered.
        do 450 i=1,nrint
          fpint(i) = 0.
          coord(i) = 0.
 450    continue
        do 490 num=1,nreg
          num1 = num-1
          lt = num1/npp
          l1 = lt+1
          lp = num1-lt*npp
          l2 = lp+1+ntt
          jrot = lt*np4+lp
          in = index(num)
 460      if(in.eq.0) go to 490
          store = 0.
          i1 = jrot
          do 480 i=1,4
            hti = spt(in,i)
            j1 = i1
            do 470 j=1,4
              j1 = j1+1
              store = store+hti*spp(in,j)*c(j1)
 470        continue
            i1 = i1+np4
 480      continue
          store = (w(in)*(r(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*teta(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*phi(in)
          in = nummer(in)
          go to 460
 490    continue
c  find the interval for which fpint is maximal on the condition that
c  there still can be added a knot.
        l1 = 1
        l2 = nrint
        if(ntest.lt.nt+1) l1=ntt+1
        if(npest.lt.np+2) l2=ntt
c  test whether we cannot further increase the number of knots.
        if(l1.gt.l2) go to 950
 500    fpmax = 0.
        l = 0
        do 510 i=l1,l2
          if(fpmax.ge.fpint(i)) go to 510
          l = i
          fpmax = fpint(i)
 510    continue
        if(l.eq.0) go to 930
c  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
c  test in what direction the new knot is going to be added.
        if(l.gt.ntt) go to 530
c  addition in the teta-direction
        l4 = l+4
        fpint(l) = 0.
        fac1 = tt(l4)-arg
        fac2 = arg-tt(l4-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500
        j = nt
        do 520 i=l4,nt
          tt(j+1) = tt(j)
          j = j-1
 520    continue
        tt(l4) = arg
        nt = nt+1
        go to 570
c  addition in the phi-direction
 530    l4 = l+4-ntt
        if(arg.lt.pi) go to 540
        arg = arg-pi
        l4 = l4-nrr
 540    fpint(l) = 0.
        fac1 = tp(l4)-arg
        fac2 = arg-tp(l4-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500
        ll = nrr+4
        j = ll
        do 550 i=l4,ll
          tp(j+1) = tp(j)
          j = j-1
 550    continue
        tp(l4) = arg
        np = np+2
        nrr = nrr+1
        do 560 i=5,ll
          j = i+nrr
          tp(j) = tp(i)+pi
 560    continue
c  restart the computations with the new set of knots.
 570  continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spherical spline.             c
c ********************************************************             c
c we have determined the number of knots and their position. we now    c
c compute the coefficients of the smoothing spline sp(teta,phi).       c
c the observation matrix a is extended by the rows of a matrix, expres-c
c sing that sp(teta,phi) must be a constant function in the variable   c
c phi and a cubic polynomial in the variable teta. the corresponding   c
c weights of these additional rows are set to 1/(p). iteratively       c
c we than have to determine the value of p such that f(p) = sum((w(i)* c
c (r(i)-sp(teta(i),phi(i))))**2)  be = s.                              c
c we already know that the least-squares polynomial corresponds to p=0,c
c and that the least-squares spherical spline corresponds to p=infin.  c
c the iteration process makes use of rational interpolation. since f(p)c
c is a convex and strictly decreasing function of p, it can be approx- c
c imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c
c three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c
c f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c
c of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<0.c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jumps of the 3-th order derivative of
c  the b-splines at the knots tt(l),l=5,...,nt-4.
 580  call fpdisc(tt,nt,5,bt,ntest)
c  evaluate the discontinuity jumps of the 3-th order derivative of
c  the b-splines at the knots tp(l),l=5,...,np-4.
      call fpdisc(tp,np,5,bp,npest)
c  initial value for p.
      p1 = 0.
      f1 = sup-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 585 i=1,ncof
        p = p+a(i,1)
 585  continue
      rn = ncof
      p = rn/p
c  find the bandwidth of the extended observation matrix.
      iband4 = iband+3
      if(ntt.le.4) iband4 = ncof
      iband3 = iband4 -1
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 920 iter=1,maxit
        pinv = one/p
c  store the triangularized observation matrix into q.
        do 600 i=1,ncof
          ff(i) = f(i)
          do 590 j=1,iband4
            q(i,j) = 0.
 590      continue
          do 600 j=1,iband
            q(i,j) = a(i,j)
 600    continue
c  extend the observation matrix with the rows of a matrix, expressing
c  that for teta=cst. sp(teta,phi) must be a constant function.
        nt6 = nt-6
        do 720 i=5,np4
          ii = i-4
          do 610 l=1,npp
             row(l) = 0.
 610      continue
          ll = ii
          do 620  l=1,5
             if(ll.gt.npp) ll=1
             row(ll) = row(ll)+bp(ii,l)
             ll = ll+1
 620      continue
          facc = 0.
          facs = 0.
          do 630 l=1,npp
             facc = facc+row(l)*coco(l)
             facs = facs+row(l)*cosi(l)
 630      continue
          do 720 j=1,nt6
c  initialize the new row.
            do 640 l=1,iband
              h(l) = 0.
 640        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            jrot = 4+(j-2)*npp
            if(j.gt.1 .and. j.lt.nt6) go to 650
            h(1) = facc
            h(2) = facs
            if(j.eq.1) jrot = 2
            go to 670
 650        do 660 l=1,npp
               h(l)=row(l)
 660        continue
 670        do 675 l=1,iband
               h(l) = h(l)*pinv
 675        continue
            ri = 0.
c  rotate the new row into triangle by givens transformations.
            do 710 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if(piv.eq.0.) then
                if (i2.le.0) go to 720
                go to 690
              endif
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
c  apply that givens transformation to the right hand side.
              call fprota(co,si,ri,ff(irot))
              if(i2.eq.0) go to 720
c  apply that givens transformation to the left hand side.
              do 680 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 680          continue
 690          do 700 l=1,i2
                h(l) = h(l+1)
 700          continue
              h(i2+1) = 0.
 710        continue
 720    continue
c  extend the observation matrix with the rows of a matrix expressing
c  that for phi=cst. sp(teta,phi) must be a cubic polynomial.
        do 810 i=5,nt4
          ii = i-4
          do 810 j=1,npp
c  initialize the new row
            do 730 l=1,iband4
              h(l) = 0.
 730        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            j1 = 1
            do 760 l=1,5
               il = ii+l
               ij = npp
               if(il.ne.3 .and. il.ne.nt4) go to 750
               j1 = j1+3-j
               j2 = j1-2
               ij = 0
               if(il.ne.3) go to 740
               j1 = 1
               j2 = 2
               ij = j+2
 740           h(j2) = bt(ii,l)*coco(j)
               h(j2+1) = bt(ii,l)*cosi(j)
 750           h(j1) = h(j1)+bt(ii,l)
               j1 = j1+ij
 760        continue
            do 765 l=1,iband4
               h(l) = h(l)*pinv
 765        continue
            ri = 0.
            jrot = 1
            if(ii.gt.2) jrot = 3+j+(ii-3)*npp
c  rotate the new row into triangle by givens transformations.
            do 800 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if(piv.eq.0.) then
                if (i2.le.0) go to 810
                go to 780
              endif
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
c  apply that givens transformation to the right hand side.
              call fprota(co,si,ri,ff(irot))
              if(i2.eq.0) go to 810
c  apply that givens transformation to the left hand side.
              do 770 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 770          continue
 780          do 790 l=1,i2
                h(l) = h(l+1)
 790          continue
              h(i2+1) = 0.
 800        continue
 810    continue
c  find dmax, the maximum value for the diagonal elements in the
c  reduced triangle.
        dmax = 0.
        do 820 i=1,ncof
          if(q(i,1).le.dmax) go to 820
          dmax = q(i,1)
 820    continue
c  check whether the matrix is rank deficient.
        sigma = eps*dmax
        do 830 i=1,ncof
          if(q(i,1).le.sigma) go to 840
 830    continue
c  backward substitution in case of full rank.
        call fpback(q,ff,ncof,iband4,c,ncc)
        rank = ncof
        go to 845
c  in case of rank deficiency, find the minimum norm solution.
 840    lwest = ncof*iband4+ncof+iband4
        if(lwrk.lt.lwest) go to 925
        lf = 1
        lh = lf+ncof
        la = lh+iband4
        call fprank(q,ff,ncof,iband4,ncc,sigma,c,sq,rank,wrk(la),
     *   wrk(lf),wrk(lh))
 845    do 850 i=1,ncof
           q(i,1) = q(i,1)/dmax
 850    continue
c  find the coefficients in the standard b-spline representation of
c  the spherical spline.
        call fprpsp(nt,np,coco,cosi,c,ff,ncoff)
c  compute f(p).
        fp = 0.
        do 890 num = 1,nreg
          num1 = num-1
          lt = num1/npp
          lp = num1-lt*npp
          jrot = lt*np4+lp
          in = index(num)
 860      if(in.eq.0) go to 890
          store = 0.
          i1 = jrot
          do 880 i=1,4
            hti = spt(in,i)
            j1 = i1
            do 870 j=1,4
              j1 = j1+1
              store = store+hti*spp(in,j)*c(j1)
 870        continue
            i1 = i1+np4
 880      continue
          fp = fp+(w(in)*(r(in)-store))**2
          in = nummer(in)
          go to 860
 890    continue
c  test whether the approximation sp(teta,phi) is an acceptable solution
        fpms = fp-s
        if(abs(fpms).le.acc) go to 980
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 940
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 900
        if((f2-f3).gt.acc) go to 895
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 920
 895    if(f2.lt.0.) ich3 = 1
 900    if(ich1.ne.0) go to 910
        if((f1-f2).gt.acc) go to 905
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 920
        if(p.ge.p3) p = p2*con1 +p3*con9
        go to 920
 905    if(f2.gt.0.) ich1 = 1
c  test whether the iteration process proceeds as theoretically
c  expected.
 910    if(f2.ge.f1 .or. f2.le.f3) go to 945
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 920  continue
c  error codes and messages.
 925  ier = lwest
      go to 990
 930  ier = 5
      go to 990
 935  ier = 4
      go to 990
 940  ier = 3
      go to 990
 945  ier = 2
      go to 990
 950  ier = 1
      go to 990
 960  ier = -2
      go to 990
 970  ier = -1
      fp = 0.
 980  if(ncof.ne.rank) ier = -rank
 990  return
      end
