      subroutine fppogr(iopt,ider,u,mu,v,mv,z,mz,z0,r,s,nuest,nvest,
     * tol,maxit,nc,nu,tu,nv,tv,c,fp,fp0,fpold,reducu,reducv,fpintu,
     * fpintv,dz,step,lastdi,nplusu,nplusv,lasttu,nru,nrv,nrdatu,
     * nrdatv,wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      integer mu,mv,mz,nuest,nvest,maxit,nc,nu,nv,lastdi,nplusu,nplusv,
     * lasttu,lwrk,ier
      real*8 z0,r,s,tol,fp,fp0,fpold,reducu,reducv,step
c  ..array arguments..
      integer iopt(3),ider(2),nrdatu(nuest),nrdatv(nvest),nru(mu),
     * nrv(mv)
      real*8 u(mu),v(mv),z(mz),tu(nuest),tv(nvest),c(nc),fpintu(nuest),
     * fpintv(nvest),dz(3),wrk(lwrk)
c  ..local scalars..
      real*8 acc,fpms,f1,f2,f3,p,per,pi,p1,p2,p3,vb,ve,zmax,zmin,rn,one,
     *
     * con1,con4,con9
      integer i,ich1,ich3,ifbu,ifbv,ifsu,ifsv,istart,iter,i1,i2,j,ju,
     * ktu,l,l1,l2,l3,l4,mpm,mumin,mu0,mu1,nn,nplu,nplv,npl1,nrintu,
     * nrintv,nue,numax,nve,nvmax
c  ..local arrays..
      integer idd(2)
      real*8 dzz(3)
c  ..function references..
      real*8 abs,datan2,fprati
      integer max0,min0
c  ..subroutine references..
c    fpknot,fpopdi
c  ..
c   set constants
      one = 1d0
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
c   initialization
      ifsu = 0
      ifsv = 0
      ifbu = 0
      ifbv = 0
      p = -one
      mumin = 4-iopt(3)
      if(ider(1).ge.0) mumin = mumin-1
      if(iopt(2).eq.1 .and. ider(2).eq.1) mumin = mumin-1
      pi = datan2(0d0,-one)
      per = pi+pi
      vb = v(1)
      ve = vb+per
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c  given a set of knots we compute the least-squares spline sinf(u,v)  c
c  and the corresponding sum of squared residuals fp = f(p=inf).       c
c  if iopt(1)=-1  sinf(u,v) is the requested approximation.            c
c  if iopt(1)>=0  we check whether we can accept the knots:            c
c    if fp <= s we will continue with the current set of knots.        c
c    if fp >  s we will increase the number of knots and compute the   c
c       corresponding least-squares spline until finally fp <= s.      c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c     knots in the u-direction equals nu=numax=mu+5+iopt(2)+iopt(3)    c
c     and in the v-direction nv=nvmax=mv+7.                            c
c    if s>0 and                                                        c
c      iopt(1)=0 we first compute the least-squares polynomial,i.e. a  c
c       spline without interior knots : nu=8 ; nv=8.                   c
c      iopt(1)=1 we start with the set of knots found at the last call c
c       of the routine, except for the case that s > fp0; then we      c
c       compute the least-squares polynomial directly.                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(iopt(1).lt.0) go to 120
c  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  numax and nvmax denote the number of knots needed for interpolation.
      numax = mu+5+iopt(2)+iopt(3)
      nvmax = mv+7
      nue = min0(numax,nuest)
      nve = min0(nvmax,nvest)
      if(s.gt.0.) go to 100
c  if s = 0, s(u,v) is an interpolating spline.
      nu = numax
      nv = nvmax
c  test whether the required storage space exceeds the available one.
      if(nu.gt.nuest .or. nv.gt.nvest) go to 420
c  find the position of the knots in the v-direction.
      do 10 l=1,mv
        tv(l+3) = v(l)
  10  continue
      tv(mv+4) = ve
      l1 = mv-2
      l2 = mv+5
      do 20 i=1,3
         tv(i) = v(l1)-per
         tv(l2) = v(i+1)+per
         l1 = l1+1
         l2 = l2+1
  20  continue
c  if not all the derivative values g(i,j) are given, we will first
c  estimate these values by computing a least-squares spline
      idd(1) = ider(1)
      if(idd(1).eq.0) idd(1) = 1
      if(idd(1).gt.0) dz(1) = z0
      idd(2) = ider(2)
      if(ider(1).lt.0) go to 30
      if(iopt(2).eq.0 .or. ider(2).ne.0) go to 70
c we set up the knots in the u-direction for computing the least-squares
c spline.
  30  i1 = 3
      i2 = mu-2
      nu = 4
      do 40 i=1,mu
         if(i1.gt.i2) go to 50
         nu = nu+1
         tu(nu) = u(i1)
         i1 = i1+2
  40  continue
  50  do 60 i=1,4
         tu(i) = 0.
         nu = nu+1
         tu(nu) = r
  60  continue
c we compute the least-squares spline for estimating the derivatives.
      call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,iopt,idd,
     *  tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,
     *  wrk,lwrk)
      ifsu = 0
c if all the derivatives at the origin are known, we compute the
c interpolating spline.
c we set up the knots in the u-direction, needed for interpolation.
  70  nn = numax-8
      if(nn.eq.0) go to 95
      ju = 2-iopt(2)
      do 80 l=1,nn
        tu(l+4) = u(ju)
        ju = ju+1
  80  continue
      nu = numax
      l = nu
      do 90 i=1,4
         tu(i) = 0.
         tu(l) = r
         l = l-1
  90  continue
c we compute the interpolating spline.
  95  call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,iopt,idd,
     *  tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,
     *  wrk,lwrk)
      go to 430
c  if s>0 our initial choice of knots depends on the value of iopt(1).
 100  ier = 0
      if(iopt(1).eq.0) go to 115
      step = -step
      if(fp0.le.s) go to 115
c  if iopt(1)=1 and fp0 > s we start computing the least-squares spline
c  according to the set of knots found at the last call of the routine.
c  we determine the number of grid coordinates u(i) inside each knot
c  interval (tu(l),tu(l+1)).
      l = 5
      j = 1
      nrdatu(1) = 0
      mu0 = 2-iopt(2)
      mu1 = mu-2+iopt(3)
      do 105 i=mu0,mu1
        nrdatu(j) = nrdatu(j)+1
        if(u(i).lt.tu(l)) go to 105
        nrdatu(j) = nrdatu(j)-1
        l = l+1
        j = j+1
        nrdatu(j) = 0
 105  continue
c  we determine the number of grid coordinates v(i) inside each knot
c  interval (tv(l),tv(l+1)).
      l = 5
      j = 1
      nrdatv(1) = 0
      do 110 i=2,mv
        nrdatv(j) = nrdatv(j)+1
        if(v(i).lt.tv(l)) go to 110
        nrdatv(j) = nrdatv(j)-1
        l = l+1
        j = j+1
        nrdatv(j) = 0
 110  continue
      idd(1) = ider(1)
      idd(2) = ider(2)
      go to 120
c  if iopt(1)=0 or iopt(1)=1 and s >= fp0,we start computing the least-
c  squares polynomial (which is a spline without interior knots).
 115  ier = -2
      idd(1) = ider(1)
      idd(2) = 1
      nu = 8
      nv = 8
      nrdatu(1) = mu-3+iopt(2)+iopt(3)
      nrdatv(1) = mv-1
      lastdi = 0
      nplusu = 0
      nplusv = 0
      fp0 = 0.
      fpold = 0.
      reducu = 0.
      reducv = 0.
c  main loop for the different sets of knots.mpm=mu+mv is a save upper
c  bound for the number of trials.
 120  mpm = mu+mv
      do 270 iter=1,mpm
c  find nrintu (nrintv) which is the number of knot intervals in the
c  u-direction (v-direction).
        nrintu = nu-7
        nrintv = nv-7
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(u,v).
        i = nu
        do 130 j=1,4
          tu(j) = 0.
          tu(i) = r
          i = i-1
 130    continue
        l1 = 4
        l2 = l1
        l3 = nv-3
        l4 = l3
        tv(l2) = vb
        tv(l3) = ve
        do 140 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tv(l2) = tv(l4)-per
          tv(l3) = tv(l1)+per
 140    continue
c  find an estimate of the range of possible values for the optimal
c  derivatives at the origin.
        ktu = nrdatu(1)+2-iopt(2)
        if(nrintu.eq.1) ktu = mu
        if(ktu.lt.mumin) ktu = mumin
        if(ktu.eq.lasttu) go to 150
         zmin = z0
         zmax = z0
         l = mv*ktu
         do 145 i=1,l
            if(z(i).lt.zmin) zmin = z(i)
            if(z(i).gt.zmax) zmax = z(i)
 145     continue
         step = zmax-zmin
         lasttu = ktu
c  find the least-squares spline sinf(u,v).
 150    call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,iopt,idd,
     *   tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,
     *   wrk,lwrk)
        if(step.lt.0.) step = -step
        if(ier.eq.(-2)) fp0 = fp
c  test whether the least-squares spline is an acceptable solution.
        if(iopt(1).lt.0) go to 440
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
c  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
c  if nu=numax and nv=nvmax, sinf(u,v) is an interpolating spline
        if(nu.eq.numax .and. nv.eq.nvmax) go to 430
c  increase the number of knots.
c  if nu=nue and nv=nve we cannot further increase the number of knots
c  because of the storage capacity limitation.
        if(nu.eq.nue .and. nv.eq.nve) go to 420
        if(ider(1).eq.0) fpintu(1) = fpintu(1)+(z0-c(1))**2
        ier = 0
c  adjust the parameter reducu or reducv according to the direction
c  in which the last added knots were located.
        if (lastdi.lt.0) go to 160
        if (lastdi.eq.0) go to 155
        go to 170
 155     nplv = 3
         idd(2) = ider(2)
         fpold = fp
         go to 230
 160    reducu = fpold-fp
        go to 175
 170    reducv = fpold-fp
c  store the sum of squared residuals for the current set of knots.
 175    fpold = fp
c  find nplu, the number of knots we should add in the u-direction.
        nplu = 1
        if(nu.eq.8) go to 180
        npl1 = nplusu*2
        rn = nplusu
        if(reducu.gt.acc) npl1 = rn*fpms/reducu
        nplu = min0(nplusu*2,max0(npl1,nplusu/2,1))
c  find nplv, the number of knots we should add in the v-direction.
 180    nplv = 3
        if(nv.eq.8) go to 190
        npl1 = nplusv*2
        rn = nplusv
        if(reducv.gt.acc) npl1 = rn*fpms/reducv
        nplv = min0(nplusv*2,max0(npl1,nplusv/2,1))
c  test whether we are going to add knots in the u- or v-direction.
 190    if (nplu.lt.nplv) go to 210
        if (nplu.eq.nplv) go to 200
        go to 230
 200    if(lastdi.lt.0) go to 230
 210    if(nu.eq.nue) go to 230
c  addition in the u-direction.
        lastdi = -1
        nplusu = nplu
        ifsu = 0
        istart = 0
        if(iopt(2).eq.0) istart = 1
        do 220 l=1,nplusu
c  add a new knot in the u-direction
          call fpknot(u,mu,tu,nu,fpintu,nrdatu,nrintu,nuest,istart)
c  test whether we cannot further increase the number of knots in the
c  u-direction.
          if(nu.eq.nue) go to 270
 220    continue
        go to 270
 230    if(nv.eq.nve) go to 210
c  addition in the v-direction.
        lastdi = 1
        nplusv = nplv
        ifsv = 0
        do 240 l=1,nplusv
c  add a new knot in the v-direction.
          call fpknot(v,mv,tv,nv,fpintv,nrdatv,nrintv,nvest,1)
c  test whether we cannot further increase the number of knots in the
c  v-direction.
          if(nv.eq.nve) go to 270
 240    continue
c  restart the computations with the new set of knots.
 270  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 300  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(u,v)                c
c *****************************************************                c
c  we have determined the number of knots and their position. we now   c
c  compute the b-spline coefficients of the smoothing spline sp(u,v).  c
c  this smoothing spline depends on the parameter p in such a way that c
c    f(p) = sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)             c
c  is a continuous, strictly decreasing function of p. moreover the    c
c  least-squares polynomial corresponds to p=0 and the least-squares   c
c  spline to p=infinity. then iteratively we have to determine the     c
c  positive value of p such that f(p)=s. the process which is proposed c
c  here makes use of rational interpolation. f(p) is approximated by a c
c  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
c  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
c  are used to calculate the new value of p such that r(p)=s.          c
c  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      dzz(1) = dz(1)
      dzz(2) = dz(2)
      dzz(3) = dz(3)
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
c  find the smoothing spline sp(u,v) and the corresponding sum f(p).
        call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dzz,iopt,idd,
     *   tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,
     *   wrk,lwrk)
c  test whether the approximation sp(u,v) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 320
        if((f2-f3).gt.acc) go to 310
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2.lt.0.) ich3 = 1
 320    if(ich1.ne.0) go to 340
        if((f1-f2).gt.acc) go to 330
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 350
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 350
c  test whether the iteration process proceeds as theoretically
c  expected.
 330    if(f2.gt.0.) ich1 = 1
 340    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 350  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end
