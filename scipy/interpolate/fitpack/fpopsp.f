      recursive subroutine fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,
     * mr,r0,r1,dr,iopt,ider,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,
     * fp,fpu,fpv,nru,nrv,wrk,lwrk)
      implicit none
c  given the set of function values r(i,j) defined on the rectangular
c  grid (u(i),v(j)),i=1,2,...,mu;j=1,2,...,mv, fpopsp determines a
c  smooth bicubic spline approximation with given knots tu(i),i=1,..,nu
c  in the u-direction and tv(j),j=1,2,...,nv in the v-direction. this
c  spline sp(u,v) will be periodic in the variable v and will satisfy
c  the following constraints
c
c     s(tu(1),v) = dr(1) , tv(4) <=v<= tv(nv-3)
c
c     s(tu(nu),v) = dr(4) , tv(4) <=v<= tv(nv-3)
c
c  and (if iopt(2) = 1)
c
c     d s(tu(1),v)
c     ------------ =  dr(2)*cos(v)+dr(3)*sin(v) , tv(4) <=v<= tv(nv-3)
c     d u
c
c  and (if iopt(3) = 1)
c
c     d s(tu(nu),v)
c     ------------- =  dr(5)*cos(v)+dr(6)*sin(v) , tv(4) <=v<= tv(nv-3)
c     d u
c
c  where the parameters dr(i) correspond to the derivative values at the
c  poles as defined in subroutine spgrid.
c
c  the b-spline coefficients of sp(u,v) are determined as the least-
c  squares solution  of an overdetermined linear system which depends
c  on the value of p and on the values dr(i),i=1,...,6. the correspond-
c  ing sum of squared residuals sq is a simple quadratic function in
c  the variables dr(i). these may or may not be provided. the values
c  dr(i) which are not given will be determined so as to minimize the
c  resulting sum of squared residuals sq. in that case the user must
c  provide some initial guess dr(i) and some estimate (dr(i)-step,
c  dr(i)+step) of the range of possible values for these latter.
c
c  sp(u,v) also depends on the parameter p (p>0) in such a way that
c    - if p tends to infinity, sp(u,v) becomes the least-squares spline
c      with given knots, satisfying the constraints.
c    - if p tends to zero, sp(u,v) becomes the least-squares polynomial,
c      satisfying the constraints.
c    - the function  f(p)=sumi=1,mu(sumj=1,mv((r(i,j)-sp(u(i),v(j)))**2)
c      is continuous and strictly decreasing for p>0.
c
c  ..scalar arguments..
      integer ifsu,ifsv,ifbu,ifbv,mu,mv,mr,nu,nv,nuest,nvest,
     * nc,lwrk
      real*8 r0,r1,p,fp
c  ..array arguments..
      integer ider(4),nru(mu),nrv(mv),iopt(3)
      real*8 u(mu),v(mv),r(mr),dr(6),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv)
     *,
     * wrk(lwrk),step(2)
c  ..local scalars..
      real*8 sq,sqq,sq0,sq1,step1,step2,three
      integer i,id0,iop0,iop1,i1,j,l,lau,lav1,lav2,la0,la1,lbu,lbv,lb0,
     * lb1,lc0,lc1,lcs,lq,lri,lsu,lsv,l1,l2,mm,mvnu,number, id1
c  ..local arrays..
      integer nr(6)
      real*8 delta(6),drr(6),sum(6),a(6,6),g(6)
c  ..function references..
      integer max0
c  ..subroutine references..
c    fpgrsp,fpsysy
c  ..
c  set constant
      three = 3
c  we partition the working space
      lsu = 1
      lsv = lsu+4*mu
      lri = lsv+4*mv
      mm = max0(nuest,mv+nvest)
      lq = lri+mm
      mvnu = nuest*(mv+nvest-8)
      lau = lq+mvnu
      lav1 = lau+5*nuest
      lav2 = lav1+6*nvest
      lbu = lav2+4*nvest
      lbv = lbu+5*nuest
      la0 = lbv+5*nvest
      la1 = la0+2*mv
      lb0 = la1+2*mv
      lb1 = lb0+2*nvest
      lc0 = lb1+2*nvest
      lc1 = lc0+nvest
      lcs = lc1+nvest
c  we calculate the smoothing spline sp(u,v) according to the input
c  values dr(i),i=1,...,6.
      iop0 = iopt(2)
      iop1 = iopt(3)
      id0 = ider(1)
      id1 = ider(3)
      call fpgrsp(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,r,mr,dr,
     * iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,
     * wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     * wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),
     * wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
      sq0 = 0.
      sq1 = 0.
      if(id0.eq.0) sq0 = (r0-dr(1))**2
      if(id1.eq.0) sq1 = (r1-dr(4))**2
      sq = sq+sq0+sq1
c in case all derivative values dr(i) are given (step<=0) or in case
c we have spline interpolation, we accept this spline as a solution.
      if(sq.le.0.) return
      if(step(1).le.0. .and. step(2).le.0.) return
      do 10 i=1,6
        drr(i) = dr(i)
  10  continue
c number denotes the number of derivative values dr(i) that still must
c be optimized. let us denote these parameters by g(j),j=1,...,number.
      number = 0
      if(id0.gt.0) go to 20
      number = 1
      nr(1) = 1
      delta(1) = step(1)
  20  if(iop0.eq.0) go to 30
      if(ider(2).ne.0) go to 30
      step2 = step(1)*three/(tu(5)-tu(4))
      nr(number+1) = 2
      nr(number+2) = 3
      delta(number+1) = step2
      delta(number+2) = step2
      number = number+2
  30  if(id1.gt.0) go to 40
      number = number+1
      nr(number) = 4
      delta(number) = step(2)
  40  if(iop1.eq.0) go to 50
      if(ider(4).ne.0) go to 50
      step2 = step(2)*three/(tu(nu)-tu(nu-4))
      nr(number+1) = 5
      nr(number+2) = 6
      delta(number+1) = step2
      delta(number+2) = step2
      number = number+2
  50  if(number.eq.0) return
c the sum of squared residulas sq is a quadratic polynomial in the
c parameters g(j). we determine the unknown coefficients of this
c polymomial by calculating (number+1)*(number+2)/2 different splines
c according to specific values for g(j).
      do 60 i=1,number
         l = nr(i)
         step1 = delta(i)
         drr(l) = dr(l)+step1
         call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr,
     *    iop0,iop1,tu,nu,tv,nv,p,c,nc,sum(i),fp,fpu,fpv,mm,mvnu,
     *    wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *    wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),
     *    wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
         if(id0.eq.0) sq0 = (r0-drr(1))**2
         if(id1.eq.0) sq1 = (r1-drr(4))**2
         sum(i) = sum(i)+sq0+sq1
         drr(l) = dr(l)-step1
         call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr,
     *    iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu,
     *    wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *    wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),
     *    wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
         if(id0.eq.0) sq0 = (r0-drr(1))**2
         if(id1.eq.0) sq1 = (r1-drr(4))**2
         sqq = sqq+sq0+sq1
         drr(l) = dr(l)
         a(i,i) = (sum(i)+sqq-sq-sq)/step1**2
         if(a(i,i).le.0.) go to 110
         g(i) = (sqq-sum(i))/(step1+step1)
  60  continue
      if(number.eq.1) go to 90
      do 80 i=2,number
         l1 = nr(i)
         step1 = delta(i)
         drr(l1) = dr(l1)+step1
         i1 = i-1
         do 70 j=1,i1
            l2 = nr(j)
            step2 = delta(j)
            drr(l2) = dr(l2)+step2
            call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr,
     *       iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu,
     *       wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *       wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),
     *       wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
            if(id0.eq.0) sq0 = (r0-drr(1))**2
            if(id1.eq.0) sq1 = (r1-drr(4))**2
            sqq = sqq+sq0+sq1
            a(i,j) = (sq+sqq-sum(i)-sum(j))/(step1*step2)
            drr(l2) = dr(l2)
  70     continue
         drr(l1) = dr(l1)
  80  continue
c the optimal values g(j) are found as the solution of the system
c d (sq) / d (g(j)) = 0 , j=1,...,number.
  90  call fpsysy(a,number,g)
      do 100 i=1,number
         l = nr(i)
         dr(l) = dr(l)+g(i)
 100  continue
c we determine the spline sp(u,v) according to the optimal values g(j).
 110  call fpgrsp(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,r,mr,dr,
     * iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,
     * wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     * wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),
     * wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
      if(id0.eq.0) sq0 = (r0-dr(1))**2
      if(id1.eq.0) sq1 = (r1-dr(4))**2
      sq = sq+sq0+sq1
      return
      end
