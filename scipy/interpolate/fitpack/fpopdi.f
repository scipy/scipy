      recursive subroutine fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,
     * mz,z0,dz,iopt,ider,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,
     * fpu,fpv,nru,nrv,wrk,lwrk)
      implicit none
c  given the set of function values z(i,j) defined on the rectangular
c  grid (u(i),v(j)),i=1,2,...,mu;j=1,2,...,mv, fpopdi determines a
c  smooth bicubic spline approximation with given knots tu(i),i=1,..,nu
c  in the u-direction and tv(j),j=1,2,...,nv in the v-direction. this
c  spline sp(u,v) will be periodic in the variable v and will satisfy
c  the following constraints
c
c     s(tu(1),v) = dz(1) , tv(4) <=v<= tv(nv-3)
c
c  and (if iopt(2) = 1)
c
c     d s(tu(1),v)
c     ------------ =  dz(2)*cos(v)+dz(3)*sin(v) , tv(4) <=v<= tv(nv-3)
c     d u
c
c  and (if iopt(3) = 1)
c
c     s(tu(nu),v)  =  0   tv(4) <=v<= tv(nv-3)
c
c  where the parameters dz(i) correspond to the derivative values g(i,j)
c  as defined in subroutine pogrid.
c
c  the b-spline coefficients of sp(u,v) are determined as the least-
c  squares solution  of an overdetermined linear system which depends
c  on the value of p and on the values dz(i),i=1,2,3. the correspond-
c  ing sum of squared residuals sq is a simple quadratic function in
c  the variables dz(i). these may or may not be provided. the values
c  dz(i) which are not given will be determined so as to minimize the
c  resulting sum of squared residuals sq. in that case the user must
c  provide some initial guess dz(i) and some estimate (dz(i)-step,
c  dz(i)+step) of the range of possible values for these latter.
c
c  sp(u,v) also depends on the parameter p (p>0) in such a way that
c    - if p tends to infinity, sp(u,v) becomes the least-squares spline
c      with given knots, satisfying the constraints.
c    - if p tends to zero, sp(u,v) becomes the least-squares polynomial,
c      satisfying the constraints.
c    - the function  f(p)=sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)
c      is continuous and strictly decreasing for p>0.
c
c  ..scalar arguments..
      integer ifsu,ifsv,ifbu,ifbv,mu,mv,mz,nu,nv,nuest,nvest,
     * nc,lwrk
      real*8 z0,p,step,fp
c  ..array arguments..
      integer ider(2),nru(mu),nrv(mv),iopt(3)
      real*8 u(mu),v(mv),z(mz),dz(3),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv)
     *,
     * wrk(lwrk)
c  ..local scalars..
      real*8 res,sq,sqq,step1,step2,three
      integer i,id0,iop0,iop1,i1,j,l,laa,lau,lav1,lav2,lbb,lbu,lbv,
     * lcc,lcs,lq,lri,lsu,lsv,l1,l2,mm,mvnu,number
c  ..local arrays..
      integer nr(3)
      real*8 delta(3),dzz(3),sum(3),a(6,6),g(6)
c  ..function references..
      integer max0
c  ..subroutine references..
c    fpgrdi,fpsysy
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
      laa = lbv+5*nvest
      lbb = laa+2*mv
      lcc = lbb+2*nvest
      lcs = lcc+nvest
c  we calculate the smoothing spline sp(u,v) according to the input
c  values dz(i),i=1,2,3.
      iop0 = iopt(2)
      iop1 = iopt(3)
      call fpgrdi(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,z,mz,dz,
     * iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,
     * wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     * wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),
     * wrk(lcc),wrk(lcs),nru,nrv)
      id0 = ider(1)
      if(id0.ne.0) go to 5
      res = (z0-dz(1))**2
      fp = fp+res
      sq = sq+res
c in case all derivative values dz(i) are given (step<=0) or in case
c we have spline interpolation, we accept this spline as a solution.
  5   if(step.le.0. .or. sq.le.0.) return
      dzz(1) = dz(1)
      dzz(2) = dz(2)
      dzz(3) = dz(3)
c number denotes the number of derivative values dz(i) that still must
c be optimized. let us denote these parameters by g(j),j=1,...,number.
      number = 0
      if(id0.gt.0) go to 10
      number = 1
      nr(1) = 1
      delta(1) = step
  10  if(iop0.eq.0) go to 20
      if(ider(2).ne.0) go to 20
      step2 = step*three/tu(5)
      nr(number+1) = 2
      nr(number+2) = 3
      delta(number+1) = step2
      delta(number+2) = step2
      number = number+2
  20  if(number.eq.0) return
c the sum of squared residuals sq is a quadratic polynomial in the
c parameters g(j). we determine the unknown coefficients of this
c polymomial by calculating (number+1)*(number+2)/2 different splines
c according to specific values for g(j).
      do 30 i=1,number
         l = nr(i)
         step1 = delta(i)
         dzz(l) = dz(l)+step1
         call fpgrdi(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,z,mz,dzz,
     *    iop0,iop1,tu,nu,tv,nv,p,c,nc,sum(i),fp,fpu,fpv,mm,mvnu,
     *    wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *    wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),
     *    wrk(lcc),wrk(lcs),nru,nrv)
         if(id0.eq.0) sum(i) = sum(i)+(z0-dzz(1))**2
         dzz(l) = dz(l)-step1
         call fpgrdi(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,z,mz,dzz,
     *    iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu,
     *    wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *    wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),
     *    wrk(lcc),wrk(lcs),nru,nrv)
         if(id0.eq.0) sqq = sqq+(z0-dzz(1))**2
         a(i,i) = (sum(i)+sqq-sq-sq)/step1**2
         if(a(i,i).le.0.) go to 80
         g(i) = (sqq-sum(i))/(step1+step1)
         dzz(l) = dz(l)
  30  continue
      if(number.eq.1) go to 60
      do 50 i=2,number
         l1 = nr(i)
         step1 = delta(i)
         dzz(l1) = dz(l1)+step1
         i1 = i-1
         do 40 j=1,i1
            l2 = nr(j)
            step2 = delta(j)
            dzz(l2) = dz(l2)+step2
            call fpgrdi(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,z,mz,dzz,
     *       iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu,
     *       wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *       wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),
     *       wrk(lcc),wrk(lcs),nru,nrv)
            if(id0.eq.0) sqq = sqq+(z0-dzz(1))**2
            a(i,j) = (sq+sqq-sum(i)-sum(j))/(step1*step2)
            dzz(l2) = dz(l2)
  40     continue
         dzz(l1) = dz(l1)
  50  continue
c the optimal values g(j) are found as the solution of the system
c d (sq) / d (g(j)) = 0 , j=1,...,number.
  60  call fpsysy(a,number,g)
      do 70 i=1,number
         l = nr(i)
         dz(l) = dz(l)+g(i)
  70  continue
c we determine the spline sp(u,v) according to the optimal values g(j).
  80  call fpgrdi(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,z,mz,dz,
     * iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,
     * wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     * wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),
     * wrk(lcc),wrk(lcs),nru,nrv)
      if(id0.eq.0) fp = fp+(z0-dz(1))**2
      return
      end
