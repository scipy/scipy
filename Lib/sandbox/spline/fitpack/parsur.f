      subroutine parsur(iopt,ipar,idim,mu,u,mv,v,f,s,nuest,nvest,
     * nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c  given the set of ordered points f(i,j) in the idim-dimensional space,
c  corresponding to grid values (u(i),v(j)) ,i=1,...,mu ; j=1,...,mv,
c  parsur determines a smooth approximating spline surface s(u,v) , i.e.
c    f1 = s1(u,v)
c      ...                u(1) <= u <= u(mu) ; v(1) <= v <= v(mv)
c    fidim = sidim(u,v)
c  with sl(u,v), l=1,2,...,idim bicubic spline functions with common
c  knots tu(i),i=1,...,nu in the u-variable and tv(j),j=1,...,nv in the
c  v-variable.
c  in addition, these splines will be periodic in the variable u if
c  ipar(1) = 1 and periodic in the variable v if ipar(2) = 1.
c  if iopt=-1, parsur determines the least-squares bicubic spline
c  surface according to a given set of knots.
c  if iopt>=0, the number of knots of s(u,v) and their position
c  is chosen automatically by the routine. the smoothness of s(u,v) is
c  achieved by minimalizing the discontinuity jumps of the derivatives
c  of the splines at the knots. the amount of smoothness of s(u,v) is
c  determined by the condition that
c  fp=sumi=1,mu(sumj=1,mv(dist(f(i,j)-s(u(i),v(j)))**2))<=s,
c  with s a given non-negative constant.
c  the fit s(u,v) is given in its b-spline representation and can be
c  evaluated by means of routine surev.
c
c calling sequence:
c     call parsur(iopt,ipar,idim,mu,u,mv,v,f,s,nuest,nvest,nu,tu,
c    *  nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer flag. unchanged on exit.
c          on entry iopt must specify whether a least-squares surface
c          (iopt=-1) or a smoothing surface (iopt=0 or 1)must be
c          determined.
c          if iopt=0 the routine will start with the initial set of
c          knots needed for determining the least-squares polynomial
c          surface.
c          if iopt=1 the routine will continue with the set of knots
c          found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately
c          preceded by another call with iopt = 1 or iopt = 0.
c  ipar  : integer array of dimension 2. unchanged on exit.
c          on entry ipar(1) must specify whether (ipar(1)=1) or not
c          (ipar(1)=0) the splines must be periodic in the variable u.
c          on entry ipar(2) must specify whether (ipar(2)=1) or not
c          (ipar(2)=0) the splines must be periodic in the variable v.
c  idim  : integer. on entry idim must specify the dimension of the
c          surface. 1 <= idim <= 3. unchanged on exit.
c  mu    : integer. on entry mu must specify the number of grid points
c          along the u-axis. unchanged on exit.
c          mu >= mumin where mumin=4-2*ipar(1)
c  u     : real array of dimension at least (mu). before entry, u(i)
c          must be set to the u-co-ordinate of the i-th grid point
c          along the u-axis, for i=1,2,...,mu. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c  mv    : integer. on entry mv must specify the number of grid points
c          along the v-axis. unchanged on exit.
c          mv >= mvmin where mvmin=4-2*ipar(2)
c  v     : real array of dimension at least (mv). before entry, v(j)
c          must be set to the v-co-ordinate of the j-th grid point
c          along the v-axis, for j=1,2,...,mv. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c  f     : real array of dimension at least (mu*mv*idim).
c          before entry, f(mu*mv*(l-1)+mv*(i-1)+j) must be set to the
c          l-th co-ordinate of the data point corresponding to the
c          the grid point (u(i),v(j)) for l=1,...,idim ,i=1,...,mu
c          and j=1,...,mv. unchanged on exit.
c          if ipar(1)=1 it is expected that f(mu*mv*(l-1)+mv*(mu-1)+j)
c          = f(mu*mv*(l-1)+j), l=1,...,idim ; j=1,...,mv
c          if ipar(2)=1 it is expected that f(mu*mv*(l-1)+mv*(i-1)+mv)
c          = f(mu*mv*(l-1)+mv*(i-1)+1), l=1,...,idim ; i=1,...,mu
c  s     : real. on entry (if iopt>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nuest : integer. unchanged on exit.
c  nvest : integer. unchanged on exit.
c          on entry, nuest and nvest must specify an upper bound for the
c          number of knots required in the u- and v-directions respect.
c          these numbers will also determine the storage space needed by
c          the routine. nuest >= 8, nvest >= 8.
c          in most practical situation nuest = mu/2, nvest=mv/2, will
c          be sufficient. always large enough are nuest=mu+4+2*ipar(1),
c          nvest = mv+4+2*ipar(2), the number of knots needed for
c          interpolation (s=0). see also further comments.
c  nu    : integer.
c          unless ier=10 (in case iopt>=0), nu will contain the total
c          number of knots with respect to the u-variable, of the spline
c          surface returned. if the computation mode iopt=1 is used,
c          the value of nu should be left unchanged between subsequent
c          calls. in case iopt=-1, the value of nu should be specified
c          on entry.
c  tu    : real array of dimension at least (nuest).
c          on succesful exit, this array will contain the knots of the
c          splines with respect to the u-variable, i.e. the position of
c          the interior knots tu(5),...,tu(nu-4) as well as the position
c          of the additional knots tu(1),...,tu(4) and tu(nu-3),...,
c          tu(nu) needed for the b-spline representation.
c          if the computation mode iopt=1 is used,the values of tu(1)
c          ...,tu(nu) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tu(5),
c          ...tu(nu-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  nv    : integer.
c          unless ier=10 (in case iopt>=0), nv will contain the total
c          number of knots with respect to the v-variable, of the spline
c          surface returned. if the computation mode iopt=1 is used,
c          the value of nv should be left unchanged between subsequent
c          calls. in case iopt=-1, the value of nv should be specified
c          on entry.
c  tv    : real array of dimension at least (nvest).
c          on succesful exit, this array will contain the knots of the
c          splines with respect to the v-variable, i.e. the position of
c          the interior knots tv(5),...,tv(nv-4) as well as the position
c          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
c          tv(nv) needed for the b-spline representation.
c          if the computation mode iopt=1 is used,the values of tv(1)
c          ...,tv(nv) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tv(5),
c          ...tv(nv-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : real array of dimension at least (nuest-4)*(nvest-4)*idim.
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(u,v)
c  fp    : real. unless ier=10, fp contains the sum of squared
c          residuals of the spline surface returned.
c  wrk   : real array of dimension (lwrk). used as workspace.
c          if the computation mode iopt=1 is used the values of
c          wrk(1),...,wrk(4) should be left unchanged between subsequent
c          calls.
c  lwrk  : integer. on entry lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.
c          lwrk must not be too small.
c           lwrk >= 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))+
c           4*(mu+mv)+q*idim where q is the larger of mv and nuest.
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c          if the computation mode iopt=1 is used the values of
c          iwrk(1),.,iwrk(3) should be left unchanged between subsequent
c          calls.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= 3+mu+mv+nuest+nvest.
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the surface returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline surface returned is an
c            interpolating surface (fp=0).
c   ier=-2 : normal return. the surface returned is the least-squares
c            polynomial surface. in this extreme case fp gives the
c            upper bound for the smoothing factor s.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nuest and
c            nvest.
c            probably causes : nuest or nvest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the least-squares surface
c            according to the current set of knots. the parameter fp
c            gives the corresponding sum of squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing surface with
c            fp = s. probably causes : s too small.
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing surface
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt<=1, 0<=ipar(1)<=1, 0<=ipar(2)<=1, 1 <=idim<=3
c            mu >= 4-2*ipar(1),mv >= 4-2*ipar(2), nuest >=8, nvest >= 8,
c            kwrk>=3+mu+mv+nuest+nvest,
c            lwrk >= 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))
c             +4*(mu+mv)+max(nuest,mv)*idim
c            u(i-1)<u(i),i=2,..,mu, v(i-1)<v(i),i=2,...,mv
c            if iopt=-1: 8<=nu<=min(nuest,mu+4+2*ipar(1))
c                        u(1)<tu(5)<tu(6)<...<tu(nu-4)<u(mu)
c                        8<=nv<=min(nvest,mv+4+2*ipar(2))
c                        v(1)<tv(5)<tv(6)<...<tv(nv-4)<v(mv)
c                    the schoenberg-whitney conditions, i.e. there must
c                    be subset of grid co-ordinates uu(p) and vv(q) such
c                    that   tu(p) < uu(p) < tu(p+4) ,p=1,...,nu-4
c                           tv(q) < vv(q) < tv(q+4) ,q=1,...,nv-4
c                     (see fpchec or fpchep)
c            if iopt>=0: s>=0
c                       if s=0: nuest>=mu+4+2*ipar(1)
c                               nvest>=mv+4+2*ipar(2)
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c
c further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the surface will be too smooth and signal will be
c   lost ; if s is too small the surface will pick up too much noise. in
c   the extreme cases the program will return an interpolating surface
c   if s=0 and the constrained least-squares polynomial surface if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the accuracy of the data values.
c   if the user has an idea of the statistical errors on the data, he
c   can also find a proper estimate for s. for, by assuming that, if he
c   specifies the right s, parsur will return a surface s(u,v) which
c   exactly reproduces the surface underlying the data he can evaluate
c   the sum(dist(f(i,j)-s(u(i),v(j)))**2) to find a good estimate for s.
c   for example, if he knows that the statistical errors on his f(i,j)-
c   values is not greater than 0.1, he may expect that a good s should
c   have a value not larger than mu*mv*(0.1)**2.
c   if nothing is known about the statistical error in f(i,j), s must
c   be determined by trial and error, taking account of the comments
c   above. the best is then to start with a very large value of s (to
c   determine the le-sq polynomial surface and the corresponding upper
c   bound fp0 for s) and then to progressively decrease the value of s
c   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
c   and more carefully as the approximation shows more detail) to
c   obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt = 1 the program will continue with the knots found at
c   the last call of the routine. this will save a lot of computation
c   time if parsur is called repeatedly for different values of s.
c   the number of knots of the surface returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   surface underlying the data. if the computation mode iopt = 1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1,the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   parsur once more with the chosen value for s but now with iopt=0.
c   indeed, parsur may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nuest and
c   nvest. indeed, if at a certain stage in parsur the number of knots
c   in one direction (say nu) has reached the value of its upper bound
c   (nuest), then from that moment on all subsequent knots are added
c   in the other (v) direction. this may indicate that the value of
c   nuest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting nuest=8 (the lowest allowable value for
c   nuest), the user can indicate that he wants an approximation with
c   splines which are simple cubic polynomials in the variable u.
c
c  other subroutines required:
c    fppasu,fpchec,fpchep,fpknot,fprati,fpgrpa,fptrnp,fpback,
c    fpbacp,fpbspl,fptrpe,fpdisc,fpgivs,fprota
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      real*8 s,fp
      integer iopt,idim,mu,mv,nuest,nvest,nu,nv,lwrk,kwrk,ier
c  ..array arguments..
      real*8 u(mu),v(mv),f(mu*mv*idim),tu(nuest),tv(nvest),
     * c((nuest-4)*(nvest-4)*idim),wrk(lwrk)
      integer ipar(2),iwrk(kwrk)
c  ..local scalars..
      real*8 tol,ub,ue,vb,ve,peru,perv
      integer i,j,jwrk,kndu,kndv,knru,knrv,kwest,l1,l2,l3,l4,
     * lfpu,lfpv,lwest,lww,maxit,nc,mf,mumin,mvmin
c  ..function references..
      integer max0
c  ..subroutine references..
c    fppasu,fpchec,fpchep
c  ..
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 200
      if(ipar(1).lt.0 .or. ipar(1).gt.1) go to 200
      if(ipar(2).lt.0 .or. ipar(2).gt.1) go to 200
      if(idim.le.0 .or. idim.gt.3) go to 200
      mumin = 4-2*ipar(1)
      if(mu.lt.mumin .or. nuest.lt.8) go to 200
      mvmin = 4-2*ipar(2)
      if(mv.lt.mvmin .or. nvest.lt.8) go to 200
      mf = mu*mv
      nc = (nuest-4)*(nvest-4)
      lwest = 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))+
     * 4*(mu+mv)+max0(nuest,mv)*idim
      kwest = 3+mu+mv+nuest+nvest
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 200
      do 10 i=2,mu
        if(u(i-1).ge.u(i)) go to 200
  10  continue
      do 20 i=2,mv
        if(v(i-1).ge.v(i)) go to 200
  20  continue
      if(iopt.ge.0) go to 100
      if(nu.lt.8 .or. nu.gt.nuest) go to 200
      ub = u(1)
      ue = u(mu)
      if (ipar(1).ne.0) go to 40
      j = nu
      do 30 i=1,4
        tu(i) = ub
        tu(j) = ue
        j = j-1
  30  continue
      call fpchec(u,mu,tu,nu,3,ier)
      if(ier.ne.0) go to 200
      go to 60
  40  l1 = 4
      l2 = l1
      l3 = nu-3
      l4 = l3
      peru = ue-ub
      tu(l2) = ub
      tu(l3) = ue
      do 50 j=1,3
        l1 = l1+1
        l2 = l2-1
        l3 = l3+1
        l4 = l4-1
        tu(l2) = tu(l4)-peru
        tu(l3) = tu(l1)+peru
  50  continue
      call fpchep(u,mu,tu,nu,3,ier)
      if(ier.ne.0) go to 200
  60  if(nv.lt.8 .or. nv.gt.nvest) go to 200
      vb = v(1)
      ve = v(mv)
      if (ipar(2).ne.0) go to 80
      j = nv
      do 70 i=1,4
        tv(i) = vb
        tv(j) = ve
        j = j-1
  70  continue
      call fpchec(v,mv,tv,nv,3,ier)
      if(ier.ne.0) go to 200
      go to 150
  80  l1 = 4
      l2 = l1
      l3 = nv-3
      l4 = l3
      perv = ve-vb
      tv(l2) = vb
      tv(l3) = ve
      do 90 j=1,3
        l1 = l1+1
        l2 = l2-1
        l3 = l3+1
        l4 = l4-1
        tv(l2) = tv(l4)-perv
        tv(l3) = tv(l1)+perv
  90  continue
      call fpchep(v,mv,tv,nv,3,ier)
      if (ier.eq.0) go to 150
      go to 200
 100  if(s.lt.0.) go to 200
      if(s.eq.0. .and. (nuest.lt.(mu+4+2*ipar(1)) .or.
     * nvest.lt.(mv+4+2*ipar(2))) )go to 200
      ier = 0
c  we partition the working space and determine the spline approximation
 150  lfpu = 5
      lfpv = lfpu+nuest
      lww = lfpv+nvest
      jwrk = lwrk-4-nuest-nvest
      knru = 4
      knrv = knru+mu
      kndu = knrv+mv
      kndv = kndu+nuest
      call fppasu(iopt,ipar,idim,u,mu,v,mv,f,mf,s,nuest,nvest,
     * tol,maxit,nc,nu,tu,nv,tv,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),
     * wrk(lfpu),wrk(lfpv),iwrk(1),iwrk(2),iwrk(3),iwrk(knru),
     * iwrk(knrv),iwrk(kndu),iwrk(kndv),wrk(lww),jwrk,ier)
 200  return
      end

