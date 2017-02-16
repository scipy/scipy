      subroutine sphere(iopt,m,teta,phi,r,w,s,ntest,npest,eps,
     *  nt,tt,np,tp,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
      implicit none
c  subroutine sphere determines a smooth bicubic spherical spline
c  approximation s(teta,phi), 0 <= teta <= pi ; 0 <= phi <= 2*pi
c  to a given set of data points (teta(i),phi(i),r(i)),i=1,2,...,m.
c  such a spline has the following specific properties
c
c    (1) s(0,phi)  = constant   0 <=phi<= 2*pi.
c
c    (2) s(pi,phi) = constant   0 <=phi<= 2*pi
c
c         j             j
c        d s(teta,0)   d s(teta,2*pi)
c    (3) ----------- = ------------   0 <=teta<=pi, j=0,1,2
c             j             j
c        d phi         d phi
c
c        d s(0,phi)    d s(0,0)             d s(0,pi/2)
c    (4) ----------  = -------- *cos(phi) + ----------- *sin(phi)
c        d teta        d teta               d teta
c
c        d s(pi,phi)   d s(pi,0)            d s(pi,pi/2)
c    (5) ----------- = ---------*cos(phi) + ------------*sin(phi)
c        d teta        d teta               d teta
c
c  if iopt =-1 sphere calculates a weighted least-squares spherical
c  spline according to a given set of knots in teta- and phi- direction.
c  if iopt >=0, the number of knots in each direction and their position
c  tt(j),j=1,2,...,nt ; tp(j),j=1,2,...,np are chosen automatically by
c  the routine. the smoothness of s(teta,phi) is then achieved by mini-
c  malizing the discontinuity jumps of the derivatives of the spline
c  at the knots. the amount of smoothness of s(teta,phi) is determined
c  by the condition that fp = sum((w(i)*(r(i)-s(teta(i),phi(i))))**2)
c  be <= s, with s a given non-negative constant.
c  the spherical spline is given in the standard b-spline representation
c  of bicubic splines and can be evaluated by means of subroutine bispev
c
c calling sequence:
c     call sphere(iopt,m,teta,phi,r,w,s,ntest,npest,eps,
c    *  nt,tt,np,tp,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer flag. on entry iopt must specify whether a weighted
c          least-squares spherical spline (iopt=-1) or a smoothing
c          spherical spline (iopt=0 or 1) must be determined.
c          if iopt=0 the routine will start with an initial set of knots
c          tt(i)=0,tt(i+4)=pi,i=1,...,4;tp(i)=0,tp(i+4)=2*pi,i=1,...,4.
c          if iopt=1 the routine will continue with the set of knots
c          found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately pre-
c                     ceded by another call with iopt=1 or iopt=0.
c          unchanged on exit.
c  m     : integer. on entry m must specify the number of data points.
c          m >= 2. unchanged on exit.
c  teta  : real array of dimension at least (m).
c  phi   : real array of dimension at least (m).
c  r     : real array of dimension at least (m).
c          before entry,teta(i),phi(i),r(i) must be set to the spherical
c          co-ordinates of the i-th data point, for i=1,...,m.the order
c          of the data points is immaterial. unchanged on exit.
c  w     : real array of dimension at least (m). before entry, w(i) must
c          be set to the i-th value in the set of weights. the w(i) must
c          be strictly positive. unchanged on exit.
c  s     : real. on entry (in case iopt>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  ntest : integer. unchanged on exit.
c  npest : integer. unchanged on exit.
c          on entry, ntest and npest must specify an upper bound for the
c          number of knots required in the teta- and phi-directions.
c          these numbers will also determine the storage space needed by
c          the routine. ntest >= 8, npest >= 8.
c          in most practical situation ntest = npest = 8+sqrt(m/2) will
c          be sufficient. see also further comments.
c  eps   : real.
c          on entry, eps must specify a threshold for determining the
c          effective rank of an over-determined linear system of equat-
c          ions. 0 < eps < 1.  if the number of decimal digits in the
c          computer representation of a real number is q, then 10**(-q)
c          is a suitable value for eps in most practical applications.
c          unchanged on exit.
c  nt    : integer.
c          unless ier=10 (in case iopt >=0), nt will contain the total
c          number of knots with respect to the teta-variable, of the
c          spline approximation returned. if the computation mode iopt=1
c          is used, the value of nt should be left unchanged between
c          subsequent calls.
c          in case iopt=-1, the value of nt should be specified on entry
c  tt    : real array of dimension at least ntest.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the teta-variable, i.e. the position
c          of the interior knots tt(5),...,tt(nt-4) as well as the
c          position of the additional knots tt(1)=...=tt(4)=0 and
c          tt(nt-3)=...=tt(nt)=pi needed for the b-spline representation
c          if the computation mode iopt=1 is used, the values of tt(1),
c          ...,tt(nt) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tt(5),
c          ...tt(nt-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  np    : integer.
c          unless ier=10 (in case iopt >=0), np will contain the total
c          number of knots with respect to the phi-variable, of the
c          spline approximation returned. if the computation mode iopt=1
c          is used, the value of np should be left unchanged between
c          subsequent calls.
c          in case iopt=-1, the value of np (>=9) should be specified
c          on entry.
c  tp    : real array of dimension at least npest.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the phi-variable, i.e. the position of
c          the interior knots tp(5),...,tp(np-4) as well as the position
c          of the additional knots tp(1),...,tp(4) and tp(np-3),...,
c          tp(np) needed for the b-spline representation.
c          if the computation mode iopt=1 is used, the values of tp(1),
c          ...,tp(np) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tp(5),
c          ...tp(np-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : real array of dimension at least (ntest-4)*(npest-4).
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(teta,phi).
c  fp    : real. unless ier=10, fp contains the weighted sum of
c          squared residuals of the spline approximation returned.
c  wrk1  : real array of dimension (lwrk1). used as workspace.
c          if the computation mode iopt=1 is used the value of wrk1(1)
c          should be left unchanged between subsequent calls.
c          on exit wrk1(2),wrk1(3),...,wrk1(1+ncof) will contain the
c          values d(i)/max(d(i)),i=1,...,ncof=6+(np-7)*(nt-8)
c          with d(i) the i-th diagonal element of the reduced triangular
c          matrix for calculating the b-spline coefficients. it includes
c          those elements whose square is less than eps,which are treat-
c          ed as 0 in the case of presumed rank deficiency (ier<-2).
c  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
c          the array wrk1 as declared in the calling (sub)program.
c          lwrk1 must not be too small. let
c            u = ntest-7, v = npest-7, then
c          lwrk1 >= 185+52*v+10*u+14*u*v+8*(u-1)*v**2+8*m
c  wrk2  : real array of dimension (lwrk2). used as workspace, but
c          only in the case a rank deficient system is encountered.
c  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
c          the array wrk2 as declared in the calling (sub)program.
c          lwrk2 > 0 . a save upper bound  for lwrk2 = 48+21*v+7*u*v+
c          4*(u-1)*v**2 where u,v are as above. if there are enough data
c          points, scattered uniformly over the approximation domain
c          and if the smoothing factor s is not too small, there is a
c          good chance that this extra workspace is not needed. a lot
c          of memory might therefore be saved by setting lwrk2=1.
c          (see also ier > 10)
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= m+(ntest-7)*(npest-7).
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is a spherical
c            interpolating spline (fp=0).
c   ier=-2 : normal return. the spline returned is the weighted least-
c            squares constrained polynomial . in this extreme case
c            fp gives the upper bound for the smoothing factor s.
c   ier<-2 : warning. the coefficients of the spline returned have been
c            computed as the minimal norm least-squares solution of a
c            (numerically) rank deficient system. (-ier) gives the rank.
c            especially if the rank deficiency which can be computed as
c            6+(nt-8)*(np-7)+ier, is large the results may be inaccurate
c            they could also seriously depend on the value of eps.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters ntest and
c            npest.
c            probably causes : ntest or npest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the weighted least-squares
c            spherical spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing spline with
c            fp = s. probably causes : s too small or badly chosen eps.
c            there is an approximation returned but the corresponding
c            weighted sum of squared residuals does not satisfy the
c            condition abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing spline
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            weighted sum of squared residuals does not satisfy the
c            condition abs(fp-s)/s < tol.
c   ier=4  : error. no more knots can be added because the dimension
c            of the spherical spline 6+(nt-8)*(np-7) already exceeds
c            the number of data points m.
c            probably causes : either s or m too small.
c            the approximation returned is the weighted least-squares
c            spherical spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=5  : error. no more knots can be added because the additional
c            knot would (quasi) coincide with an old one.
c            probably causes : s too small or too large a weight to an
c            inaccurate data point.
c            the approximation returned is the weighted least-squares
c            spherical spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt<=1,  m>=2, ntest>=8 ,npest >=8, 0<eps<1,
c            0<=teta(i)<=pi, 0<=phi(i)<=2*pi, w(i)>0, i=1,...,m
c            lwrk1 >= 185+52*v+10*u+14*u*v+8*(u-1)*v**2+8*m
c            kwrk >= m+(ntest-7)*(npest-7)
c            if iopt=-1: 8<=nt<=ntest , 9<=np<=npest
c                        0<tt(5)<tt(6)<...<tt(nt-4)<pi
c                        0<tp(5)<tp(6)<...<tp(np-4)<2*pi
c            if iopt>=0: s>=0
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
c            space for computing the minimal least-squares solution of
c            a rank deficient system of linear equations. ier gives the
c            requested value for lwrk2. there is no approximation re-
c            turned but, having saved the information contained in nt,
c            np,tt,tp,wrk1, and having adjusted the value of lwrk2 and
c            the dimension of the array wrk2 accordingly, the user can
c            continue at the point the program was left, by calling
c            sphere with iopt=1.
c
c further comments:
c  by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the constrained weighted least-squares polynomial if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   r(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in r(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to choose s very small is strongly discouraged. this considerably
c   increases computation time and memory requirements. it may also
c   cause rank-deficiency (ier<-2) and endager numerical stability.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if sphere is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   sphere once more with the selected value for s but now with iopt=0.
c   indeed, sphere may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds ntest and
c   npest. indeed, if at a certain stage in sphere the number of knots
c   in one direction (say nt) has reached the value of its upper bound
c   (ntest), then from that moment on all subsequent knots are added
c   in the other (phi) direction. this may indicate that the value of
c   ntest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting ntest=8 (the lowest allowable value for
c   ntest), the user can indicate that he wants an approximation which
c   is a cubic polynomial in the variable teta.
c
c  other subroutines required:
c    fpback,fpbspl,fpsphe,fpdisc,fpgivs,fprank,fprati,fprota,fporde,
c    fprpsp
c
c  references:
c   dierckx p. : algorithms for smoothing data on the sphere with tensor
c                product splines, computing 32 (1984) 319-342.
c   dierckx p. : algorithms for smoothing data on the sphere with tensor
c                product splines, report tw62, dept. computer science,
c                k.u.leuven, 1983.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : july 1983
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      real*8 s,eps,fp
      integer iopt,m,ntest,npest,nt,np,lwrk1,lwrk2,kwrk,ier
c  ..array arguments..
      real*8 teta(m),phi(m),r(m),w(m),tt(ntest),tp(npest),
     * c((ntest-4)*(npest-4)),wrk1(lwrk1),wrk2(lwrk2)
      integer iwrk(kwrk)
c  ..local scalars..
      real*8 tol,pi,pi2,one
      integer i,ib1,ib3,ki,kn,kwest,la,lbt,lcc,lcs,lro,j,
     * lbp,lco,lf,lff,lfp,lh,lq,lst,lsp,lwest,maxit,ncest,ncc,ntt,
     * npp,nreg,nrint,ncof,nt4,np4
c  ..function references..
      real*8 atan
c  ..subroutine references..
c    fpsphe
c  ..
c  set constants
      one = 0.1e+01
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid,control is immediately repassed to the calling program.
      ier = 10
      if(eps.le.0. .or. eps.ge.1.) go to 80
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 80
      if(m.lt.2) go to 80
      if(ntest.lt.8 .or. npest.lt.8) go to 80
      nt4 = ntest-4
      np4 = npest-4
      ncest = nt4*np4
      ntt = ntest-7
      npp = npest-7
      ncc = 6+npp*(ntt-1)
      nrint = ntt+npp
      nreg = ntt*npp
      ncof = 6+3*npp
      ib1 = 4*npp
      ib3 = ib1+3
      if(ncof.gt.ib1) ib1 = ncof
      if(ncof.gt.ib3) ib3 = ncof
      lwest = 185+52*npp+10*ntt+14*ntt*npp+8*(m+(ntt-1)*npp**2)
      kwest = m+nreg
      if(lwrk1.lt.lwest .or. kwrk.lt.kwest) go to 80
      if(iopt.gt.0) go to 60
      pi = atan(one)*4
      pi2 = pi+pi
      do 20 i=1,m
        if(w(i).le.0.) go to 80
        if(teta(i).lt.0. .or. teta(i).gt.pi) go to 80
        if(phi(i) .lt.0. .or. phi(i).gt.pi2) go to 80
  20  continue
      if(iopt.eq.0) go to 60
      ntt = nt-8
      if(ntt.lt.0 .or. nt.gt.ntest) go to 80
      if(ntt.eq.0) go to 40
      tt(4) = 0.
      do 30 i=1,ntt
         j = i+4
         if(tt(j).le.tt(j-1) .or. tt(j).ge.pi) go to 80
  30  continue
  40  npp = np-8
      if(npp.lt.1 .or. np.gt.npest) go to 80
      tp(4) = 0.
      do 50 i=1,npp
         j = i+4
         if(tp(j).le.tp(j-1) .or. tp(j).ge.pi2) go to 80
  50  continue
      go to 70
  60  if(s.lt.0.) go to 80
  70  ier = 0
c  we partition the working space and determine the spline approximation
      kn = 1
      ki = kn+m
      lq = 2
      la = lq+ncc*ib3
      lf = la+ncc*ib1
      lff = lf+ncc
      lfp = lff+ncest
      lco = lfp+nrint
      lh = lco+nrint
      lbt = lh+ib3
      lbp = lbt+5*ntest
      lro = lbp+5*npest
      lcc = lro+npest
      lcs = lcc+npest
      lst = lcs+npest
      lsp = lst+m*4
      call fpsphe(iopt,m,teta,phi,r,w,s,ntest,npest,eps,tol,maxit,
     * ib1,ib3,ncest,ncc,nrint,nreg,nt,tt,np,tp,c,fp,wrk1(1),wrk1(lfp),
     * wrk1(lco),wrk1(lf),wrk1(lff),wrk1(lro),wrk1(lcc),wrk1(lcs),
     * wrk1(la),wrk1(lq),wrk1(lbt),wrk1(lbp),wrk1(lst),wrk1(lsp),
     * wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
  80  return
      end

