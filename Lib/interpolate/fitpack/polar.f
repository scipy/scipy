      subroutine polar(iopt,m,x,y,z,w,rad,s,nuest,nvest,eps,nu,tu,
     *  nv,tv,u,v,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c  subroutine polar fits a smooth function f(x,y) to a set of data
c  points (x(i),y(i),z(i)) scattered arbitrarily over an approximation
c  domain  x**2+y**2 <= rad(atan(y/x))**2. through the transformation
c    x = u*rad(v)*cos(v) , y = u*rad(v)*sin(v)
c  the approximation problem is reduced to the determination of a bi-
c  cubic spline s(u,v) fitting a corresponding set of data points
c  (u(i),v(i),z(i)) on the rectangle 0<=u<=1,-pi<=v<=pi.
c  in order to have continuous partial derivatives
c              i+j
c             d   f(0,0)
c    g(i,j) = ----------
c                i   j
c              dx  dy
c
c  s(u,v)=f(x,y) must satisfy the following conditions
c
c    (1) s(0,v) = g(0,0)   -pi <=v<= pi.
c
c        d s(0,v)
c    (2) -------- = rad(v)*(cos(v)*g(1,0)+sin(v)*g(0,1))
c        d u
c                                                    -pi <=v<= pi
c         2
c        d s(0,v)         2       2             2
c    (3) -------- = rad(v)*(cos(v)*g(2,0)+sin(v)*g(0,2)+sin(2*v)*g(1,1))
c           2
c        d u                                         -pi <=v<= pi
c
c  moreover, s(u,v) must be periodic in the variable v, i.e.
c
c         j            j
c        d s(u,-pi)   d s(u,pi)
c    (4) ---------- = ---------   0 <=u<= 1, j=0,1,2
c           j           j
c        d v         d v
c
c  if iopt(1) < 0 circle calculates a weighted least-squares spline
c  according to a given set of knots in u- and v- direction.
c  if iopt(1) >=0, the number of knots in each direction and their pos-
c  ition tu(j),j=1,2,...,nu ; tv(j),j=1,2,...,nv are chosen automatical-
c  ly by the routine. the smoothness of s(u,v) is then achieved by mini-
c  malizing the discontinuity jumps of the derivatives of the spline
c  at the knots. the amount of smoothness of s(u,v) is determined  by
c  the condition that fp = sum((w(i)*(z(i)-s(u(i),v(i))))**2) be <= s,
c  with s a given non-negative constant.
c  the bicubic spline is given in its standard b-spline representation
c  and the corresponding function f(x,y) can be evaluated by means of
c  function program evapol.
c
c calling sequence:
c     call polar(iopt,m,x,y,z,w,rad,s,nuest,nvest,eps,nu,tu,
c    *  nv,tv,u,v,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer array of dimension 3, specifying different options.
c          unchanged on exit.
c  iopt(1):on entry iopt(1) must specify whether a weighted
c          least-squares polar spline (iopt(1)=-1) or a smoothing
c          polar spline (iopt(1)=0 or 1) must be determined.
c          if iopt(1)=0 the routine will start with an initial set of
c          knots tu(i)=0,tu(i+4)=1,i=1,...,4;tv(i)=(2*i-9)*pi,i=1,...,8.
c          if iopt(1)=1 the routine will continue with the set of knots
c          found at the last call of the routine.
c          attention: a call with iopt(1)=1 must always be immediately
c          preceded by another call with iopt(1) = 1 or iopt(1) = 0.
c  iopt(2):on entry iopt(2) must specify the requested order of conti-
c          nuity for f(x,y) at the origin.
c          if iopt(2)=0 only condition (1) must be fulfilled,
c          if iopt(2)=1 conditions (1)+(2) must be fulfilled and
c          if iopt(2)=2 conditions (1)+(2)+(3) must be fulfilled.
c  iopt(3):on entry iopt(3) must specify whether (iopt(3)=1) or not
c          (iopt(3)=0) the approximation f(x,y) must vanish at the
c          boundary of the approximation domain.
c  m     : integer. on entry m must specify the number of data points.
c          m >= 4-iopt(2)-iopt(3) unchanged on exit.
c  x     : real array of dimension at least (m).
c  y     : real array of dimension at least (m).
c  z     : real array of dimension at least (m).
c          before entry, x(i),y(i),z(i) must be set to the co-ordinates
c          of the i-th data point, for i=1,...,m. the order of the data
c          points is immaterial. unchanged on exit.
c  w     : real array of dimension at least (m). before entry, w(i) must
c          be set to the i-th value in the set of weights. the w(i) must
c          be strictly positive. unchanged on exit.
c  rad   : real function subprogram defining the boundary of the approx-
c          imation domain, i.e   x = rad(v)*cos(v) , y = rad(v)*sin(v),
c          -pi <= v <= pi.
c          must be declared external in the calling (sub)program.
c  s     : real. on entry (in case iopt(1) >=0) s must specify the
c          smoothing factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nuest : integer. unchanged on exit.
c  nvest : integer. unchanged on exit.
c          on entry, nuest and nvest must specify an upper bound for the
c          number of knots required in the u- and v-directions resp.
c          these numbers will also determine the storage space needed by
c          the routine. nuest >= 8, nvest >= 8.
c          in most practical situation nuest = nvest = 8+sqrt(m/2) will
c          be sufficient. see also further comments.
c  eps   : real.
c          on entry, eps must specify a threshold for determining the
c          effective rank of an over-determined linear system of equat-
c          ions. 0 < eps < 1.  if the number of decimal digits in the
c          computer representation of a real number is q, then 10**(-q)
c          is a suitable value for eps in most practical applications.
c          unchanged on exit.
c  nu    : integer.
c          unless ier=10 (in case iopt(1) >=0),nu will contain the total
c          number of knots with respect to the u-variable, of the spline
c          approximation returned. if the computation mode iopt(1)=1
c          is used, the value of nu should be left unchanged between
c          subsequent calls.
c          in case iopt(1)=-1,the value of nu must be specified on entry
c  tu    : real array of dimension at least nuest.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the u-variable, i.e. the position
c          of the interior knots tu(5),...,tu(nu-4) as well as the
c          position of the additional knots tu(1)=...=tu(4)=0 and
c          tu(nu-3)=...=tu(nu)=1 needed for the b-spline representation
c          if the computation mode iopt(1)=1 is used,the values of
c          tu(1),...,tu(nu) should be left unchanged between subsequent
c          calls. if the computation mode iopt(1)=-1 is used,the values
c          tu(5),...tu(nu-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  nv    : integer.
c          unless ier=10 (in case iopt(1)>=0), nv will contain the total
c          number of knots with respect to the v-variable, of the spline
c          approximation returned. if the computation mode iopt(1)=1
c          is used, the value of nv should be left unchanged between
c          subsequent calls. in case iopt(1)=-1, the value of nv should
c          be specified on entry.
c  tv    : real array of dimension at least nvest.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the v-variable, i.e. the position of
c          the interior knots tv(5),...,tv(nv-4) as well as the position
c          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
c          tv(nv) needed for the b-spline representation.
c          if the computation mode iopt(1)=1 is used, the values of
c          tv(1),...,tv(nv) should be left unchanged between subsequent
c          calls. if the computation mode iopt(1)=-1 is used,the values
c          tv(5),...tv(nv-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  u     : real array of dimension at least (m).
c  v     : real array of dimension at least (m).
c          on succesful exit, u(i),v(i) contains the co-ordinates of
c          the i-th data point with respect to the transformed rectan-
c          gular approximation domain, for i=1,2,...,m.
c          if the computation mode iopt(1)=1 is used the values of
c          u(i),v(i) should be left unchanged between subsequent calls.
c  c     : real array of dimension at least (nuest-4)*(nvest-4).
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(u,v).
c  fp    : real. unless ier=10, fp contains the weighted sum of
c          squared residuals of the spline approximation returned.
c  wrk1  : real array of dimension (lwrk1). used as workspace.
c          if the computation mode iopt(1)=1 is used the value of
c          wrk1(1) should be left unchanged between subsequent calls.
c          on exit wrk1(2),wrk1(3),...,wrk1(1+ncof) will contain the
c          values d(i)/max(d(i)),i=1,...,ncof=1+iopt(2)*(iopt(2)+3)/2+
c          (nv-7)*(nu-5-iopt(2)-iopt(3)) with d(i) the i-th diagonal el-
c          ement of the triangular matrix for calculating the b-spline
c          coefficients.it includes those elements whose square is < eps
c          which are treated as 0 in the case of rank deficiency(ier=-2)
c  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
c          the array wrk1 as declared in the calling (sub)program.
c          lwrk1 must not be too small. let
c            k = nuest-7, l = nvest-7, p = 1+iopt(2)*(iopt(2)+3)/2,
c            q = k+2-iopt(2)-iopt(3) then
c          lwrk1 >= 129+10*k+21*l+k*l+(p+l*q)*(1+8*l+p)+8*m
c  wrk2  : real array of dimension (lwrk2). used as workspace, but
c          only in the case a rank deficient system is encountered.
c  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
c          the array wrk2 as declared in the calling (sub)program.
c          lwrk2 > 0 . a save upper bound  for lwrk2 = (p+l*q+1)*(4*l+p)
c          +p+l*q where p,l,q are as above. if there are enough data
c          points, scattered uniformly over the approximation domain
c          and if the smoothing factor s is not too small, there is a
c          good chance that this extra workspace is not needed. a lot
c          of memory might therefore be saved by setting lwrk2=1.
c          (see also ier > 10)
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= m+(nuest-7)*(nvest-7).
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is an interpolating
c            spline (fp=0).
c   ier=-2 : normal return. the spline returned is the weighted least-
c            squares constrained polynomial . in this extreme case
c            fp gives the upper bound for the smoothing factor s.
c   ier<-2 : warning. the coefficients of the spline returned have been
c            computed as the minimal norm least-squares solution of a
c            (numerically) rank deficient system. (-ier) gives the rank.
c            especially if the rank deficiency which can be computed as
c            1+iopt(2)*(iopt(2)+3)/2+(nv-7)*(nu-5-iopt(2)-iopt(3))+ier
c            is large the results may be inaccurate.
c            they could also seriously depend on the value of eps.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nuest and
c            nvest.
c            probably causes : nuest or nvest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the weighted least-squares
c            polar spline according to the current set of knots.
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
c            of the spline 1+iopt(2)*(iopt(2)+3)/2+(nv-7)*(nu-5-iopt(2)
c            -iopt(3)) already exceeds the number of data points m.
c            probably causes : either s or m too small.
c            the approximation returned is the weighted least-squares
c            polar spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=5  : error. no more knots can be added because the additional
c            knot would (quasi) coincide with an old one.
c            probably causes : s too small or too large a weight to an
c            inaccurate data point.
c            the approximation returned is the weighted least-squares
c            polar spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt(1)<=1 , 0<=iopt(2)<=2 , 0<=iopt(3)<=1 ,
c            m>=4-iopt(2)-iopt(3) , nuest>=8 ,nvest >=8, 0<eps<1,
c            0<=teta(i)<=pi, 0<=phi(i)<=2*pi, w(i)>0, i=1,...,m
c            lwrk1 >= 129+10*k+21*l+k*l+(p+l*q)*(1+8*l+p)+8*m
c            kwrk >= m+(nuest-7)*(nvest-7)
c            if iopt(1)=-1:9<=nu<=nuest,9+iopt(2)*(iopt(2)+1)<=nv<=nvest
c                          0<tu(5)<tu(6)<...<tu(nu-4)<1
c                          -pi<tv(5)<tv(6)<...<tv(nv-4)<pi
c            if iopt(1)>=0: s>=0
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
c            space for computing the minimal least-squares solution of
c            a rank deficient system of linear equations. ier gives the
c            requested value for lwrk2. there is no approximation re-
c            turned but, having saved the information contained in nu,
c            nv,tu,tv,wrk1,u,v and having adjusted the value of lwrk2
c            and the dimension of the array wrk2 accordingly, the user
c            can continue at the point the program was left, by calling
c            polar with iopt(1)=1.
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
c   z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in z(i)
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
c   must set iopt(1)=0.
c   if iopt(1)=1 the program will continue with the set of knots found
c   at the last call of the routine. this will save a lot of computation
c   time if polar is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt(1)=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt(1)=1,the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   polar once more with the selected value for s but now with iopt(1)=0
c   indeed, polar may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nuest and
c   nvest. indeed, if at a certain stage in polar the number of knots
c   in one direction (say nu) has reached the value of its upper bound
c   (nuest), then from that moment on all subsequent knots are added
c   in the other (v) direction. this may indicate that the value of
c   nuest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c
c  other subroutines required:
c    fpback,fpbspl,fppola,fpdisc,fpgivs,fprank,fprati,fprota,fporde,
c    fprppo
c
c  references:
c   dierckx p.: an algorithm for fitting data over a circle using tensor
c               product splines,j.comp.appl.maths 15 (1986) 161-173.
c   dierckx p.: an algorithm for fitting data on a circle using tensor
c               product splines, report tw68, dept. computer science,
c               k.u.leuven, 1984.
c   dierckx p.: curve and surface fitting with splines, monographs on
c               numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : june 1984
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      real*8 s,eps,fp
      integer m,nuest,nvest,nu,nv,lwrk1,lwrk2,kwrk,ier
c  ..array arguments..
      real*8 x(m),y(m),z(m),w(m),tu(nuest),tv(nvest),u(m),v(m),
     * c((nuest-4)*(nvest-4)),wrk1(lwrk1),wrk2(lwrk2)
      integer iopt(3),iwrk(kwrk)
c  ..user specified function
      real*8 rad
c  ..local scalars..
      real*8 tol,pi,dist,r,one
      integer i,ib1,ib3,ki,kn,kwest,la,lbu,lcc,lcs,lro,j
     * lbv,lco,lf,lff,lfp,lh,lq,lsu,lsv,lwest,maxit,ncest,ncc,nuu,
     * nvv,nreg,nrint,nu4,nv4,iopt1,iopt2,iopt3,ipar,nvmin
c  ..function references..
      real*8 datan2,sqrt
      external rad
c  ..subroutine references..
c    fppola
c  ..
c  set up constants
      one = 1d0
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid,control is immediately repassed to the calling program.
      ier = 10
      if(eps.le.0. .or. eps.ge.1.) go to 60
      iopt1 = iopt(1)
      if(iopt1.lt.(-1) .or. iopt1.gt.1) go to 60
      iopt2 = iopt(2)
      if(iopt2.lt.0 .or. iopt2.gt.2) go to 60
      iopt3 = iopt(3)
      if(iopt3.lt.0 .or. iopt3.gt.1) go to 60
      if(m.lt.(4-iopt2-iopt3)) go to 60
      if(nuest.lt.8 .or. nvest.lt.8) go to 60
      nu4 = nuest-4
      nv4 = nvest-4
      ncest = nu4*nv4
      nuu = nuest-7
      nvv = nvest-7
      ipar = 1+iopt2*(iopt2+3)/2
      ncc = ipar+nvv*(nuest-5-iopt2-iopt3)
      nrint = nuu+nvv
      nreg = nuu*nvv
      ib1 = 4*nvv
      ib3 = ib1+ipar
      lwest = ncc*(1+ib1+ib3)+2*nrint+ncest+m*8+ib3+5*nuest+12*nvest
      kwest = m+nreg
      if(lwrk1.lt.lwest .or. kwrk.lt.kwest) go to 60
      if(iopt1.gt.0) go to 40
      do 10 i=1,m
        if(w(i).le.0.) go to 60
        dist = x(i)**2+y(i)**2
        u(i) = 0.
        v(i) = 0.
        if(dist.le.0.) go to 10
        v(i) = datan2(y(i),x(i))
        r = rad(v(i))
        if(r.le.0.) go to 60
        u(i) = sqrt(dist)/r
        if(u(i).gt.one) go to 60
  10  continue
      if(iopt1.eq.0) go to 40
      nuu = nu-8
      if(nuu.lt.1 .or. nu.gt.nuest) go to 60
      tu(4) = 0.
      do 20 i=1,nuu
         j = i+4
         if(tu(j).le.tu(j-1) .or. tu(j).ge.one) go to 60
  20  continue
      nvv = nv-8
      nvmin = 9+iopt2*(iopt2+1)
      if(nv.lt.nvmin .or. nv.gt.nvest) go to 60
      pi = datan2(0d0,-one)
      tv(4) = -pi
      do 30 i=1,nvv
         j = i+4
         if(tv(j).le.tv(j-1) .or. tv(j).ge.pi) go to 60
  30  continue
      go to 50
  40  if(s.lt.0.) go to 60
  50  ier = 0
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
      lbu = lh+ib3
      lbv = lbu+5*nuest
      lro = lbv+5*nvest
      lcc = lro+nvest
      lcs = lcc+nvest
      lsu = lcs+nvest*5
      lsv = lsu+m*4
      call fppola(iopt1,iopt2,iopt3,m,u,v,z,w,rad,s,nuest,nvest,eps,tol,
     *
     * maxit,ib1,ib3,ncest,ncc,nrint,nreg,nu,tu,nv,tv,c,fp,wrk1(1),
     * wrk1(lfp),wrk1(lco),wrk1(lf),wrk1(lff),wrk1(lro),wrk1(lcc),
     * wrk1(lcs),wrk1(la),wrk1(lq),wrk1(lbu),wrk1(lbv),wrk1(lsu),
     * wrk1(lsv),wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
  60  return
      end

