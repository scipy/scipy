      recursive subroutine parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,
     * nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier)
      implicit none
c  given the ordered set of m points x(i) in the idim-dimensional space
c  and given also a corresponding set of strictly increasing values u(i)
c  and the set of positive numbers w(i),i=1,2,...,m, subroutine parcur
c  determines a smooth approximating spline curve s(u), i.e.
c      x1 = s1(u)
c      x2 = s2(u)       ub <= u <= ue
c      .........
c      xidim = sidim(u)
c  with sj(u),j=1,2,...,idim spline functions of degree k with common
c  knots t(j),j=1,2,...,n.
c  if ipar=1 the values ub,ue and u(i),i=1,2,...,m must be supplied by
c  the user. if ipar=0 these values are chosen automatically by parcur
c  as  v(1) = 0
c      v(i) = v(i-1) + dist(x(i),x(i-1)) ,i=2,3,...,m
c      u(i) = v(i)/v(m) ,i=1,2,...,m
c      ub = u(1) = 0, ue = u(m) = 1.
c  if iopt=-1 parcur calculates the weighted least-squares spline curve
c  according to a given set of knots.
c  if iopt>=0 the number of knots of the splines sj(u) and the position
c  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
c  ness of s(u) is then achieved by minimalizing the discontinuity
c  jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
c  n-k-1. the amount of smoothness is determined by the condition that
c  f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
c  negative constant, called the smoothing factor.
c  the fit s(u) is given in the b-spline representation and can be
c  evaluated by means of subroutine curev.
c
c  calling sequence:
c     call parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,nc,c,
c    * fp,wrk,lwrk,iwrk,ier)
c
c  parameters:
c   iopt  : integer flag. on entry iopt must specify whether a weighted
c           least-squares spline curve (iopt=-1) or a smoothing spline
c           curve (iopt=0 or 1) must be determined.if iopt=0 the routine
c           will start with an initial set of knots t(i)=ub,t(i+k+1)=ue,
c           i=1,2,...,k+1. if iopt=1 the routine will continue with the
c           knots found at the last call of the routine.
c           attention: a call with iopt=1 must always be immediately
c           preceded by another call with iopt=1 or iopt=0.
c           unchanged on exit.
c   ipar  : integer flag. on entry ipar must specify whether (ipar=1)
c           the user will supply the parameter values u(i),ub and ue
c           or whether (ipar=0) these values are to be calculated by
c           parcur. unchanged on exit.
c   idim  : integer. on entry idim must specify the dimension of the
c           curve. 0 < idim < 11.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > k. unchanged on exit.
c   u     : real array of dimension at least (m). in case ipar=1,before
c           entry, u(i) must be set to the i-th value of the parameter
c           variable u for i=1,2,...,m. these values must then be
c           supplied in strictly ascending order and will be unchanged
c           on exit. in case ipar=0, on exit,array u will contain the
c           values u(i) as determined by parcur.
c   mx    : integer. on entry mx must specify the actual dimension of
c           the array x as declared in the calling (sub)program. mx must
c           not be too small (see x). unchanged on exit.
c   x     : real array of dimension at least idim*m.
c           before entry, x(idim*(i-1)+j) must contain the j-th coord-
c           inate of the i-th data point for i=1,2,...,m and j=1,2,...,
c           idim. unchanged on exit.
c   w     : real array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. unchanged on exit.
c           see also further comments.
c   ub,ue : real values. on entry (in case ipar=1) ub and ue must
c           contain the lower and upper bound for the parameter u.
c           ub <=u(1), ue>= u(m). if ipar = 0 these values will
c           automatically be set to 0 and 1 by parcur.
c   k     : integer. on entry k must specify the degree of the splines.
c           1<=k<=5. it is recommended to use cubic splines (k=3).
c           the user is strongly dissuaded from choosing k even,together
c           with a small s-value. unchanged on exit.
c   s     : real.on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the splines returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is nest=m+k+1, the number of knots
c           needed for interpolation (s=0). unchanged on exit.
c   n     : integer.
c           unless ier = 10 (in case iopt >=0), n will contain the
c           total number of knots of the smoothing spline curve returned
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : real array of dimension at least (nest).
c           on successful exit, this array will contain the knots of the
c           spline curve,i.e. the position of the interior knots t(k+2),
c           t(k+3),..,t(n-k-1) as well as the position of the additional
c           t(1)=t(2)=...=t(k+1)=ub and t(n-k)=...=t(n)=ue needed for
c           the b-spline representation.
c           if the computation mode iopt=1 is used, the values of t(1),
c           t(2),...,t(n) should be left unchanged between subsequent
c           calls. if the computation mode iopt=-1 is used, the values
c           t(k+2),...,t(n-k-1) must be supplied by the user, before
c           entry. see also the restrictions (ier=10).
c   nc    : integer. on entry nc must specify the actual dimension of
c           the array c as declared in the calling (sub)program. nc
c           must not be too small (see c). unchanged on exit.
c   c     : real array of dimension at least (nest*idim).
c           on successful exit, this array will contain the coefficients
c           in the b-spline representation of the spline curve s(u),i.e.
c           the b-spline coefficients of the spline sj(u) will be given
c           in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
c   fp    : real. unless ier = 10, fp contains the weighted sum of
c           squared residuals of the spline curve returned.
c   wrk   : real array of dimension at least m*(k+1)+nest*(6+idim+3*k).
c           used as working space. if the computation mode iopt=1 is
c           used, the values wrk(1),...,wrk(n) should be left unchanged
c           between subsequent calls.
c   lwrk  : integer. on entry,lwrk must specify the actual dimension of
c           the array wrk as declared in the calling (sub)program. lwrk
c           must not be too small (see wrk). unchanged on exit.
c   iwrk  : integer array of dimension at least (nest).
c           used as working space. if the computation mode iopt=1 is
c           used,the values iwrk(1),...,iwrk(n) should be left unchanged
c           between subsequent calls.
c   ier   : integer. unless the routine detects an error, ier contains a
c           non-positive value on exit, i.e.
c    ier=0  : normal return. the curve returned has a residual sum of
c             squares fp such that abs(fp-s)/s <= tol with tol a relat-
c             ive tolerance set to 0.001 by the program.
c    ier=-1 : normal return. the curve returned is an interpolating
c             spline curve (fp=0).
c    ier=-2 : normal return. the curve returned is the weighted least-
c             squares polynomial curve of degree k.in this extreme case
c             fp gives the upper bound fp0 for the smoothing factor s.
c    ier=1  : error. the required storage space exceeds the available
c             storage space, as specified by the parameter nest.
c             probably causes : nest too small. if nest is already
c             large (say nest > m/2), it may also indicate that s is
c             too small
c             the approximation returned is the least-squares spline
c             curve according to the knots t(1),t(2),...,t(n). (n=nest)
c             the parameter fp gives the corresponding weighted sum of
c             squared residuals (fp>s).
c    ier=2  : error. a theoretically impossible result was found during
c             the iteration process for finding a smoothing spline curve
c             with fp = s. probably causes : s too small.
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=3  : error. the maximal number of iterations maxit (set to 20
c             by the program) allowed for finding a smoothing curve
c             with fp=s has been reached. probably causes : s too small
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=10 : error. on entry, the input data are controlled on validity
c             the following restrictions must be satisfied.
c             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
c             0<=ipar<=1, 0<idim<=10, lwrk>=(k+1)*m+nest*(6+idim+3*k),
c             nc>=nest*idim
c             if ipar=0: sum j=1,idim (x(idim*i+j)-x(idim*(i-1)+j))**2>0
c                        i=1,2,...,m-1.
c             if ipar=1: ub<=u(1)<u(2)<...<u(m)<=ue
c             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
c                         ub<t(k+2)<t(k+3)<...<t(n-k-1)<ue
c                            (ub=0 and ue=1 in case ipar=0)
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points uu(j) such that
c                         t(j) < uu(j) < t(j+k+1), j=1,2,...,n-k-1
c             if iopt>=0: s>=0
c                         if s=0 : nest >= m+k+1
c             if one of these conditions is found to be violated,control
c             is immediately repassed to the calling program. in that
c             case there is no approximation returned.
c
c  further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the curve will be too smooth and signal will be
c   lost ; if s is too small the curve will pick up too much noise. in
c   the extreme cases the program will return an interpolating curve if
c   s=0 and the least-squares polynomial curve of degree k if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   x(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in x(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial curve and the upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximating curve shows more detail) to obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if parcur is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   curve underlying the data. but, if the computation mode iopt=1 is
c   used, the knots returned may also depend on the s-values at previous
c   calls (if these were smaller). therefore, if after a number of
c   trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   parcur once more with the selected value for s but now with iopt=0.
c   indeed, parcur may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c   the form of the approximating curve can strongly be affected by
c   the choice of the parameter values u(i). if there is no physical
c   reason for choosing a particular parameter u, often good results
c   will be obtained with the choice of parcur (in case ipar=0), i.e.
c        v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m
c   where
c        q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 )
c   other possibilities for q(i) are
c        q(i)= sum j=1,idim (xj(i)-xj(i-1))**2
c        q(i)= sum j=1,idim abs(xj(i)-xj(i-1))
c        q(i)= max j=1,idim abs(xj(i)-xj(i-1))
c        q(i)= 1
c
c  other subroutines required:
c    fpback,fpbspl,fpchec,fppara,fpdisc,fpgivs,fpknot,fprati,fprota
c
c  references:
c   dierckx p. : algorithms for smoothing data with periodic and
c                parametric splines, computer graphics and image
c                processing 20 (1982) 171-184.
c   dierckx p. : algorithms for smoothing data with periodic and param-
c                etric splines, report tw55, dept. computer science,
c                k.u.leuven, 1981.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      real*8 ub,ue,s,fp
      integer iopt,ipar,idim,m,mx,k,nest,n,nc,lwrk,ier
c  ..array arguments..
      real*8 u(m),x(mx),w(m),t(nest),c(nc),wrk(lwrk)
      integer iwrk(nest)
c  ..local scalars..
      real*8 tol,dist
      integer i,ia,ib,ifp,ig,iq,iz,i1,i2,j,k1,k2,lwest,maxit,nmin,ncc
c ..function references
      real*8 sqrt
c  ..
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 90
      if(ipar.lt.0 .or. ipar.gt.1) go to 90
      if(idim.le.0 .or. idim.gt.10) go to 90
      if(k.le.0 .or. k.gt.5) go to 90
      k1 = k+1
      k2 = k1+1
      nmin = 2*k1
      if(m.lt.k1 .or. nest.lt.nmin) go to 90
      ncc = nest*idim
      if(mx.lt.m*idim .or. nc.lt.ncc) go to 90
      lwest = m*k1+nest*(6+idim+3*k)
      if(lwrk.lt.lwest) go to 90
      if(ipar.ne.0 .or. iopt.gt.0) go to 40
      i1 = 0
      i2 = idim
      u(1) = 0.
      do 20 i=2,m
         dist = 0.
         do 10 j=1,idim
            i1 = i1+1
            i2 = i2+1
            dist = dist+(x(i2)-x(i1))**2
  10     continue
         u(i) = u(i-1)+sqrt(dist)
  20  continue
      if(u(m).le.0.) go to 90
      do 30 i=2,m
         u(i) = u(i)/u(m)
  30  continue
      ub = 0.
      ue = 1.
      u(m) = ue
  40  if(ub.gt.u(1) .or. ue.lt.u(m) .or. w(1).le.0.) go to 90
      do 50 i=2,m
         if(u(i-1).ge.u(i) .or. w(i).le.0.) go to 90
  50  continue
      if(iopt.ge.0) go to 70
      if(n.lt.nmin .or. n.gt.nest) go to 90
      j = n
      do 60 i=1,k1
         t(i) = ub
         t(j) = ue
         j = j-1
  60  continue
      call fpchec(u,m,t,n,k,ier)
      if (ier.eq.0) go to 80
      go to 90
  70  if(s.lt.0.) go to 90
      if(s.eq.0. .and. nest.lt.(m+k1)) go to 90
      ier = 0
c we partition the working space and determine the spline curve.
  80  ifp = 1
      iz = ifp+nest
      ia = iz+ncc
      ib = ia+nest*k1
      ig = ib+nest*k2
      iq = ig+nest*k2
      call fppara(iopt,idim,m,u,mx,x,w,ub,ue,k,s,nest,tol,maxit,k1,k2,
     * n,t,ncc,c,fp,wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),
     * iwrk,ier)
  90  return
      end
