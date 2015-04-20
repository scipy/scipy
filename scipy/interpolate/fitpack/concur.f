      subroutine concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb,ie,de,ne,k,s,
     * nest,n,t,nc,c,np,cp,fp,wrk,lwrk,iwrk,ier)
      implicit none
c  given the ordered set of m points x(i) in the idim-dimensional space
c  and given also a corresponding set of strictly increasing values u(i)
c  and the set of positive numbers w(i),i=1,2,...,m, subroutine concur
c  determines a smooth approximating spline curve s(u), i.e.
c      x1 = s1(u)
c      x2 = s2(u)      ub = u(1) <= u <= u(m) = ue
c      .........
c      xidim = sidim(u)
c  with sj(u),j=1,2,...,idim spline functions of odd degree k with
c  common knots t(j),j=1,2,...,n.
c  in addition these splines will satisfy the following boundary
c  constraints        (l)
c      if ib > 0 :  sj   (u(1)) = db(idim*l+j) ,l=0,1,...,ib-1
c  and                (l)
c      if ie > 0 :  sj   (u(m)) = de(idim*l+j) ,l=0,1,...,ie-1.
c  if iopt=-1 concur calculates the weighted least-squares spline curve
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
c     call concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb,ie,de,ne,k,s,nest,n,
c    * t,nc,c,np,cp,fp,wrk,lwrk,iwrk,ier)
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
c   idim  : integer. on entry idim must specify the dimension of the
c           curve. 0 < idim < 11.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > k-max(ib-1,0)-max(ie-1,0). unchanged on exit.
c   u     : real array of dimension at least (m). before entry,
c           u(i) must be set to the i-th value of the parameter variable
c           u for i=1,2,...,m. these values must be supplied in
c           strictly ascending order and will be unchanged on exit.
c   mx    : integer. on entry mx must specify the actual dimension of
c           the arrays x and xx as declared in the calling (sub)program
c           mx must not be too small (see x). unchanged on exit.
c   x     : real array of dimension at least idim*m.
c           before entry, x(idim*(i-1)+j) must contain the j-th coord-
c           inate of the i-th data point for i=1,2,...,m and j=1,2,...,
c           idim. unchanged on exit.
c   xx    : real array of dimension at least idim*m.
c           used as working space. on exit xx contains the coordinates
c           of the data points to which a spline curve with zero deriv-
c           ative constraints has been determined.
c           if the computation mode iopt =1 is used xx should be left
c           unchanged between calls.
c   w     : real array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. unchanged on exit.
c           see also further comments.
c   ib    : integer. on entry ib must specify the number of derivative
c           constraints for the curve at the begin point. 0<=ib<=(k+1)/2
c           unchanged on exit.
c   db    : real array of dimension nb. before entry db(idim*l+j) must
c           contain the l-th order derivative of sj(u) at u=u(1) for
c           j=1,2,...,idim and l=0,1,...,ib-1 (if ib>0).
c           unchanged on exit.
c   nb    : integer, specifying the dimension of db. nb>=max(1,idim*ib)
c           unchanged on exit.
c   ie    : integer. on entry ie must specify the number of derivative
c           constraints for the curve at the end point. 0<=ie<=(k+1)/2
c           unchanged on exit.
c   de    : real array of dimension ne. before entry de(idim*l+j) must
c           contain the l-th order derivative of sj(u) at u=u(m) for
c           j=1,2,...,idim and l=0,1,...,ie-1 (if ie>0).
c           unchanged on exit.
c   ne    : integer, specifying the dimension of de. ne>=max(1,idim*ie)
c           unchanged on exit.
c   k     : integer. on entry k must specify the degree of the splines.
c           k=1,3 or 5.
c           unchanged on exit.
c   s     : real.on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the splines returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is nest=m+k+1+max(0,ib-1)+max(0,ie-1),
c           the number of knots needed for interpolation (s=0).
c           unchanged on exit.
c   n     : integer.
c           unless ier = 10 (in case iopt >=0), n will contain the
c           total number of knots of the smoothing spline curve returned
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the knots of the
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
c           on succesful exit, this array will contain the coefficients
c           in the b-spline representation of the spline curve s(u),i.e.
c           the b-spline coefficients of the spline sj(u) will be given
c           in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
c   cp    : real array of dimension at least 2*(k+1)*idim.
c           on exit cp will contain the b-spline coefficients of a
c           polynomial curve which satisfies the boundary constraints.
c           if the computation mode iopt =1 is used cp should be left
c           unchanged between calls.
c   np    : integer. on entry np must specify the actual dimension of
c           the array cp as declared in the calling (sub)program. np
c           must not be too small (see cp). unchanged on exit.
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
c             spline curve, satisfying the constraints (fp=0).
c    ier=-2 : normal return. the curve returned is the weighted least-
c             squares polynomial curve of degree k, satisfying the
c             constraints. in this extreme case fp gives the upper
c             bound fp0 for the smoothing factor s.
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
c             the iteration proces for finding a smoothing spline curve
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
c             -1<=iopt<=1, k = 1,3 or 5, m>k-max(0,ib-1)-max(0,ie-1),
c             nest>=2k+2, 0<idim<=10, lwrk>=(k+1)*m+nest*(6+idim+3*k),
c             nc >=nest*idim ,u(1)<u(2)<...<u(m),w(i)>0 i=1,2,...,m,
c             mx>=idim*m,0<=ib<=(k+1)/2,0<=ie<=(k+1)/2,nb>=1,ne>=1,
c             nb>=ib*idim,ne>=ib*idim,np>=2*(k+1)*idim,
c             if iopt=-1:2*k+2<=n<=min(nest,mmax) with mmax = m+k+1+
c                        max(0,ib-1)+max(0,ie-1)
c                        u(1)<t(k+2)<t(k+3)<...<t(n-k-1)<u(m)
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points uu(j) such that
c                         t(j) < uu(j) < t(j+k+1), j=1+max(0,ib-1),...
c                                                   ,n+k-1-max(0,ie-1)
c             if iopt>=0: s>=0
c                         if s=0 : nest >=mmax (see above)
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
c   time if concur is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   curve underlying the data. but, if the computation mode iopt=1 is
c   used, the knots returned may also depend on the s-values at previous
c   calls (if these were smaller). therefore, if after a number of
c   trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   concur once more with the selected value for s but now with iopt=0.
c   indeed, concur may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c   the form of the approximating curve can strongly be affected by
c   the choice of the parameter values u(i). if there is no physical
c   reason for choosing a particular parameter u, often good results
c   will be obtained with the choice
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
c    fpback,fpbspl,fpched,fpcons,fpdisc,fpgivs,fpknot,fprati,fprota
c    curev,fppocu,fpadpo,fpinst
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
      real*8 s,fp
      integer iopt,idim,m,mx,ib,nb,ie,ne,k,nest,n,nc,np,lwrk,ier
c  ..array arguments..
      real*8 u(m),x(mx),xx(mx),db(nb),de(ne),w(m),t(nest),c(nc),wrk(lwrk
     *)
      real*8 cp(np)
      integer iwrk(nest)
c  ..local scalars..
      real*8 tol
      integer i,ib1,ie1,ja,jb,jfp,jg,jq,jz,j,k1,k2,lwest,maxit,nmin,
     * ncc,kk,mmin,nmax,mxx
c ..function references
      integer max0
c  ..
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 90
      if(idim.le.0 .or. idim.gt.10) go to 90
      if(k.le.0 .or. k.gt.5) go to 90
      k1 = k+1
      kk = k1/2
      if(kk*2.ne.k1) go to 90
      k2 = k1+1
      if(ib.lt.0 .or. ib.gt.kk) go to 90
      if(ie.lt.0 .or. ie.gt.kk) go to 90
      nmin = 2*k1
      ib1 = max0(0,ib-1)
      ie1 = max0(0,ie-1)
      mmin = k1-ib1-ie1
      if(m.lt.mmin .or. nest.lt.nmin) go to 90
      if(nb.lt.(idim*ib) .or. ne.lt.(idim*ie)) go to 90
      if(np.lt.(2*k1*idim)) go to 90
      mxx = m*idim
      ncc = nest*idim
      if(mx.lt.mxx .or. nc.lt.ncc) go to 90
      lwest = m*k1+nest*(6+idim+3*k)
      if(lwrk.lt.lwest) go to 90
      if(w(1).le.0.) go to 90
      do 10 i=2,m
         if(u(i-1).ge.u(i) .or. w(i).le.0.) go to 90
  10  continue
      if(iopt.ge.0) go to 30
      if(n.lt.nmin .or. n.gt.nest) go to 90
      j = n
      do 20 i=1,k1
         t(i) = u(1)
         t(j) = u(m)
         j = j-1
  20  continue
      call fpched(u,m,t,n,k,ib,ie,ier)
      if (ier.eq.0) go to 40
      go to 90
  30  if(s.lt.0.) go to 90
      nmax = m+k1+ib1+ie1
      if(s.eq.0. .and. nest.lt.nmax) go to 90
      ier = 0
      if(iopt.gt.0) go to 70
c  we determine a polynomial curve satisfying the boundary constraints.
  40  call fppocu(idim,k,u(1),u(m),ib,db,nb,ie,de,ne,cp,np)
c  we generate new data points which will be approximated by a spline
c  with zero derivative constraints.
      j = nmin
      do 50 i=1,k1
        wrk(i) = u(1)
        wrk(j) = u(m)
        j = j-1
  50  continue
c  evaluate the polynomial curve
      call curev(idim,wrk,nmin,cp,np,k,u,m,xx,mxx,ier)
c  substract from the old data, the values of the polynomial curve
      do 60 i=1,mxx
        xx(i) = x(i)-xx(i)
  60  continue
c we partition the working space and determine the spline curve.
  70  jfp = 1
      jz = jfp+nest
      ja = jz+ncc
      jb = ja+nest*k1
      jg = jb+nest*k2
      jq = jg+nest*k2
      call fpcons(iopt,idim,m,u,mxx,xx,w,ib,ie,k,s,nest,tol,maxit,k1,
     * k2,n,t,ncc,c,fp,wrk(jfp),wrk(jz),wrk(ja),wrk(jb),wrk(jg),wrk(jq),
     *
     * iwrk,ier)
c  add the polynomial curve to the calculated spline.
      call fpadpo(idim,t,n,c,ncc,k,cp,np,wrk(jz),wrk(ja),wrk(jb))
  90  return
      end
