      recursive subroutine concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,
     * n,t,c,sq,sx,bind,wrk,lwrk,iwrk,kwrk,ier)
      implicit none
c  given the set of data points (x(i),y(i)) and the set of positive
c  numbers w(i), i=1,2,...,m,subroutine concon determines a cubic spline
c  approximation s(x) which satisfies the following local convexity
c  constraints  s''(x(i))*v(i) <= 0, i=1,2,...,m.
c  the number of knots n and the position t(j),j=1,2,...n is chosen
c  automatically by the routine in a way that
c       sq = sum((w(i)*(y(i)-s(x(i))))**2) be <= s.
c  the fit is given in the b-spline representation (b-spline coef-
c  ficients c(j),j=1,2,...n-4) and can be evaluated by means of
c  subroutine splev.
c
c  calling sequence:
c
c     call concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,
c    * sx,bind,wrk,lwrk,iwrk,kwrk,ier)
c
c  parameters:
c    iopt: integer flag.
c          if iopt=0, the routine will start with the minimal number of
c          knots to guarantee that the convexity conditions will be
c          satisfied. if iopt=1, the routine will continue with the set
c          of knots found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately
c          preceded by another call with iopt=1 or iopt=0.
c          unchanged on exit.
c    m   : integer. on entry m must specify the number of data points.
c          m > 3. unchanged on exit.
c    x   : real array of dimension at least (m). before entry, x(i)
c          must be set to the i-th value of the independent variable x,
c          for i=1,2,...,m. these values must be supplied in strictly
c          ascending order. unchanged on exit.
c    y   : real array of dimension at least (m). before entry, y(i)
c          must be set to the i-th value of the dependent variable y,
c          for i=1,2,...,m. unchanged on exit.
c    w   : real array of dimension at least (m). before entry, w(i)
c          must be set to the i-th value in the set of weights. the
c          w(i) must be strictly positive. unchanged on exit.
c    v   : real array of dimension at least (m). before entry, v(i)
c          must be set to 1 if s(x) must be locally concave at x(i),
c          to (-1) if s(x) must be locally convex at x(i) and to 0
c          if no convexity constraint is imposed at x(i).
c    s   : real. on entry s must specify an over-estimate for the
c          the weighted sum of squared residuals sq of the requested
c          spline. s >=0. unchanged on exit.
c   nest : integer. on entry nest must contain an over-estimate of the
c          total number of knots of the spline returned, to indicate
c          the storage space available to the routine. nest >=8.
c          in most practical situation nest=m/2 will be sufficient.
c          always large enough is  nest=m+4. unchanged on exit.
c  maxtr : integer. on entry maxtr must contain an over-estimate of the
c          total number of records in the used tree structure, to indic-
c          ate the storage space available to the routine. maxtr >=1
c          in most practical situation maxtr=100 will be sufficient.
c          always large enough is
c                         nest-5      nest-6
c              maxtr =  (       ) + (        )  with l the greatest
c                           l          l+1
c          integer <= (nest-6)/2 . unchanged on exit.
c  maxbin: integer. on entry maxbin must contain an over-estimate of the
c          number of knots where s(x) will have a zero second derivative
c          maxbin >=1. in most practical situation maxbin = 10 will be
c          sufficient. always large enough is maxbin=nest-6.
c          unchanged on exit.
c    n   : integer.
c          on exit with ier <=0, n will contain the total number of
c          knots of the spline approximation returned. if the comput-
c          ation mode iopt=1 is used this value of n should be left
c          unchanged between subsequent calls.
c    t   : real array of dimension at least (nest).
c          on exit with ier<=0, this array will contain the knots of the
c          spline,i.e. the position of the interior knots t(5),t(6),...,
c          t(n-4) as well as the position of the additional knots
c          t(1)=t(2)=t(3)=t(4)=x(1) and t(n-3)=t(n-2)=t(n-1)=t(n)=x(m)
c          needed for the b-spline representation.
c          if the computation mode iopt=1 is used, the values of t(1),
c          t(2),...,t(n) should be left unchanged between subsequent
c          calls.
c    c   : real array of dimension at least (nest).
c          on successful exit, this array will contain the coefficients
c          c(1),c(2),..,c(n-4) in the b-spline representation of s(x)
c    sq  : real. unless ier>0 , sq contains the weighted sum of
c          squared residuals of the spline approximation returned.
c    sx  : real array of dimension at least m. on exit with ier<=0
c          this array will contain the spline values s(x(i)),i=1,...,m
c          if the computation mode iopt=1 is used, the values of sx(1),
c          sx(2),...,sx(m) should be left unchanged between subsequent
c          calls.
c    bind: logical array of dimension at least nest. on exit with ier<=0
c          this array will indicate the knots where s''(x)=0, i.e.
c                s''(t(j+3)) .eq. 0 if  bind(j) = .true.
c                s''(t(j+3)) .ne. 0 if  bind(j) = .false., j=1,2,...,n-6
c          if the computation mode iopt=1 is used, the values of bind(1)
c          ,...,bind(n-6) should be left unchanged between subsequent
c          calls.
c   wrk  : real array of dimension at least (m*4+nest*8+maxbin*(maxbin+
c          nest+1)). used as working space.
c   lwrk : integer. on entry,lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.lwrk
c          must not be too small (see wrk). unchanged on exit.
c   iwrk : integer array of dimension at least (maxtr*4+2*(maxbin+1))
c          used as working space.
c   kwrk : integer. on entry,kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program. kwrk
c          must not be too small (see iwrk). unchanged on exit.
c   ier   : integer. error flag
c      ier=0 : normal return, s(x) satisfies the concavity/convexity
c              constraints and sq <= s.
c      ier<0 : abnormal termination: s(x) satisfies the concavity/
c              convexity constraints but sq > s.
c        ier=-3 : the requested storage space exceeds the available
c                 storage space as specified by the parameter nest.
c                 probably causes: nest too small. if nest is already
c                 large (say nest > m/2), it may also indicate that s
c                 is too small.
c                 the approximation returned is the least-squares cubic
c                 spline according to the knots t(1),...,t(n) (n=nest)
c                 which satisfies the convexity constraints.
c        ier=-2 : the maximal number of knots n=m+4 has been reached.
c                 probably causes: s too small.
c        ier=-1 : the number of knots n is less than the maximal number
c                 m+4 but concon finds that adding one or more knots
c                 will not further reduce the value of sq.
c                 probably causes : s too small.
c      ier>0 : abnormal termination: no approximation is returned
c        ier=1  : the number of knots where s''(x)=0 exceeds maxbin.
c                 probably causes : maxbin too small.
c        ier=2  : the number of records in the tree structure exceeds
c                 maxtr.
c                 probably causes : maxtr too small.
c        ier=3  : the algorithm finds no solution to the posed quadratic
c                 programming problem.
c                 probably causes : rounding errors.
c        ier=4  : the minimum number of knots (given by n) to guarantee
c                 that the concavity/convexity conditions will be
c                 satisfied is greater than nest.
c                 probably causes: nest too small.
c        ier=5  : the minimum number of knots (given by n) to guarantee
c                 that the concavity/convexity conditions will be
c                 satisfied is greater than m+4.
c                 probably causes: strongly alternating convexity and
c                 concavity conditions. normally the situation can be
c                 coped with by adding n-m-4 extra data points (found
c                 by linear interpolation e.g.) with a small weight w(i)
c                 and a v(i) number equal to zero.
c        ier=10 : on entry, the input data are controlled on validity.
c                 the following restrictions must be satisfied
c                   0<=iopt<=1, m>3, nest>=8, s>=0, maxtr>=1, maxbin>=1,
c                   kwrk>=maxtr*4+2*(maxbin+1), w(i)>0, x(i) < x(i+1),
c                   lwrk>=m*4+nest*8+maxbin*(maxbin+nest+1)
c                 if one of these restrictions is found to be violated
c                 control is immediately repassed to the calling program
c
c  further comments:
c    as an example of the use of the computation mode iopt=1, the
c    following program segment will cause concon to return control
c    each time a spline with a new set of knots has been computed.
c     .............
c     iopt = 0
c     s = 0.1e+60  (s very large)
c     do 10 i=1,m
c       call concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx,
c    *  bind,wrk,lwrk,iwrk,kwrk,ier)
c       ......
c       s = sq
c       iopt=1
c 10  continue
c     .............
c
c  other subroutines required:
c    fpcoco,fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno
c
c  references:
c   dierckx p. : an algorithm for cubic spline fitting with convexity
c                constraints, computing 24 (1980) 349-371.
c   dierckx p. : an algorithm for least-squares cubic spline fitting
c                with convexity and concavity constraints, report tw39,
c                dept. computer science, k.u.leuven, 1978.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c   p. dierckx
c   dept. computer science, k.u.leuven
c   celestijnenlaan 200a, b-3001 heverlee, belgium.
c   e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : march 1978
c  latest update : march 1987.
c
c  ..
c  ..scalar arguments..
      real*8 s,sq
      integer iopt,m,nest,maxtr,maxbin,n,lwrk,kwrk,ier
c  ..array arguments..
      real*8 x(m),y(m),w(m),v(m),t(nest),c(nest),sx(m),wrk(lwrk)
      integer iwrk(kwrk)
      logical bind(nest)
c  ..local scalars..
      integer i,lwest,kwest,ie,iw,lww
      real*8 one
c  ..
c  set constant
      one = 0.1e+01
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt.lt.0 .or. iopt.gt.1) go to 30
      if(m.lt.4 .or. nest.lt.8) go to 30
      if(s.lt.0.) go to 30
      if(maxtr.lt.1 .or. maxbin.lt.1) go to 30
      lwest = 8*nest+m*4+maxbin*(1+nest+maxbin)
      kwest = 4*maxtr+2*(maxbin+1)
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 30
      if(iopt.gt.0) go to 20
      if(w(1).le.0.) go to 30
      if(v(1).gt.0.) v(1) = one
      if(v(1).lt.0.) v(1) = -one
      do 10 i=2,m
         if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 30
         if(v(i).gt.0.) v(i) = one
         if(v(i).lt.0.) v(i) = -one
  10  continue
  20  ier = 0
c  we partition the working space and determine the spline approximation
      ie = 1
      iw = ie+nest
      lww = lwrk-nest
      call fpcoco(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx,
     * bind,wrk(ie),wrk(iw),lww,iwrk,kwrk,ier)
  30  return
      end
