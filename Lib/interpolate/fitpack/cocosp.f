      subroutine cocosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,wrk,
     * lwrk,iwrk,kwrk,ier)
c  given the set of data points (x(i),y(i)) and the set of positive
c  numbers w(i),i=1,2,...,m, subroutine cocosp determines the weighted
c  least-squares cubic spline s(x) with given knots t(j),j=1,2,...,n
c  which satisfies the following concavity/convexity conditions
c      s''(t(j+3))*e(j) <= 0, j=1,2,...n-6
c  the fit is given in the b-spline representation( b-spline coef-
c  ficients c(j),j=1,2,...n-4) and can be evaluated by means of
c  subroutine splev.
c
c  calling sequence:
c     call cocosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,wrk,
c    * lwrk,iwrk,kwrk,ier)
c
c  parameters:
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
c    n   : integer. on entry n must contain the total number of knots
c          of the cubic spline. m+4>=n>=8. unchanged on exit.
c    t   : real array of dimension at least (n). before entry, this
c          array must contain the knots of the spline, i.e. the position
c          of the interior knots t(5),t(6),...,t(n-4) as well as the
c          position of the boundary knots t(1),t(2),t(3),t(4) and t(n-3)
c          t(n-2),t(n-1),t(n) needed for the b-spline representation.
c          unchanged on exit. see also the restrictions (ier=10).
c    e   : real array of dimension at least (n). before entry, e(j)
c          must be set to 1 if s(x) must be locally concave at t(j+3),
c          to (-1) if s(x) must be locally convex at t(j+3) and to 0
c          if no convexity constraint is imposed at t(j+3),j=1,2,..,n-6.
c          e(n-5),...,e(n) are not used. unchanged on exit.
c  maxtr : integer. on entry maxtr must contain an over-estimate of the
c          total number of records in the used tree structure, to indic-
c          ate the storage space available to the routine. maxtr >=1
c          in most practical situation maxtr=100 will be sufficient.
c          always large enough is
c                         n-5       n-6
c              maxtr =  (     ) + (     )  with l the greatest
c                          l        l+1
c          integer <= (n-6)/2 . unchanged on exit.
c  maxbin: integer. on entry maxbin must contain an over-estimate of the
c          number of knots where s(x) will have a zero second derivative
c          maxbin >=1. in most practical situation maxbin = 10 will be
c          sufficient. always large enough is maxbin=n-6.
c          unchanged on exit.
c    c   : real array of dimension at least (n).
c          on succesful exit, this array will contain the coefficients
c          c(1),c(2),..,c(n-4) in the b-spline representation of s(x)
c    sq  : real. on succesful exit, sq contains the weighted sum of
c          squared residuals of the spline approximation returned.
c    sx  : real array of dimension at least m. on succesful exit
c          this array will contain the spline values s(x(i)),i=1,...,m
c   bind : logical array of dimension at least (n). on succesful exit
c          this array will indicate the knots where s''(x)=0, i.e.
c                s''(t(j+3)) .eq. 0 if  bind(j) = .true.
c                s''(t(j+3)) .ne. 0 if  bind(j) = .false., j=1,2,...,n-6
c   wrk  : real array of dimension at least  m*4+n*7+maxbin*(maxbin+n+1)
c          used as working space.
c   lwrk : integer. on entry,lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.lwrk
c          must not be too small (see wrk). unchanged on exit.
c   iwrk : integer array of dimension at least (maxtr*4+2*(maxbin+1))
c          used as working space.
c   kwrk : integer. on entry,kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program. kwrk
c          must not be too small (see iwrk). unchanged on exit.
c   ier   : integer. error flag
c      ier=0 : succesful exit.
c      ier>0 : abnormal termination: no approximation is returned
c        ier=1  : the number of knots where s''(x)=0 exceeds maxbin.
c                 probably causes : maxbin too small.
c        ier=2  : the number of records in the tree structure exceeds
c                 maxtr.
c                 probably causes : maxtr too small.
c        ier=3  : the algoritm finds no solution to the posed quadratic
c                 programming problem.
c                 probably causes : rounding errors.
c        ier=10 : on entry, the input data are controlled on validity.
c                 the following restrictions must be satisfied
c                   m>3, maxtr>=1, maxbin>=1, 8<=n<=m+4,w(i) > 0,
c                   x(1)<x(2)<...<x(m), t(1)<=t(2)<=t(3)<=t(4)<=x(1),
c                   x(1)<t(5)<t(6)<...<t(n-4)<x(m)<=t(n-3)<=...<=t(n),
c                   kwrk>=maxtr*4+2*(maxbin+1),
c                   lwrk>=m*4+n*7+maxbin*(maxbin+n+1),
c                   the schoenberg-whitney conditions, i.e. there must
c                   be a subset of data points xx(j) such that
c                     t(j) < xx(j) < t(j+4), j=1,2,...,n-4
c                 if one of these restrictions is found to be violated
c                 control is immediately repassed to the calling program
c
c
c  other subroutines required:
c    fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno,fpchec
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
      real*8 sq
      integer m,n,maxtr,maxbin,lwrk,kwrk,ier
c  ..array arguments..
      real*8 x(m),y(m),w(m),t(n),e(n),c(n),sx(m),wrk(lwrk)
      integer iwrk(kwrk)
      logical bind(n)
c  ..local scalars..
      integer i,ia,ib,ic,iq,iu,iz,izz,ji,jib,jjb,jl,jr,ju,kwest,
     * lwest,mb,nm,n6
      real*8 one
c  ..
c  set constant
      one = 0.1e+01
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(m.lt.4 .or. n.lt.8) go to 40
      if(maxtr.lt.1 .or. maxbin.lt.1) go to 40
      lwest = 7*n+m*4+maxbin*(1+n+maxbin)
      kwest = 4*maxtr+2*(maxbin+1)
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 40
      if(w(1).le.0.) go to 40
      do 10 i=2,m
         if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 40
  10  continue
      call fpchec(x,m,t,n,3,ier)
      if(ier) 40,20,40
c  set numbers e(i)
  20  n6 = n-6
      do 30 i=1,n6
        if(e(i).gt.0.) e(i) = one
        if(e(i).lt.0.) e(i) = -one
  30  continue
c  we partition the working space and determine the spline approximation
      nm = n+maxbin
      mb = maxbin+1
      ia = 1
      ib = ia+4*n
      ic = ib+nm*maxbin
      iz = ic+n
      izz = iz+n
      iu = izz+n
      iq = iu+maxbin
      ji = 1
      ju = ji+maxtr
      jl = ju+maxtr
      jr = jl+maxtr
      jjb = jr+maxtr
      jib = jjb+mb
      call fpcosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,nm,mb,wrk(ia),
     *
     * wrk(ib),wrk(ic),wrk(iz),wrk(izz),wrk(iu),wrk(iq),iwrk(ji),
     * iwrk(ju),iwrk(jl),iwrk(jr),iwrk(jjb),iwrk(jib),ier)
  40  return
      end
