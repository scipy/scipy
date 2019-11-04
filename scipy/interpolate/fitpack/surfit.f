      subroutine surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c given the set of data points (x(i),y(i),z(i)) and the set of positive
c numbers w(i),i=1,...,m, subroutine surfit determines a smooth bivar-
c iate spline approximation s(x,y) of degrees kx and ky on the rect-
c angle xb <= x <= xe, yb <= y <= ye.
c if iopt = -1 surfit calculates the weighted least-squares spline
c according to a given set of knots.
c if iopt >= 0 the total numbers nx and ny of these knots and their
c position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
c ally by the routine. the smoothness of s(x,y) is then achieved by
c minimalizing the discontinuity jumps in the derivatives of s(x,y)
c across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
c the amounth of smoothness is determined by the condition that f(p) =
c sum ((w(i)*(z(i)-s(x(i),y(i))))**2) be <= s, with s a given non-neg-
c ative constant, called the smoothing factor.
c the fit is given in the b-spline representation (b-spline coefficients
c c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
c uated by means of subroutine bispev.
c
c calling sequence:
c     call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
c    *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer flag. on entry iopt must specify whether a weighted
c          least-squares spline (iopt=-1) or a smoothing spline (iopt=0
c          or 1) must be determined.
c          if iopt=0 the routine will start with an initial set of knots
c          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
c          1,...,ky+1. if iopt=1 the routine will continue with the set
c          of knots found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately pre-
c                     ceded by another call with iopt=1 or iopt=0.
c          unchanged on exit.
c  m     : integer. on entry m must specify the number of data points.
c          m >= (kx+1)*(ky+1). unchanged on exit.
c  x     : real array of dimension at least (m).
c  y     : real array of dimension at least (m).
c  z     : real array of dimension at least (m).
c          before entry, x(i),y(i),z(i) must be set to the co-ordinates
c          of the i-th data point, for i=1,...,m. the order of the data
c          points is immaterial. unchanged on exit.
c  w     : real array of dimension at least (m). before entry, w(i) must
c          be set to the i-th value in the set of weights. the w(i) must
c          be strictly positive. unchanged on exit.
c  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
c  yb,ye   aries of the rectangular approximation domain.
c          xb<=x(i)<=xe,yb<=y(i)<=ye,i=1,...,m. unchanged on exit.
c  kx,ky : integer values. on entry kx and ky must specify the degrees
c          of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
c          (kx=ky=3) splines. unchanged on exit.
c  s     : real. on entry (in case iopt>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nxest : integer. unchanged on exit.
c  nyest : integer. unchanged on exit.
c          on entry, nxest and nyest must specify an upper bound for the
c          number of knots required in the x- and y-directions respect.
c          these numbers will also determine the storage space needed by
c          the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
c          in most practical situation nxest = kx+1+sqrt(m/2), nyest =
c          ky+1+sqrt(m/2) will be sufficient. see also further comments.
c  nmax  : integer. on entry nmax must specify the actual dimension of
c          the arrays tx and ty. nmax >= nxest, nmax >=nyest.
c          unchanged on exit.
c  eps   : real.
c          on entry, eps must specify a threshold for determining the
c          effective rank of an over-determined linear system of equat-
c          ions. 0 < eps < 1.  if the number of decimal digits in the
c          computer representation of a real number is q, then 10**(-q)
c          is a suitable value for eps in most practical applications.
c          unchanged on exit.
c  nx    : integer.
c          unless ier=10 (in case iopt >=0), nx will contain the total
c          number of knots with respect to the x-variable, of the spline
c          approximation returned. if the computation mode iopt=1 is
c          used, the value of nx should be left unchanged between sub-
c          sequent calls.
c          in case iopt=-1, the value of nx should be specified on entry
c  tx    : real array of dimension nmax.
c          on successful exit, this array will contain the knots of the
c          spline with respect to the x-variable, i.e. the position of
c          the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
c          position of the additional knots tx(1)=...=tx(kx+1)=xb and
c          tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
c          if the computation mode iopt=1 is used, the values of tx(1),
c          ...,tx(nx) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tx(kx+2),
c          ...tx(nx-kx-1) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  ny    : integer.
c          unless ier=10 (in case iopt >=0), ny will contain the total
c          number of knots with respect to the y-variable, of the spline
c          approximation returned. if the computation mode iopt=1 is
c          used, the value of ny should be left unchanged between sub-
c          sequent calls.
c          in case iopt=-1, the value of ny should be specified on entry
c  ty    : real array of dimension nmax.
c          on successful exit, this array will contain the knots of the
c          spline with respect to the y-variable, i.e. the position of
c          the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
c          position of the additional knots ty(1)=...=ty(ky+1)=yb and
c          ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
c          if the computation mode iopt=1 is used, the values of ty(1),
c          ...,ty(ny) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values ty(ky+2),
c          ...ty(ny-ky-1) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
c          on successful exit, c contains the coefficients of the spline
c          approximation s(x,y)
c  fp    : real. unless ier=10, fp contains the weighted sum of
c          squared residuals of the spline approximation returned.
c  wrk1  : real array of dimension (lwrk1). used as workspace.
c          if the computation mode iopt=1 is used the value of wrk1(1)
c          should be left unchanged between subsequent calls.
c          on exit wrk1(2),wrk1(3),...,wrk1(1+(nx-kx-1)*(ny-ky-1)) will
c          contain the values d(i)/max(d(i)),i=1,...,(nx-kx-1)*(ny-ky-1)
c          with d(i) the i-th diagonal element of the reduced triangular
c          matrix for calculating the b-spline coefficients. it includes
c          those elements whose square is less than eps,which are treat-
c          ed as 0 in the case of presumed rank deficiency (ier<-2).
c  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
c          the array wrk1 as declared in the calling (sub)program.
c          lwrk1 must not be too small. let
c            u = nxest-kx-1, v = nyest-ky-1, km = max(kx,ky)+1,
c            ne = max(nxest,nyest), bx = kx*v+ky+1, by = ky*u+kx+1,
c            if(bx.le.by) b1 = bx, b2 = b1+v-ky
c            if(bx.gt.by) b1 = by, b2 = b1+u-kx  then
c          lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
c  wrk2  : real array of dimension (lwrk2). used as workspace, but
c          only in the case a rank deficient system is encountered.
c  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
c          the array wrk2 as declared in the calling (sub)program.
c          lwrk2 > 0 . a save upper boundfor lwrk2 = u*v*(b2+1)+b2
c          where u,v and b2 are as above. if there are enough data
c          points, scattered uniformly over the approximation domain
c          and if the smoothing factor s is not too small, there is a
c          good chance that this extra workspace is not needed. a lot
c          of memory might therefore be saved by setting lwrk2=1.
c          (see also ier > 10)
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1).
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is an interpolating
c            spline (fp=0).
c   ier=-2 : normal return. the spline returned is the weighted least-
c            squares polynomial of degrees kx and ky. in this extreme
c            case fp gives the upper bound for the smoothing factor s.
c   ier<-2 : warning. the coefficients of the spline returned have been
c            computed as the minimal norm least-squares solution of a
c            (numerically) rank deficient system. (-ier) gives the rank.
c            especially if the rank deficiency which can be computed as
c            (nx-kx-1)*(ny-ky-1)+ier, is large the results may be inac-
c            curate. they could also seriously depend on the value of
c            eps.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nxest and
c            nyest.
c            probably causes : nxest or nyest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the weighted least-squares
c            spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration process for finding a smoothing spline with
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
c   ier=4  : error. no more knots can be added because the number of
c            b-spline coefficients (nx-kx-1)*(ny-ky-1) already exceeds
c            the number of data points m.
c            probably causes : either s or m too small.
c            the approximation returned is the weighted least-squares
c            spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=5  : error. no more knots can be added because the additional
c            knot would (quasi) coincide with an old one.
c            probably causes : s too small or too large a weight to an
c            inaccurate data point.
c            the approximation returned is the weighted least-squares
c            spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt<=1, 1<=kx,ky<=5, m>=(kx+1)*(ky+1), nxest>=2*kx+2,
c            nyest>=2*ky+2, 0<eps<1, nmax>=nxest, nmax>=nyest,
c            xb<=x(i)<=xe, yb<=y(i)<=ye, w(i)>0, i=1,...,m
c            lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
c            kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1)
c            if iopt=-1: 2*kx+2<=nx<=nxest
c                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
c                        2*ky+2<=ny<=nyest
c                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
c            if iopt>=0: s>=0
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
c            space for computing the minimal least-squares solution of
c            a rank deficient system of linear equations. ier gives the
c            requested value for lwrk2. there is no approximation re-
c            turned but, having saved the information contained in nx,
c            ny,tx,ty,wrk1, and having adjusted the value of lwrk2 and
c            the dimension of the array wrk2 accordingly, the user can
c            continue at the point the program was left, by calling
c            surfit with iopt=1.
c
c further comments:
c  by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the weighted least-squares polynomial (degrees kx,ky)if s is
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
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if surfit is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   surfit once more with the selected value for s but now with iopt=0.
c   indeed, surfit may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nxest and
c   nyest. indeed, if at a certain stage in surfit the number of knots
c   in one direction (say nx) has reached the value of its upper bound
c   (nxest), then from that moment on all subsequent knots are added
c   in the other (y) direction. this may indicate that the value of
c   nxest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting nxest=2*kx+2 (the lowest allowable value for
c   nxest), the user can indicate that he wants an approximation which
c   is a simple polynomial of degree kx in the variable x.
c
c  other subroutines required:
c    fpback,fpbspl,fpsurf,fpdisc,fpgivs,fprank,fprati,fprota,fporde
c
c  references:
c   dierckx p. : an algorithm for surface fitting with spline functions
c                ima j. numer. anal. 1 (1981) 267-283.
c   dierckx p. : an algorithm for surface fitting with spline functions
c                report tw50, dept. computer science,k.u.leuven, 1980.
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
      real*8 xb,xe,yb,ye,s,eps,fp
      integer iopt,m,kx,ky,nxest,nyest,nmax,nx,ny,lwrk1,lwrk2,kwrk,ier
c  ..array arguments..
      real*8 x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),
     * c((nxest-kx-1)*(nyest-ky-1)),wrk1(lwrk1),wrk2(lwrk2)
      integer iwrk(kwrk)
c  ..local scalars..
      real*8 tol
      integer i,ib1,ib3,jb1,ki,kmax,km1,km2,kn,kwest,kx1,ky1,la,lbx,
     * lby,lco,lf,lff,lfp,lh,lq,lsx,lsy,lwest,maxit,ncest,nest,nek,
     * nminx,nminy,nmx,nmy,nreg,nrint,nxk,nyk
c  ..function references..
      integer max0
c  ..subroutine references..
c    fpsurf
c  ..
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid,control is immediately repassed to the calling program.
      ier = 10
      if(eps.le.0. .or. eps.ge.1.) go to 71
      if(kx.le.0 .or. kx.gt.5) go to 71
      kx1 = kx+1
      if(ky.le.0 .or. ky.gt.5) go to 71
      ky1 = ky+1
      kmax = max0(kx,ky)
      km1 = kmax+1
      km2 = km1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 71
      if(m.lt.(kx1*ky1)) go to 71
      nminx = 2*kx1
      if(nxest.lt.nminx .or. nxest.gt.nmax) go to 71
      nminy = 2*ky1
      if(nyest.lt.nminy .or. nyest.gt.nmax) go to 71
      nest = max0(nxest,nyest)
      nxk = nxest-kx1
      nyk = nyest-ky1
      ncest = nxk*nyk
      nmx = nxest-nminx+1
      nmy = nyest-nminy+1
      nrint = nmx+nmy
      nreg = nmx*nmy
      ib1 = kx*nyk+ky1
      jb1 = ky*nxk+kx1
      ib3 = kx1*nyk+1
      if(ib1.le.jb1) go to 10
      ib1 = jb1
      ib3 = ky1*nxk+1
  10  lwest = ncest*(2+ib1+ib3)+2*(nrint+nest*km2+m*km1)+ib3
      kwest = m+nreg
      if(lwrk1.lt.lwest .or. kwrk.lt.kwest) go to 71
      if(xb.ge.xe .or. yb.ge.ye) go to 71
      do 20 i=1,m
        if(w(i).le.0.) go to 70
        if(x(i).lt.xb .or. x(i).gt.xe) go to 71
        if(y(i).lt.yb .or. y(i).gt.ye) go to 71
  20  continue
      if(iopt.ge.0) go to 50
      if(nx.lt.nminx .or. nx.gt.nxest) go to 71
      nxk = nx-kx1
      tx(kx1) = xb
      tx(nxk+1) = xe
      do 30 i=kx1,nxk
        if(tx(i+1).le.tx(i)) go to 72
  30  continue
      if(ny.lt.nminy .or. ny.gt.nyest) go to 71
      nyk = ny-ky1
      ty(ky1) = yb
      ty(nyk+1) = ye
      do 40 i=ky1,nyk
        if(ty(i+1).le.ty(i)) go to 73
  40  continue
      go to 60
  50  if(s.lt.0.) go to 71
  60  ier = 0
c  we partition the working space and determine the spline approximation
      kn = 1
      ki = kn+m
      lq = 2
      la = lq+ncest*ib3
      lf = la+ncest*ib1
      lff = lf+ncest
      lfp = lff+ncest
      lco = lfp+nrint
      lh = lco+nrint
      lbx = lh+ib3
      nek = nest*km2
      lby = lbx+nek
      lsx = lby+nek
      lsy = lsx+m*km1
      call fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     * eps,tol,maxit,nest,km1,km2,ib1,ib3,ncest,nrint,nreg,nx,tx,
     * ny,ty,c,fp,wrk1(1),wrk1(lfp),wrk1(lco),wrk1(lf),wrk1(lff),
     * wrk1(la),wrk1(lq),wrk1(lbx),wrk1(lby),wrk1(lsx),wrk1(lsy),
     * wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
 70   return
 71   print*,"iopt,kx,ky,m=",iopt,kx,ky,m
      print*,"nxest,nyest,nmax=",nxest,nyest,nmax
      print*,"lwrk1,lwrk2,kwrk=",lwrk1,lwrk2,kwrk
      print*,"xb,xe,yb,ye=",xb,xe,yb,ye
      print*,"eps,s",eps,s
      return
 72   print*,"tx=",tx
      return
 73   print*,"ty=",ty
      return
      end
