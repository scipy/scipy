      subroutine regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,
     * nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c given the set of values z(i,j) on the rectangular grid (x(i),y(j)),
c i=1,...,mx;j=1,...,my, subroutine regrid determines a smooth bivar-
c iate spline approximation s(x,y) of degrees kx and ky on the rect-
c angle xb <= x <= xe, yb <= y <= ye.
c if iopt = -1 regrid calculates the least-squares spline according
c to a given set of knots.
c if iopt >= 0 the total numbers nx and ny of these knots and their
c position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
c ally by the routine. the smoothness of s(x,y) is then achieved by
c minimalizing the discontinuity jumps in the derivatives of s(x,y)
c across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
c the amounth of smoothness is determined by the condition that f(p) =
c sum ((z(i,j)-s(x(i),y(j))))**2) be <= s, with s a given non-negative
c constant, called the smoothing factor.
c the fit is given in the b-spline representation (b-spline coefficients
c c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
c uated by means of subroutine bispev.
c
c calling sequence:
c     call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
c    *  nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer flag. on entry iopt must specify whether a least-
c          squares spline (iopt=-1) or a smoothing spline (iopt=0 or 1)
c          must be determined.
c          if iopt=0 the routine will start with an initial set of knots
c          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
c          1,...,ky+1. if iopt=1 the routine will continue with the set
c          of knots found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately pre-
c                     ceded by another call with iopt=1 or iopt=0 and
c                     s.ne.0.
c          unchanged on exit.
c  mx    : integer. on entry mx must specify the number of grid points
c          along the x-axis. mx > kx . unchanged on exit.
c  x     : real array of dimension at least (mx). before entry, x(i)
c          must be set to the x-co-ordinate of the i-th grid point
c          along the x-axis, for i=1,2,...,mx. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c  my    : integer. on entry my must specify the number of grid points
c          along the y-axis. my > ky . unchanged on exit.
c  y     : real array of dimension at least (my). before entry, y(j)
c          must be set to the y-co-ordinate of the j-th grid point
c          along the y-axis, for j=1,2,...,my. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c  z     : real array of dimension at least (mx*my).
c          before entry, z(my*(i-1)+j) must be set to the data value at
c          the grid point (x(i),y(j)) for i=1,...,mx and j=1,...,my.
c          unchanged on exit.
c  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
c  yb,ye   aries of the rectangular approximation domain.
c          xb<=x(i)<=xe,i=1,...,mx; yb<=y(j)<=ye,j=1,...,my.
c          unchanged on exit.
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
c          in most practical situation nxest = mx/2, nyest=my/2, will
c          be sufficient. always large enough are nxest=mx+kx+1, nyest=
c          my+ky+1, the number of knots needed for interpolation (s=0).
c          see also further comments.
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
c  fp    : real. unless ier=10, fp contains the sum of squared
c          residuals of the spline approximation returned.
c  wrk   : real array of dimension (lwrk). used as workspace.
c          if the computation mode iopt=1 is used the values of wrk(1),
c          ...,wrk(4) should be left unchanged between subsequent calls.
c  lwrk  : integer. on entry lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.
c          lwrk must not be too small.
c           lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
c            my*(ky+1) +u
c           where u is the larger of my and nxest.
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c          if the computation mode iopt=1 is used the values of iwrk(1),
c          ...,iwrk(3) should be left unchanged between subsequent calls
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= 3+mx+my+nxest+nyest.
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is an interpolating
c            spline (fp=0).
c   ier=-2 : normal return. the spline returned is the least-squares
c            polynomial of degrees kx and ky. in this extreme case fp
c            gives the upper bound for the smoothing factor s.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nxest and
c            nyest.
c            probably causes : nxest or nyest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the least-squares spline
c            according to the current set of knots. the parameter fp
c            gives the corresponding sum of squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration process for finding a smoothing spline with
c            fp = s. probably causes : s too small.
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing spline
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt<=1, 1<=kx,ky<=5, mx>kx, my>ky, nxest>=2*kx+2,
c            nyest>=2*ky+2, kwrk>=3+mx+my+nxest+nyest,
c            lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
c             my*(ky+1) +max(my,nxest),
c            xb<=x(i-1)<x(i)<=xe,i=2,..,mx,yb<=y(j-1)<y(j)<=ye,j=2,..,my
c            if iopt=-1: 2*kx+2<=nx<=min(nxest,mx+kx+1)
c                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
c                        2*ky+2<=ny<=min(nyest,my+ky+1)
c                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
c                    the schoenberg-whitney conditions, i.e. there must
c                    be subset of grid co-ordinates xx(p) and yy(q) such
c                    that   tx(p) < xx(p) < tx(p+kx+1) ,p=1,...,nx-kx-1
c                           ty(q) < yy(q) < ty(q+ky+1) ,q=1,...,ny-ky-1
c            if iopt>=0: s>=0
c                        if s=0 : nxest>=mx+kx+1, nyest>=my+ky+1
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c
c further comments:
c   regrid does not allow individual weighting of the data-values.
c   so, if these were determined to widely different accuracies, then
c   perhaps the general data set routine surfit should rather be used
c   in spite of efficiency.
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the least-squares polynomial (degrees kx,ky) if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the accuracy of the data values.
c   if the user has an idea of the statistical errors on the data, he
c   can also find a proper estimate for s. for, by assuming that, if he
c   specifies the right s, regrid will return a spline s(x,y) which
c   exactly reproduces the function underlying the data he can evaluate
c   the sum((z(i,j)-s(x(i),y(j)))**2) to find a good estimate for this s
c   for example, if he knows that the statistical errors on his z(i,j)-
c   values is not greater than 0.1, he may expect that a good s should
c   have a value not larger than mx*my*(0.1)**2.
c   if nothing is known about the statistical error in z(i,j), s must
c   be determined by trial and error, taking account of the comments
c   above. the best is then to start with a very large value of s (to
c   determine the least-squares polynomial and the corresponding upper
c   bound fp0 for s) and then to progressively decrease the value of s
c   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
c   and more carefully as the approximation shows more detail) to
c   obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if regrid is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   regrid once more with the selected value for s but now with iopt=0.
c   indeed, regrid may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nxest and
c   nyest. indeed, if at a certain stage in regrid the number of knots
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
c    fpback,fpbspl,fpregr,fpdisc,fpgivs,fpgrre,fprati,fprota,fpchec,
c    fpknot
c
c  references:
c   dierckx p. : a fast algorithm for smoothing data on a rectangular
c                grid while using spline functions, siam j.numer.anal.
c                19 (1982) 1286-1304.
c   dierckx p. : a fast algorithm for smoothing data on a rectangular
c                grid while using spline functions, report tw53, dept.
c                computer science,k.u.leuven, 1980.
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
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      real*8 xb,xe,yb,ye,s,fp
      integer iopt,mx,my,kx,ky,nxest,nyest,nx,ny,lwrk,kwrk,ier
c  ..array arguments..
      real*8 x(mx),y(my),z(mx*my),tx(nxest),ty(nyest),
     * c((nxest-kx-1)*(nyest-ky-1)),wrk(lwrk)
      integer iwrk(kwrk)
c  ..local scalars..
      real*8 tol
      integer i,j,jwrk,kndx,kndy,knrx,knry,kwest,kx1,kx2,ky1,ky2,
     * lfpx,lfpy,lwest,lww,maxit,nc,nminx,nminy,mz
c  ..function references..
      integer max0
c  ..subroutine references..
c    fpregr,fpchec
c  ..
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(kx.le.0 .or. kx.gt.5) go to 70
      kx1 = kx+1
      kx2 = kx1+1
      if(ky.le.0 .or. ky.gt.5) go to 70
      ky1 = ky+1
      ky2 = ky1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 70
      nminx = 2*kx1
      if(mx.lt.kx1 .or. nxest.lt.nminx) go to 70
      nminy = 2*ky1
      if(my.lt.ky1 .or. nyest.lt.nminy) go to 70
      mz = mx*my
      nc = (nxest-kx1)*(nyest-ky1)
      lwest = 4+nxest*(my+2*kx2+1)+nyest*(2*ky2+1)+mx*kx1+
     * my*ky1+max0(nxest,my)
      kwest = 3+mx+my+nxest+nyest
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 70
      if(xb.gt.x(1) .or. xe.lt.x(mx)) go to 70
      do 10 i=2,mx
        if(x(i-1).ge.x(i)) go to 70
  10  continue
      if(yb.gt.y(1) .or. ye.lt.y(my)) go to 70
      do 20 i=2,my
        if(y(i-1).ge.y(i)) go to 70
  20  continue
      if(iopt.ge.0) go to 50
      if(nx.lt.nminx .or. nx.gt.nxest) go to 70
      j = nx
      do 30 i=1,kx1
        tx(i) = xb
        tx(j) = xe
        j = j-1
  30  continue
      call fpchec(x,mx,tx,nx,kx,ier)
      if(ier.ne.0) go to 70
      if(ny.lt.nminy .or. ny.gt.nyest) go to 70
      j = ny
      do 40 i=1,ky1
        ty(i) = yb
        ty(j) = ye
        j = j-1
  40  continue
      call fpchec(y,my,ty,ny,ky,ier)
      if (ier.eq.0) go to 60
      go to 70
  50  if(s.lt.0.) go to 70
      if(s.eq.0. .and. (nxest.lt.(mx+kx1) .or. nyest.lt.(my+ky1)) )
     * go to 70
      ier = 0
c  we partition the working space and determine the spline approximation
  60  lfpx = 5
      lfpy = lfpx+nxest
      lww = lfpy+nyest
      jwrk = lwrk-4-nxest-nyest
      knrx = 4
      knry = knrx+mx
      kndx = knry+my
      kndy = kndx+nxest
      call fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     * tol,maxit,nc,nx,tx,ny,ty,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),
     * wrk(lfpx),wrk(lfpy),iwrk(1),iwrk(2),iwrk(3),iwrk(knrx),
     * iwrk(knry),iwrk(kndx),iwrk(kndy),wrk(lww),jwrk,ier)
  70  return
      end

