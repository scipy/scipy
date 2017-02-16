      subroutine dqagpe(f,a,b,npts2,points,epsabs,epsrel,limit,result,
     *   abserr,neval,ier,alist,blist,rlist,elist,pts,iord,level,ndin,
     *   last)
c***begin prologue  dqagpe
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, general-purpose,
c             singularities at user specified points,
c             extrapolation, globally adaptive.
c***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b), hopefully
c            satisfying following claim for accuracy abs(i-result).le.
c            max(epsabs,epsrel*abs(i)). break points of the integration
c            interval, where local difficulties of the integrand may
c            occur(e.g. singularities,discontinuities),provided by user.
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            npts2  - integer
c                     number equal to two more than the number of
c                     user-supplied break points within the integration
c                     range, npts2.ge.2.
c                     if npts2.lt.2, the routine will end with ier = 6.
c
c            points - double precision
c                     vector of dimension npts2, the first (npts2-2)
c                     elements of which are the user provided break
c                     points. if these points do not constitute an
c                     ascending sequence there will be an automatic
c                     sorting.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.npts2
c                     if limit.lt.npts2, the routine will end with
c                     ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (i.e. singularity,
c                             discontinuity within the interval), it
c                             should be supplied to the routine as an
c                             element of the vector points. if necessary
c                             an appropriate special-purpose integrator
c                             must be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table. it is presumed that
c                             the requested tolerance cannot be
c                             achieved, and that the returned result is
c                             the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.gt.0.
c                         = 6 the input is invalid because
c                             npts2.lt.2 or
c                             break points are specified outside
c                             the integration range or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.npts2.
c                             result, abserr, neval, last, rlist(1),
c                             and elist(1) are set to zero. alist(1) and
c                             blist(1) are set to a and b respectively.
c
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            rlist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            pts    - double precision
c                     vector of dimension at least npts2, containing the
c                     integration limits and the break points of the
c                     interval in ascending sequence.
c
c            level  - integer
c                     vector of dimension at least limit, containing the
c                     subdivision levels of the subinterval, i.e. if
c                     (aa,bb) is a subinterval of (p1,p2) where p1 as
c                     well as p2 is a user-provided break point or
c                     integration limit, then (aa,bb) has level l if
c                     abs(bb-aa) = abs(p2-p1)*2**(-l).
c
c            ndin   - integer
c                     vector of dimension at least npts2, after first
c                     integration over the intervals (pts(i)),pts(i+1),
c                     i = 0,1, ..., npts2-2, the error estimates over
c                     some of the intervals may have been increased
c                     artificially, in order to put their subdivision
c                     forward. if this happens for the subinterval
c                     numbered k, ndin(k) is put to 1, otherwise
c                     ndin(k) = 0.
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced in the
c                     subdivisions process
c
c***references  (none)
c***routines called  d1mach,dqelg,dqk21,dqpsrt
c***end prologue  dqagpe
      double precision a,abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,b,blist,b1,b2,correc,dabs,defabs,defab1,defab2,dmax1,dmin1,
     *  dres,d1mach,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,
     *  errmax,error1,erro12,error2,errsum,ertest,f,oflow,points,pts,
     *  resa,resabs,reseps,result,res3la,rlist,rlist2,sign,temp,uflow
      integer i,id,ier,ierro,ind1,ind2,iord,ip1,iroff1,iroff2,iroff3,j,
     *  jlow,jupbnd,k,ksgn,ktmin,last,levcur,level,levmax,limit,maxerr,
     *  ndin,neval,nint,nintp1,npts,npts2,nres,nrmax,numrl2
      logical extrap,noext
c
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  level(limit),ndin(npts2),points(npts2),pts(npts2),res3la(3),
     *  rlist(limit),rlist2(52)
c
      external f
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine epsalg (rlist2 should be of dimension
c            (limexp+2) at least).
c
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2
c                       containing the part of the epsilon table which
c                       is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements in rlist2. if an appropriate
c                       approximation to the compounded integral has
c                       been obtained, it is put in rlist2(numrl2) after
c                       numrl2 has been increased by one.
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine
c                       is attempting to perform extrapolation. i.e.
c                       before subdividing the smallest interval we
c                       try to decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation is
c                       no longer allowed (true-value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqagpe
      epmach = d1mach(4)
c
c            test on validity of parameters
c            -----------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      level(1) = 0
      npts = npts2-2
      if(npts2.lt.2.or.limit.le.npts.or.(epsabs.le.0.0d+00.and.
     *  epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))) ier = 6
      if(ier.eq.6) go to 999
c
c            if any break points are provided, sort them into an
c            ascending sequence.
c
      sign = 1.0d+00
      if(a.gt.b) sign = -1.0d+00
      pts(1) = dmin1(a,b)
      if(npts.eq.0) go to 15
      do 10 i = 1,npts
        pts(i+1) = points(i)
   10 continue
   15 pts(npts+2) = dmax1(a,b)
      nint = npts+1
      a1 = pts(1)
      if(npts.eq.0) go to 40
      nintp1 = nint+1
      do 20 i = 1,nint
        ip1 = i+1
        do 20 j = ip1,nintp1
          if(pts(i).le.pts(j)) go to 20
          temp = pts(i)
          pts(i) = pts(j)
          pts(j) = temp
   20 continue
      if(pts(1).ne.dmin1(a,b).or.pts(nintp1).ne.dmax1(a,b)) ier = 6
      if(ier.eq.6) go to 999
c
c            compute first integral and error approximations.
c            ------------------------------------------------
c
   40 resabs = 0.0d+00
      do 50 i = 1,nint
        b1 = pts(i+1)
        call dqk21(f,a1,b1,area1,error1,defabs,resa)
        abserr = abserr+error1
        result = result+area1
        ndin(i) = 0
        if(error1.eq.resa.and.error1.ne.0.0d+00) ndin(i) = 1
        resabs = resabs+defabs
        level(i) = 0
        elist(i) = error1
        alist(i) = a1
        blist(i) = b1
        rlist(i) = area1
        iord(i) = i
        a1 = b1
   50 continue
      errsum = 0.0d+00
      do 55 i = 1,nint
        if(ndin(i).eq.1) elist(i) = abserr
        errsum = errsum+elist(i)
   55 continue
c
c           test on accuracy.
c
      last = nint
      neval = 21*nint
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      if(abserr.le.0.1d+03*epmach*resabs.and.abserr.gt.errbnd) ier = 2
      if(nint.eq.1) go to 80
      do 70 i = 1,npts
        jlow = i+1
        ind1 = iord(i)
        do 60 j = jlow,nint
          ind2 = iord(j)
          if(elist(ind1).gt.elist(ind2)) go to 60
          ind1 = ind2
          k = j
   60   continue
        if(ind1.eq.iord(i)) go to 70
        iord(k) = iord(i)
        iord(i) = ind1
   70 continue
      if(limit.lt.npts2) ier = 1
   80 if(ier.ne.0.or.abserr.le.errbnd) go to 210
c
c           initialization
c           --------------
c
      rlist2(1) = result
      maxerr = iord(1)
      errmax = elist(maxerr)
      area = result
      nrmax = 1
      nres = 0
      numrl2 = 1
      ktmin = 0
      extrap = .false.
      noext = .false.
      erlarg = errsum
      ertest = errbnd
      levmax = 1
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ierro = 0
      uflow = d1mach(1)
      oflow = d1mach(2)
      abserr = oflow
      ksgn = -1
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*resabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 160 last = npts2,limit
c
c           bisect the subinterval with the nrmax-th largest error
c           estimate.
c
        levcur = level(maxerr)+1
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dqk21(f,a1,b1,area1,error1,resa,defab1)
        call dqk21(f,a2,b2,area2,error2,resa,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        neval = neval+42
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 95
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     *  .or.erro12.lt.0.99d+00*errmax) go to 90
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   90   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   95   level(maxerr) = levcur
        level(last) = levcur
        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     *  (dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 100
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 110
  100   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
  110   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(errsum.le.errbnd) go to 190
c ***jump out of do-loop
        if(ier.ne.0) go to 170
        if(noext) go to 160
        erlarg = erlarg-erlast
        if(levcur+1.le.levmax) erlarg = erlarg+erro12
        if(extrap) go to 120
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(level(maxerr)+1.le.levmax) go to 160
        extrap = .true.
        nrmax = 2
  120   if(ierro.eq.3.or.erlarg.le.ertest) go to 140
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over
c           the larger intervals (erlarg) and perform extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 130 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
c ***jump out of do-loop
          if(level(maxerr)+1.le.levmax) go to 160
          nrmax = nrmax+1
  130   continue
c
c           perform extrapolation.
c
  140   numrl2 = numrl2+1
        rlist2(numrl2) = area
        if(numrl2.le.2) go to 155
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 150
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
c ***jump out of do-loop
        if(abserr.lt.ertest) go to 170
c
c           prepare bisection of the smallest interval.
c
  150   if(numrl2.eq.1) noext = .true.
        if(ier.ge.5) go to 170
  155   maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        levmax = levmax+1
        erlarg = errsum
  160 continue
c
c           set the final result.
c           ---------------------
c
c
  170 if(abserr.eq.oflow) go to 190
      if((ier+ierro).eq.0) go to 180
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00)go to 175
      if(abserr.gt.errsum)go to 190
      if(area.eq.0.0d+00) go to 210
      go to 180
  175 if(abserr/dabs(result).gt.errsum/dabs(area))go to 190
c
c           test on divergence.
c
  180 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     *  resabs*0.1d-01) go to 210
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03.or.
     *  errsum.gt.dabs(area)) ier = 6
      go to 210
c
c           compute global integral sum.
c
  190 result = 0.0d+00
      do 200 k = 1,last
        result = result+rlist(k)
  200 continue
      abserr = errsum
  210 if(ier.gt.2) ier = ier-1
      result = result*sign
  999 return
      end
