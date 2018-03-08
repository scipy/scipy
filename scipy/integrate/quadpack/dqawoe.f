      subroutine dqawoe (f,a,b,omega,integr,epsabs,epsrel,limit,icall,
     *  maxp1,result,abserr,neval,ier,last,alist,blist,rlist,elist,iord,
     *   nnlog,momcom,chebmo)
c***begin prologue  dqawoe
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             integrand with oscillatory cos or sin factor,
c             clenshaw-curtis method, (end point) singularities,
c             extrapolation, globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral
c            i = integral of f(x)*w(x) over (a,b)
c            where w(x) = cos(omega*x) or w(x)=sin(omega*x),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of oscillatory integrals
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
c            omega  - double precision
c                     parameter in the integrand weight function
c
c            integr - integer
c                     indicates which of the weight functions is to be
c                     used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1 and integr.ne.2, the routine
c                     will end with ier = 6.
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
c                     gives an upper bound on the number of subdivisions
c                     in the partition of (a,b), limit.ge.1.
c
c            icall  - integer
c                     if dqawoe is to be used only once, icall must
c                     be set to 1.  assume that during this call, the
c                     chebyshev moments (for clenshaw-curtis integration
c                     of degree 24) have been computed for intervals of
c                     lengths (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
c                     if icall.gt.1 this means that dqawoe has been
c                     called twice or more on intervals of the same
c                     length abs(b-a). the chebyshev moments already
c                     computed are then re-used in subsequent calls.
c                     if icall.lt.1, the routine will end with ier = 6.
c
c            maxp1  - integer
c                     gives an upper bound on the number of chebyshev
c                     moments which can be stored, i.e. for the
c                     intervals of lengths abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1.
c                     if maxp1.lt.1, the routine will end with ier = 6.
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
c                             routine. it is assumed that the
c                             requested accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand, in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty can
c                             be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used which is designed for
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
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved due to
c                             roundoff in the extrapolation table,
c                             and that the returned result is the
c                             best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.gt.0.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or (integr.ne.1 and integr.ne.2) or
c                             icall.lt.1 or maxp1.lt.1.
c                             result, abserr, neval, last, rlist(1),
c                             elist(1), iord(1) and nnlog(1) are set
c                             to zero. alist(1) and blist(1) are set
c                             to a and b respectively.
c
c            last  -  integer
c                     on return, last equals the number of
c                     subintervals produces in the subdivision
c                     process, which determines the number of
c                     significant elements actually in the
c                     work arrays.
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
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
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the error
c                     estimates over the subintervals,
c                     such that elist(iord(1)), ...,
c                     elist(iord(k)) form a decreasing sequence, with
c                     k = last if last.le.(limit/2+2), and
c                     k = limit+1-last otherwise.
c
c            nnlog  - integer
c                     vector of dimension at least limit, containing the
c                     subdivision levels of the subintervals, i.e.
c                     iwork(i) = l means that the subinterval
c                     numbered i is of length abs(b-a)*2**(1-l)
c
c         on entry and return
c            momcom - integer
c                     indicating that the chebyshev moments
c                     have been computed for intervals of lengths
c                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
c                     momcom.lt.maxp1
c
c            chebmo - double precision
c                     array of dimension (maxp1,25) containing the
c                     chebyshev moments
c
c***references  (none)
c***routines called  d1mach,dqc25f,dqelg,dqpsrt
c***end prologue  dqawoe
c
      double precision a,abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,b,blist,b1,b2,chebmo,correc,dabs,defab1,defab2,defabs,dmax1,
     *  domega,d1mach,dres,elist,epmach,epsabs,epsrel,erlarg,erlast,
     *  errbnd,errmax,error1,erro12,error2,errsum,ertest,f,oflow,
     *  omega,resabs,reseps,result,res3la,rlist,rlist2,small,uflow,width
      integer icall,id,ier,ierro,integr,iord,iroff1,iroff2,iroff3,
     *  jupbnd,k,ksgn,ktmin,last,limit,maxerr,maxp1,momcom,nev,neval,
     *  nnlog,nres,nrmax,nrmom,numrl2
      logical extrap,noext,extall
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit),rlist2(52),res3la(3),chebmo(maxp1,25),nnlog(limit)
c
      external f
c
c            the dimension of rlist2 is determined by  the value of
c            limexp in subroutine dqelg (rlist2 should be of
c            dimension (limexp+2) at least).
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
c                       containing the part of the epsilon table
c                       which is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
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
c                       been obtained it is put in rlist2(numrl2) after
c                       numrl2 has been increased by one
c           small     - length of the smallest interval considered
c                       up to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine is
c                       attempting to perform extrapolation, i.e. before
c                       subdividing the smallest interval we try to
c                       decrease the value of erlarg
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true  value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqawoe
      epmach = d1mach(4)
c
c         test on validity of parameters
c         ------------------------------
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
      nnlog(1) = 0
      if((integr.ne.1.and.integr.ne.2).or.(epsabs.le.0.0d+00.and.
     *  epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28)).or.icall.lt.1.or.
     *  maxp1.lt.1) ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      domega = dabs(omega)
      nrmom = 0
      if (icall.gt.1) go to 5
      momcom = 0
    5 call dqc25f(f,a,b,domega,integr,nrmom,maxp1,0,result,abserr,
     *  neval,defabs,resabs,momcom,chebmo)
c
c           test on accuracy.
c
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if(abserr.le.0.1d+03*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.abserr.le.errbnd) go to 200
c
c           initializations
c           ---------------
c
      uflow = d1mach(1)
      oflow = d1mach(2)
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
      nrmax = 1
      extrap = .false.
      noext = .false.
      ierro = 0
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ktmin = 0
      small = dabs(b-a)*0.75d+00
      nres = 0
      numrl2 = 0
      extall = .false.
      if(0.5d+00*dabs(b-a)*domega.gt.0.2d+01) go to 10
      numrl2 = 1
      extall = .true.
      rlist2(1) = result
   10 if(0.25d+00*dabs(b-a)*domega.le.0.2d+01) extall = .true.
      ksgn = -1
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 140 last = 2,limit
c
c           bisect the subinterval with the nrmax-th largest
c           error estimate.
c
        nrmom = nnlog(maxerr)+1
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dqc25f(f,a1,b1,domega,integr,nrmom,maxp1,0,
     *  area1,error1,nev,resabs,defab1,momcom,chebmo)
        neval = neval+nev
        call dqc25f(f,a2,b2,domega,integr,nrmom,maxp1,1,
     *  area2,error2,nev,resabs,defab2,momcom,chebmo)
        neval = neval+nev
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 25
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     *  .or.erro12.lt.0.99d+00*errmax) go to 20
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   20   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   25   rlist(maxerr) = area1
        rlist(last) = area2
        nnlog(maxerr) = nrmom
        nnlog(last) = nrmom
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
c           at a point of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)
     *  *(dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 30
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 40
   30   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to bisected next).
c
   40   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
      if(errsum.le.errbnd) go to 170
      if(ier.ne.0) go to 150
        if(last.eq.2.and.extall) go to 120
        if(noext) go to 140
        if(.not.extall) go to 50
        erlarg = erlarg-erlast
        if(dabs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 70
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
   50   width = dabs(blist(maxerr)-alist(maxerr))
        if(width.gt.small) go to 140
        if(extall) go to 60
c
c           test whether we can start with the extrapolation procedure
c           (we do this if we integrate over the next interval with
c           use of a gauss-kronrod rule - see subroutine dqc25f).
c
        small = small*0.5d+00
        if(0.25d+00*width*domega.gt.0.2d+01) go to 140
        extall = .true.
        go to 130
   60   extrap = .true.
        nrmax = 2
   70   if(ierro.eq.3.or.erlarg.le.ertest) go to 90
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over
c           the larger intervals (erlarg) and perform extrapolation.
c
        jupbnd = last
        if (last.gt.(limit/2+2)) jupbnd = limit+3-last
        id = nrmax
        do 80 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 140
          nrmax = nrmax+1
   80   continue
c
c           perform extrapolation.
c
   90   numrl2 = numrl2+1
        rlist2(numrl2) = area
        if(numrl2.lt.3) go to 110
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 100
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
c ***jump out of do-loop
        if(abserr.le.ertest) go to 150
c
c           prepare bisection of the smallest interval.
c
  100   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 150
  110   maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5d+00
        erlarg = errsum
        go to 140
  120   small = small*0.5d+00
        numrl2 = numrl2+1
        rlist2(numrl2) = area
  130   ertest = errbnd
        erlarg = errsum
  140 continue
c
c           set the final result.
c           ---------------------
c
  150 if(abserr.eq.oflow.or.nres.eq.0) go to 170
      if(ier+ierro.eq.0) go to 165
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00) go to 160
      if(abserr.gt.errsum) go to 170
      if(area.eq.0.0d+00) go to 190
      go to 165
  160 if(abserr/dabs(result).gt.errsum/dabs(area)) go to 170
c
c           test on divergence.
c
  165 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     * defabs*0.1d-01) go to 190
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03
     * .or.errsum.ge.dabs(area)) ier = 6
      go to 190
c
c           compute global integral sum.
c
  170 result = 0.0d+00
      do 180 k=1,last
        result = result+rlist(k)
  180 continue
      abserr = errsum
  190 if (ier.gt.2) ier=ier-1
  200 if (integr.eq.2.and.omega.lt.0.0d+00) result=-result
  999 return
      end
