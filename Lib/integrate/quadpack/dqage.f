      subroutine dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr,
     *   neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqage
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             integrand examinator, globally adaptive,
c             gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral   i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
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
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            key    - integer
c                     key for choice of local integration rule
c                     a gauss-kronrod pair is used with
c                          7 - 15 points if key.lt.2,
c                         10 - 21 points if key = 2,
c                         15 - 31 points if key = 3,
c                         20 - 41 points if key = 4,
c                         25 - 51 points if key = 5,
c                         30 - 61 points if key.gt.5.
c
c            limit  - integer
c                     gives an upperbound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1.
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
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for result and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value
c                             of limit.
c                             however, if this yields no improvement it
c                             is rather advised to analyze the integrand
c                             in order to determine the integration
c                             difficulties. if the position of a local
c                             difficulty can be determined(e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling the integrator on the
c                             subranges. if possible, an appropriate
c                             special-purpose integrator should be used
c                             which is designed for handling the type of
c                             difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             result, abserr, neval, last, rlist(1) ,
c                             elist(1) and iord(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c
c            alist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the left
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            blist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the right
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            rlist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the
c                      integral approximations on the subintervals
c
c            elist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the moduli of the
c                      absolute error estimates on the subintervals
c
c            iord    - integer
c                      vector of dimension at least limit, the first k
c                      elements of which are pointers to the
c                      error estimates over the subintervals,
c                      such that elist(iord(1)), ...,
c                      elist(iord(k)) form a decreasing sequence,
c                      with k = last if last.le.(limit/2+2), and
c                      k = limit+1-last otherwise
c
c            last    - integer
c                      number of subintervals actually produced in the
c                      subdivision process
c
c***references  (none)
c***routines called  d1mach,dqk15,dqk21,dqk31,
c                    dqk41,dqk51,dqk61,dqpsrt
c***end prologue  dqage
c
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  rlist(limit)
c
      external f
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                      (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach  is the largest relative spacing.
c           uflow  is the smallest positive magnitude.
c
c***first executable statement  dqage
      epmach = d1mach(4)
      uflow = d1mach(1)
c
c           test on validity of parameters
c           ------------------------------
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
      if(epsabs.le.0.0d+00.and.
     *  epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28)) ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      keyf = key
      if(key.le.0) keyf = 1
      if(key.ge.7) keyf = 6
      neval = 0
      if(keyf.eq.1) call dqk15(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.2) call dqk21(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.3) call dqk31(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.4) call dqk41(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.5) call dqk51(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.6) call dqk61(f,a,b,result,abserr,defabs,resabs)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
c
c           test on accuracy.
c
      errbnd = dmax1(epsabs,epsrel*dabs(result))
      if(abserr.le.0.5d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs)
     *  .or.abserr.eq.0.0d+00) go to 60
c
c           initialization
c           --------------
c
c
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
c
c           main do-loop
c           ------------
c
      do 30 last = 2,limit
c
c           bisect the subinterval with the largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        if(keyf.eq.1) call dqk15(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.2) call dqk21(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.3) call dqk31(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.4) call dqk41(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.5) call dqk51(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.6) call dqk61(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.1) call dqk15(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.2) call dqk21(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.3) call dqk31(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.4) call dqk41(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.5) call dqk51(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.6) call dqk61(f,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        neval = neval+1
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 5
        if(dabs(rlist(maxerr)-area12).le.0.1d-04*dabs(area12)
     *  .and.erro12.ge.0.99d+00*errmax) iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
    5   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
        if(errsum.le.errbnd) go to 8
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
c
c           set error flag in the case that the number of subintervals
c           equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*
     *  epmach)*(dabs(a2)+0.1d+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
    8   if(error2.gt.error1) go to 10
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 20
   10   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with the largest error estimate (to be bisected next).
c
   20   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(ier.ne.0.or.errsum.le.errbnd) go to 40
   30 continue
c
c           compute final result.
c           ---------------------
c
   40 result = 0.0d+00
      do 50 k=1,last
        result = result+rlist(k)
   50 continue
      abserr = errsum
   60 if(keyf.ne.1) neval = (10*keyf+1)*(2*neval+1)
      if(keyf.eq.1) neval = 30*neval+15
  999 return
      end
