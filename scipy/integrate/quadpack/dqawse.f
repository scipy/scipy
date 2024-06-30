      recursive subroutine dqawse(f,a,b,alfa,beta,integr,epsabs,epsrel,
     *   limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,
     *   last)
c***begin prologue  dqawse
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             algebraico-logarithmic end point singularities,
c             clenshaw-curtis method
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f*w over (a,b),
c            (where w shows a singular behaviour at the end points,
c            see parameter integr).
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        integration of functions having algebraico-logarithmic
c        end point singularities
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
c                     upper limit of integration, b.gt.a
c                     if b.le.a, the routine will end with ier = 6.
c
c            alfa   - double precision
c                     parameter in the weight function, alfa.gt.(-1)
c                     if alfa.le.(-1), the routine will end with
c                     ier = 6.
c
c            beta   - double precision
c                     parameter in the weight function, beta.gt.(-1)
c                     if beta.le.(-1), the routine will end with
c                     ier = 6.
c
c            integr - integer
c                     indicates which weight function is to be used
c                     = 1  (x-a)**alfa*(b-x)**beta
c                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
c                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
c                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
c                     if integr.lt.1 or integr.gt.4, the routine
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
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.2
c                     if limit.lt.2, the routine will end with ier = 6.
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
c                             the estimates for the integral and error
c                             are less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                         = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit. however, if this yields no
c                             improvement, it is advised to analyze the
c                             integrand in order to determine the
c                             integration difficulties which prevent the
c                             requested tolerance from being achieved.
c                             in case of a jump discontinuity or a local
c                             singularity of algebraico-logarithmic type
c                             at one or more interior points of the
c                             integration range, one should proceed by
c                             splitting up the interval at these
c                             points and calling the integrator on the
c                             subranges.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             b.le.a or alfa.le.(-1) or beta.le.(-1), or
c                             integr.lt.1 or integr.gt.4, or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             or limit.lt.2.
c                             result, abserr, neval, rlist(1), elist(1),
c                             iord(1) and last are set to zero. alist(1)
c                             and blist(1) are set to a and b
c                             respectively.
c
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
c                     vector of dimension at least limit,the first
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
c                     of which are pointers to the error
c                     estimates over the subintervals, so that
c                     elist(iord(1)), ..., elist(iord(k)) with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise form a decreasing sequence
c
c            last   - integer
c                     number of subintervals actually produced in
c                     the subdivision process
c
c***references  (none)
c***routines called  d1mach,dqc25s,dqmomo,dqpsrt
c***end prologue  dqawse
c
      double precision a,abserr,alfa,alist,area,area1,area12,area2,a1,
     *  a2,b,beta,blist,b1,b2,centre,dabs,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,erro12,error2,errsum,f,
     *  resas1,resas2,result,rg,rh,ri,rj,rlist,uflow
      integer ier,integr,iord,iroff1,iroff2,k,last,limit,maxerr,nev,
     *  neval,nrmax
c
      external f
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit),ri(25),rj(25),rh(25),rg(25)
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
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqawse
      epmach = d1mach(4)
      uflow = d1mach(1)
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 6
      neval = 0
      last = 0
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(b.le.a.or.(epsabs.eq.0.0d+00.and.
     *  epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28)).or.alfa.le.(-0.1d+01).
     *  or.beta.le.(-0.1d+01).or.integr.lt.1.or.integr.gt.4.or.
     *  limit.lt.2) go to 999
      ier = 0
c
c           compute the modified chebyshev moments.
c
      call dqmomo(alfa,beta,ri,rj,rg,rh,integr)
c
c           integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
c
      centre = 0.5d+00*(b+a)
      call dqc25s(f,a,b,a,centre,alfa,beta,ri,rj,rg,rh,area1,
     *  error1,resas1,integr,nev)
      neval = nev
      call dqc25s(f,a,b,centre,b,alfa,beta,ri,rj,rg,rh,area2,
     *  error2,resas2,integr,nev)
      last = 2
      neval = neval+nev
      result = area1+area2
      abserr = error1+error2
c
c           test on accuracy.
c
      errbnd = dmax1(epsabs,epsrel*dabs(result))
c
c           initialization
c           --------------
c
      if(error2.gt.error1) go to 10
      alist(1) = a
      alist(2) = centre
      blist(1) = centre
      blist(2) = b
      rlist(1) = area1
      rlist(2) = area2
      elist(1) = error1
      elist(2) = error2
      go to 20
   10 alist(1) = centre
      alist(2) = a
      blist(1) = b
      blist(2) = centre
      rlist(1) = area2
      rlist(2) = area1
      elist(1) = error2
      elist(2) = error1
   20 iord(1) = 1
      iord(2) = 2
      if(limit.eq.2) ier = 1
      if(abserr.le.errbnd.or.ier.eq.1) go to 999
      errmax = elist(1)
      maxerr = 1
      nrmax = 1
      area = result
      errsum = abserr
      iroff1 = 0
      iroff2 = 0
c
c            main do-loop
c            ------------
c
      do 60 last = 3,limit
c
c           bisect the subinterval with largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
c
        call dqc25s(f,a,b,a1,b1,alfa,beta,ri,rj,rg,rh,area1,
     *  error1,resas1,integr,nev)
        neval = neval+nev
        call dqc25s(f,a,b,a2,b2,alfa,beta,ri,rj,rg,rh,area2,
     *  error2,resas2,integr,nev)
        neval = neval+nev
c
c           improve previous approximations integral and error
c           and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(a.eq.a1.or.b.eq.b2) go to 30
        if(resas1.eq.error1.or.resas2.eq.error2) go to 30
c
c           test for roundoff error.
c
        if(dabs(rlist(maxerr)-area12).lt.0.1d-04*dabs(area12)
     *  .and.erro12.ge.0.99d+00*errmax) iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
   30   rlist(maxerr) = area1
        rlist(last) = area2
c
c           test on accuracy.
c
        errbnd = dmax1(epsabs,epsrel*dabs(area))
        if(errsum.le.errbnd) go to 35
c
c           set error flag in the case that the number of interval
c           bisections exceeds limit.
c
        if(last.eq.limit) ier = 1
c
c
c           set error flag in the case of roundoff error.
c
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
c
c           set error flag in the case of bad integrand behaviour
c           at interior points of integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     *  (dabs(a2)+0.1d+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
   35   if(error2.gt.error1) go to 40
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 50
   40   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with largest error estimate (to be bisected next).
c
   50   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if (ier.ne.0.or.errsum.le.errbnd) go to 70
   60 continue
c
c           compute final result.
c           ---------------------
c
   70 result = 0.0d+00
      do 80 k=1,last
        result = result+rlist(k)
   80 continue
      abserr = errsum
  999 return
      end
