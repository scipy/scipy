      recursive subroutine dqawo(f,a,b,omega,integr,epsabs,epsrel,
     *   result,abserr,neval,ier,leniw,maxp1,lenw,last,iwork,work)
c***begin prologue  dqawo
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
c            definite integral i=integral of f(x)*w(x) over (a,b)
c            where w(x) = cos(omega*x)
c            or w(x) = sin(omega*x),
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
c                     function subprogram defining the function
c                     f(x).  the actual name for f needs to be
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
c                     indicates which of the weight functions is used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1.and.integr.ne.2, the routine will
c                     end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if epsabs.le.0 and
c                     epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
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
c                     number of  integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             (= leniw/2) has been achieved. one can
c                             allow more subdivisions by increasing the
c                             value of leniw (and taking the according
c                             dimension adjustments into account).
c                             however, if this yields no improvement it
c                             is advised to analyze the integrand in
c                             order to determine the integration
c                             difficulties. if the position of a local
c                             difficulty can be determined (e.g.
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
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some interior points of the
c                             integration interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table. it is presumed that
c                             the requested tolerance cannot be achieved
c                             due to roundoff in the extrapolation
c                             table, and that the returned result is
c                             the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or (integr.ne.1 and integr.ne.2),
c                             or leniw.lt.2 or maxp1.lt.1 or
c                             lenw.lt.leniw*2+maxp1*25.
c                             result, abserr, neval, last are set to
c                             zero. except when leniw, maxp1 or lenw are
c                             invalid, work(limit*2+1), work(limit*3+1),
c                             iwork(1), iwork(limit+1) are set to zero,
c                             work(1) is set to a and work(limit+1) to
c                             b.
c
c         dimensioning parameters
c            leniw  - integer
c                     dimensioning parameter for iwork.
c                     leniw/2 equals the maximum number of subintervals
c                     allowed in the partition of the given integration
c                     interval (a,b), leniw.ge.2.
c                     if leniw.lt.2, the routine will end with ier = 6.
c
c            maxp1  - integer
c                     gives an upper bound on the number of chebyshev
c                     moments which can be stored, i.e. for the
c                     intervals of lengths abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1
c                     if maxp1.lt.1, the routine will end with ier = 6.
c
c            lenw   - integer
c                     dimensioning parameter for work
c                     lenw must be at least leniw*2+maxp1*25.
c                     if lenw.lt.(leniw*2+maxp1*25), the routine will
c                     end with ier = 6.
c
c            last   - integer
c                     on return, last equals the number of subintervals
c                     produced in the subdivision process, which
c                     determines the number of significant elements
c                     actually in the work arrays.
c
c         work arrays
c            iwork  - integer
c                     vector of dimension at least leniw
c                     on return, the first k elements of which contain
c                     pointers to the error estimates over the
c                     subintervals, such that work(limit*3+iwork(1)), ..
c                     work(limit*3+iwork(k)) form a decreasing
c                     sequence, with limit = lenw/2 , and k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise.
c                     furthermore, iwork(limit+1), ..., iwork(limit+
c                     last) indicate the subdivision levels of the
c                     subintervals, such that iwork(limit+i) = l means
c                     that the subinterval numbered i is of length
c                     abs(b-a)*2**(1-l).
c
c            work   - double precision
c                     vector of dimension at least lenw
c                     on return
c                     work(1), ..., work(last) contain the left
c                      end points of the subintervals in the
c                      partition of (a,b),
c                     work(limit+1), ..., work(limit+last) contain
c                      the right end points,
c                     work(limit*2+1), ..., work(limit*2+last) contain
c                      the integral approximations over the
c                      subintervals,
c                     work(limit*3+1), ..., work(limit*3+last)
c                      contain the error estimates.
c                     work(limit*4+1), ..., work(limit*4+maxp1*25)
c                      provide space for storing the chebyshev moments.
c                     note that limit = lenw/2.
c
c***references  (none)
c***routines called  dqawoe,xerror
c***end prologue  dqawo
c
       double precision a,abserr,b,epsabs,epsrel,f,omega,result,work
       integer ier,integr,iwork,last,limit,lenw,leniw,lvl,l1,l2,l3,l4,
     *  maxp1,momcom,neval
c
       dimension iwork(leniw),work(lenw)
c
       external f
c
c         check validity of leniw, maxp1 and lenw.
c
c***first executable statement  dqawo
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(leniw.lt.2.or.maxp1.lt.1.or.lenw.lt.(leniw*2+maxp1*25))
     *  go to 10
c
c         prepare call for dqawoe
c
      limit = leniw/2
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
      l4 = limit+l3
      call dqawoe(f,a,b,omega,integr,epsabs,epsrel,limit,1,maxp1,result,
     *   abserr,neval,ier,last,work(1),work(l1),work(l2),work(l3),
     *   iwork(1),iwork(l1),momcom,work(l4))
c
c         call error handler if necessary
c
      lvl = 0
10    if(ier.eq.6) lvl = 0
      if(ier.ne.0) call xerror('abnormal return from dqawo',26,ier,lvl)
      return
      end
