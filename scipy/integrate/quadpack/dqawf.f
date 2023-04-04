      recursive subroutine dqawf(f,a,omega,integr,epsabs,result,
     *   abserr,neval,ier,limlst,lst,leniw,maxp1,lenw,iwork,work)
c***begin prologue  dqawf
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1
c***keywords  automatic integrator, special-purpose,fourier
c             integral, integration between zeros with dqawoe,
c             convergence acceleration with dqelg
c***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            fourier integral i=integral of f(x)*w(x) over (a,infinity)
c            where w(x) = cos(omega*x) or w(x) = sin(omega*x).
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.epsabs.
c***description
c
c        computation of fourier integrals
c        standard fortran subroutine
c        double precision version
c
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
c            omega  - double precision
c                     parameter in the integrand weight function
c
c            integr - integer
c                     indicates which of the weight functions is used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1.and.integr.ne.2, the routine
c                     will end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested, epsabs.gt.0.
c                     if epsabs.le.0, the routine will end with ier = 6.
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
c                    if omega.ne.0
c                     ier = 1 maximum number of cycles allowed
c                             has been achieved, i.e. of subintervals
c                             (a+(k-1)c,a+kc) where
c                             c = (2*int(abs(omega))+1)*pi/abs(omega),
c                             for k = 1, 2, ..., lst.
c                             one can allow more cycles by increasing
c                             the value of limlst (and taking the
c                             according dimension adjustments into
c                             account). examine the array iwork which
c                             contains the error flags on the cycles, in
c                             order to look for eventual local
c                             integration difficulties.
c                             if the position of a local difficulty
c                             can be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling
c                             appropriate integrators on the subranges.
c                         = 4 the extrapolation table constructed for
c                             convergence accelaration of the series
c                             formed by the integral contributions over
c                             the cycles, does not converge to within
c                             the requested accuracy.
c                             as in the case of ier = 1, it is advised
c                             to examine the array iwork which contains
c                             the error flags on the cycles.
c                         = 6 the input is invalid because
c                             (integr.ne.1 and integr.ne.2) or
c                              epsabs.le.0 or limlst.lt.1 or
c                              leniw.lt.(limlst+2) or maxp1.lt.1 or
c                              lenw.lt.(leniw*2+maxp1*25).
c                              result, abserr, neval, lst are set to
c                              zero.
c                         = 7 bad integrand behaviour occurs within
c                             one or more of the cycles. location and
c                             type of the difficulty involved can be
c                             determined from the first lst elements of
c                             vector iwork.  here lst is the number of
c                             cycles actually needed (see below).
c                             iwork(k) = 1 the maximum number of
c                                          subdivisions (=(leniw-limlst)
c                                          /2) has been achieved on the
c                                          k th cycle.
c                                      = 2 occurrence of roundoff error
c                                          is detected and prevents the
c                                          tolerance imposed on the k th
c                                          cycle, from being achieved
c                                          on this cycle.
c                                      = 3 extremely bad integrand
c                                          behaviour occurs at some
c                                          points of the k th cycle.
c                                      = 4 the integration procedure
c                                          over the k th cycle does
c                                          not converge (to within the
c                                          required accuracy) due to
c                                          roundoff in the extrapolation
c                                          procedure invoked on this
c                                          cycle. it is assumed that the
c                                          result on this interval is
c                                          the best which can be
c                                          obtained.
c                                      = 5 the integral over the k th
c                                          cycle is probably divergent
c                                          or slowly convergent. it must
c                                          be noted that divergence can
c                                          occur with any other value of
c                                          iwork(k).
c                    if omega = 0 and integr = 1,
c                    the integral is calculated by means of dqagie,
c                    and ier = iwork(1) (with meaning as described
c                    for iwork(k),k = 1).
c
c         dimensioning parameters
c            limlst - integer
c                     limlst gives an upper bound on the number of
c                     cycles, limlst.ge.3.
c                     if limlst.lt.3, the routine will end with ier = 6.
c
c            lst    - integer
c                     on return, lst indicates the number of cycles
c                     actually needed for the integration.
c                     if omega = 0, then lst is set to 1.
c
c            leniw  - integer
c                     dimensioning parameter for iwork. on entry,
c                     (leniw-limlst)/2 equals the maximum number of
c                     subintervals allowed in the partition of each
c                     cycle, leniw.ge.(limlst+2).
c                     if leniw.lt.(limlst+2), the routine will end with
c                     ier = 6.
c
c            maxp1  - integer
c                     maxp1 gives an upper bound on the number of
c                     chebyshev moments which can be stored, i.e. for
c                     the intervals of lengths abs(b-a)*2**(-l),
c                     l = 0,1, ..., maxp1-2, maxp1.ge.1.
c                     if maxp1.lt.1, the routine will end with ier = 6.
c            lenw   - integer
c                     dimensioning parameter for work
c                     lenw must be at least leniw*2+maxp1*25.
c                     if lenw.lt.(leniw*2+maxp1*25), the routine will
c                     end with ier = 6.
c
c         work arrays
c            iwork  - integer
c                     vector of dimension at least leniw
c                     on return, iwork(k) for k = 1, 2, ..., lst
c                     contain the error flags on the cycles.
c
c            work   - double precision
c                     vector of dimension at least
c                     on return,
c                     work(1), ..., work(lst) contain the integral
c                      approximations over the cycles,
c                     work(limlst+1), ..., work(limlst+lst) contain
c                      the error extimates over the cycles.
c                     further elements of work have no specific
c                     meaning for the user.
c
c***references  (none)
c***routines called  dqawfe,xerror
c***end prologue  dqawf
c
       double precision a,abserr,epsabs,f,omega,result,work
       integer ier,integr,iwork,last,leniw,lenw,limit,limlst,ll2,lvl,
     *  lst,l1,l2,l3,l4,l5,l6,maxp1,neval
c
       dimension iwork(leniw),work(lenw)
c
       external f
c
c         check validity of limlst, leniw, maxp1 and lenw.
c
c***first executable statement  dqawf
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limlst.lt.3.or.leniw.lt.(limlst+2).or.maxp1.lt.1.or.lenw.lt.
     *   (leniw*2+maxp1*25)) go to 10
c
c         prepare call for dqawfe
c
      limit = (leniw-limlst)/2
      l1 = limlst+1
      l2 = limlst+l1
      l3 = limit+l2
      l4 = limit+l3
      l5 = limit+l4
      l6 = limit+l5
      ll2 = limit+l1
      call dqawfe(f,a,omega,integr,epsabs,limlst,limit,maxp1,result,
     *  abserr,neval,ier,work(1),work(l1),iwork(1),lst,work(l2),
     *  work(l3),work(l4),work(l5),iwork(l1),iwork(ll2),work(l6))
c
c         call error handler if necessary
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror('abnormal return from dqawf',26,ier,lvl)
      return
      end
