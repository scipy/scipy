      subroutine dqaws(f,a,b,alfa,beta,integr,epsabs,epsrel,result,
     *   abserr,neval,ier,limit,lenw,last,iwork,work)
c***begin prologue  dqaws
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             algebraico-logarithmic end-point singularities,
c             clenshaw-curtis, globally adaptive
c***author  piessens,robert,appl. math. & progr. div. -k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f*w over (a,b),
c            (where w shows a singular behaviour at the end points
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
c                     parameter in the integrand function, alfa.gt.(-1)
c                     if alfa.le.(-1), the routine will end with
c                     ier = 6.
c
c            beta   - double precision
c                     parameter in the integrand function, beta.gt.(-1)
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
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand, in order to
c                             determine the integration difficulties
c                             which prevent the requested tolerance from
c                             being achieved. in case of a jump
c                             discontinuity or a local singularity
c                             of algebraico-logarithmic type at one or
c                             more interior points of the integration
c                             range, one should proceed by splitting up
c                             the interval at these points and calling
c                             the integrator on the subranges.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             b.le.a or alfa.le.(-1) or beta.le.(-1) or
c                             or integr.lt.1 or integr.gt.4 or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.2 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero. except when lenw or limit is invalid
c                             iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit  - integer
c                     dimensioning parameter for iwork
c                     limit determines the maximum number of
c                     subintervals in the partition of the given
c                     integration interval (a,b), limit.ge.2.
c                     if limit.lt.2, the routine will end with ier = 6.
c
c            lenw   - integer
c                     dimensioning parameter for work
c                     lenw must be at least limit*4.
c                     if lenw.lt.limit*4, the routine will end
c                     with ier = 6.
c
c            last   - integer
c                     on return, last equals the number of
c                     subintervals produced in the subdivision process,
c                     which determines the significant number of
c                     elements actually in the work arrays.
c
c         work arrays
c            iwork  - integer
c                     vector of dimension limit, the first k
c                     elements of which contain pointers
c                     to the error estimates over the subintervals,
c                     such that work(limit*3+iwork(1)), ...,
c                     work(limit*3+iwork(k)) form a decreasing
c                     sequence with k = last if last.le.(limit/2+2),
c                     and k = limit+1-last otherwise
c
c            work   - double precision
c                     vector of dimension lenw
c                     on return
c                     work(1), ..., work(last) contain the left
c                      end points of the subintervals in the
c                      partition of (a,b),
c                     work(limit+1), ..., work(limit+last) contain
c                      the right end points,
c                     work(limit*2+1), ..., work(limit*2+last)
c                      contain the integral approximations over
c                      the subintervals,
c                     work(limit*3+1), ..., work(limit*3+last)
c                      contain the error estimates.
c
c***references  (none)
c***routines called  dqawse,xerror
c***end prologue  dqaws
c
      double precision a,abserr,alfa,b,beta,epsabs,epsrel,f,result,work
      integer ier,integr,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  dqaws
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limit.lt.2.or.lenw.lt.limit*4) go to 10
c
c         prepare call for dqawse.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call dqawse(f,a,b,alfa,beta,integr,epsabs,epsrel,limit,result,
     *  abserr,neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from dqaws,26,ier,lvl)
      return
      end
