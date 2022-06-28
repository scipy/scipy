      recursive subroutine dqagp(f,a,b,npts2,points,epsabs,epsrel,
     *   result,abserr,neval,ier,leniw,lenw,last,iwork,work)
c***begin prologue  dqagp
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, general-purpose,
c             singularities at user specified points,
c             extrapolation, globally adaptive
c***author  piessens,robert,appl. math. & progr. div - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            break points of the integration interval, where local
c            difficulties of the integrand may occur (e.g.
c            singularities, discontinuities), are provided by the user.
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
c                     range, npts.ge.2.
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
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved, and that
c                             the returned result is the best which
c                             can be obtained.
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
c                             result, abserr, neval, last are set to
c                             zero. except when leniw or lenw or npts2 is
c                             invalid, iwork(1), iwork(limit+1),
c                             work(limit*2+1) and work(limit*3+1)
c                             are set to zero.
c                             work(1) is set to a and work(limit+1)
c                             to b (where limit = (leniw-npts2)/2).
c
c         dimensioning parameters
c            leniw - integer
c                    dimensioning parameter for iwork
c                    leniw determines limit = (leniw-npts2)/2,
c                    which is the maximum number of subintervals in the
c                    partition of the given integration interval (a,b),
c                    leniw.ge.(3*npts2-2).
c                    if leniw.lt.(3*npts2-2), the routine will end with
c                    ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least leniw*2-npts2.
c                    if lenw.lt.leniw*2-npts2, the routine will end
c                    with ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least leniw. on return,
c                    the first k elements of which contain
c                    pointers to the error estimates over the
c                    subintervals, such that work(limit*3+iwork(1)),...,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2), and
c                    k = limit+1-last otherwise
c                    iwork(limit+1), ...,iwork(limit+last) contain the
c                     subdivision levels of the subintervals, i.e.
c                     if (aa,bb) is a subinterval of (p1,p2)
c                     where p1 as well as p2 is a user-provided
c                     break point or integration limit, then (aa,bb) has
c                     level l if abs(bb-aa) = abs(p2-p1)*2**(-l),
c                    iwork(limit*2+1), ..., iwork(limit*2+npts2) have
c                     no significance for the user,
c                    note that limit = (leniw-npts2)/2.
c
c            work  - double precision
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last)
c                     contain the corresponding error estimates,
c                    work(limit*4+1), ..., work(limit*4+npts2)
c                     contain the integration limits and the
c                     break points sorted in an ascending sequence.
c                    note that limit = (leniw-npts2)/2.
c
c***references  (none)
c***routines called  dqagpe,xerror
c***end prologue  dqagp
c
      double precision a,abserr,b,epsabs,epsrel,f,points,result,work
      integer ier,iwork,last,leniw,lenw,limit,lvl,l1,l2,l3,l4,neval,
     *  npts2
c
      dimension iwork(leniw),points(npts2),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  dqagp
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(leniw.lt.(3*npts2-2).or.lenw.lt.(leniw*2-npts2).or.npts2.lt.2)
     *  go to 10
c
c         prepare call for dqagpe.
c
      limit = (leniw-npts2)/2
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
      l4 = limit+l3
c
      call dqagpe(f,a,b,npts2,points,epsabs,epsrel,limit,result,abserr,
     *  neval,ier,work(1),work(l1),work(l2),work(l3),work(l4),
     *  iwork(1),iwork(l1),iwork(l2),last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror('abnormal return from dqagp',26,ier,lvl)
      return
      end
