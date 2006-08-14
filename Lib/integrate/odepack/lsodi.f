      subroutine lsodi (res, adda, jac, neq, y, ydoti, t, tout, itol,
     1  rtol, atol, itask, istate, iopt, rwork, lrw, iwork, liw, mf )
      external res, adda, jac
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
      double precision y, ydoti, t, tout, rtol, atol, rwork
      dimension neq(1), y(1), ydoti(1), rtol(1), atol(1), rwork(lrw),
     1          iwork(liw)
c-----------------------------------------------------------------------
c this is the march 30, 1987 version of lsodi..
c  livermore solver for ordinary differential equations (implicit form).
c this version is in double precision.
c
c lsodi solves the initial value problem for linearly implicit
c systems of first order ode-s,
c     a(t,y) * dy/dt = g(t,y) ,  where a(t,y) is a square matrix,
c or, in component form,
c     ( a   * ( dy / dt ))  + ... +  ( a     * ( dy   / dt ))  =
c        i,1      1                     i,neq      neq
c
c      =   g ( t, y , y ,..., y    )   ( i = 1,...,neq )
c           i      1   2       neq
c
c if a is singular, this is a differential-algebraic system.
c
c lsodi is a variant version of the lsode package.
c-----------------------------------------------------------------------
c reference..
c     alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c-----------------------------------------------------------------------
c authors...  jeffrey f. painter  and
c             alan c. hindmarsh
c             computing and mathematics research division, l-316
c             lawrence livermore national laboratory
c             livermore, ca 94550.
c
c-----------------------------------------------------------------------
c summary of usage.
c
c communication between the user and the lsodi package, for normal
c situations, is summarized here.  this summary describes only a subset
c of the full set of options available.  see the full description for
c details, including optional communication, nonstandard options,
c and instructions for special situations.  see also the example
c problem (with program and output) following this summary.
c
c a. first, provide a subroutine of the form..
c               subroutine res (neq, t, y, s, r, ires)
c               dimension y(neq), s(neq), r(neq)
c which computes the residual function
c     r = g(t,y)  -  a(t,y) * s ,
c as a function of t and the vectors y and s.  (s is an internally
c generated approximation to dy/dt.)  the arrays y and s are inputs
c to the res routine and should not be altered.  the residual
c vector is to be stored in the array r.  the argument ires should be
c ignored for casual use of lsodi.  (for uses of ires, see the
c paragraph on res in the full description below.)
c
c b. next, decide whether full or banded form is more economical
c for the storage of matrices.  lsodi must deal internally with the
c matrices a and dr/dy, where r is the residual function defined above.
c lsodi generates a linear combination of these two matrices, and
c this is treated in either full or banded form.
c     the matrix structure is communicated by a method flag mf,
c which is 21 or 22 for the full case, and 24 or 25 in the band case.
c     in the banded case, lsodi requires two half-bandwidth
c parameters ml and mu.  these are, respectively, the widths of the
c lower and upper parts of the band, excluding the main diagonal.
c thus the band consists of the locations (i,j) with
c i-ml .le. j .le. i+mu, and the full bandwidth is ml+mu+1.
c note that the band must accommodate the nonzero elements of
c a(t,y), dg/dy, and d(a*s)/dy (s fixed).  alternatively, one
c can define a band that encloses only the elements that are relatively
c large in magnitude, and gain some economy in storage and possibly
c also efficiency, although the appropriate threshhold for
c retaining matrix elements is highly problem-dependent.
c
c c. you must also provide a subroutine of the form..
c               subroutine adda (neq, t, y, ml, mu, p, nrowp)
c               dimension y(neq), p(nrowp,neq)
c which adds the matrix a = a(t,y) to the contents of the array p.
c t and the y array are input and should not be altered.
c     in the full matrix case, this routine should add elements of
c to p in the usual order.  i.e., add a(i,j) to p(i,j).  (ignore the
c ml and mu arguments in this case.)
c     in the band matrix case, this routine should add element a(i,j)
c to p(i-j+mu+1,j).  i.e., add the diagonal lines of a to the rows of
c p from the top down (the top line of a added to the first row of p).
c
c d. for the sake of efficiency, you are encouraged to supply the
c jacobian matrix dr/dy in closed form, where r = g(t,y) - a(t,y)*s
c (s = a fixed vector) as above.  if dr/dy is being supplied,
c use mf = 21 or 24, and provide a subroutine of the form..
c               subroutine jac (neq, t, y, s, ml, mu, p, nrowp)
c               dimension y(neq), s(neq), p(nrowp,neq)
c which computes dr/dy as a function of t, y, and s.  here t, y, and
c s are inputs, and the routine is to load dr/dy into p as follows..
c     in the full matrix case (mf = 21), load p(i,j) with dr(i)/dy(j),
c the partial derivative of r(i) with respect to y(j).  (ignore the
c ml and mu arguments in this case.)
c     in the band matrix case (mf = 24), load p(i-j+mu+1,j) with
c dr(i)/dy(j), i.e. load the diagonal lines of dr/dy into the rows of
c p from the top down.
c     in either case, only nonzero elements need be loaded, and the
c indexing of p is the same as in the adda routine.
c     note that if a is independent of y (or this dependence
c is weak enough to be ignored) then jac is to compute dg/dy.
c     if it is not feasible to provide a jac routine, use
c mf = 22 or 25, and lsodi will compute an approximate jacobian
c internally by difference quotients.
c
c e. next decide whether or not to provide the initial value of the
c derivative vector dy/dt.  if the initial value of a(t,y) is
c nonsingular (and not too ill-conditioned), you may let lsodi compute
c this vector (istate = 0).  (lsodi will solve the system a*s = g for
c s, with initial values of a and g.)  if a(t,y) is initially
c singular, then the system is a differential-algebraic system, and
c you must make use of the particular form of the system to compute the
c initial values of y and dy/dt.  in that case, use istate = 1 and
c load the initial value of dy/dt into the array ydoti.
c the input array ydoti and the initial y array must be consistent with
c the equations a*dy/dt = g.  this implies that the initial residual
c r = g(t,y) - a(t,y)*ydoti   must be approximately zero.
c
c f. write a main program which calls subroutine lsodi once for
c each point at which answers are desired.  this should also provide
c for possible use of logical unit 6 for output of error messages
c by lsodi.  on the first call to lsodi, supply arguments as follows..
c res    = name of user subroutine for residual function r.
c adda   = name of user subroutine for computing and adding a(t,y).
c jac    = name of user subroutine for jacobian matrix dr/dy
c          (mf = 21 or 24).  if not used, pass a dummy name.
c note.. the names for the res and adda routines and (if used) the
c        jac routine must be declared external in the calling program.
c neq    = number of scalar equations in the system.
c y      = array of initial values, of length neq.
c ydoti  = array of length neq (containing initial dy/dt if istate = 1).
c t      = the initial value of the independent variable.
c tout   = first point where output is desired (.ne. t).
c itol   = 1 or 2 according as atol (below) is a scalar or array.
c rtol   = relative tolerance parameter (scalar).
c atol   = absolute tolerance parameter (scalar or array).
c          the estimated local error in y(i) will be controlled so as
c          to be roughly less (in magnitude) than
c             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
c             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
c          thus the local error test passes if, in each component,
c          either the absolute error is less than atol (or atol(i)),
c          or the relative error is less than rtol.
c          use rtol = 0.0 for pure absolute error control, and
c          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
c          control.  caution.. actual (global) errors may exceed these
c          local tolerances, so choose them conservatively.
c itask  = 1 for normal computation of output values of y at t = tout.
c istate = integer flag (input and output).  set istate = 1 if the
c          initial dy/dt is supplied, and 0 otherwise.
c iopt   = 0 to indicate no optional inputs used.
c rwork  = real work array of length at least..
c             22 +  9*neq + neq**2           for mf = 21 or 22,
c             22 + 10*neq + (2*ml + mu)*neq  for mf = 24 or 25.
c lrw    = declared length of rwork (in user-s dimension).
c iwork  = integer work array of length at least 20 + neq.
c          if mf = 24 or 25, input in iwork(1),iwork(2) the lower
c          and upper half-bandwidths ml,mu.
c liw    = declared length of iwork (in user-s dimension).
c mf     = method flag.  standard values are..
c          21 for a user-supplied full jacobian.
c          22 for an internally generated full jacobian.
c          24 for a user-supplied banded jacobian.
c          25 for an internally generated banded jacobian.
c          for other choices of mf, see the paragraph on mf in
c          the full description below.
c note that the main program must declare arrays y, ydoti, rwork, iwork,
c and possibly atol.
c
c g. the output from the first call (or any call) is..
c      y = array of computed values of y(t) vector.
c      t = corresponding value of independent variable (normally tout).
c istate = 2  if lsodi was successful, negative otherwise.
c          -1 means excess work done on this call (check all inputs).
c          -2 means excess accuracy requested (tolerances too small).
c          -3 means illegal input detected (see printed message).
c          -4 means repeated error test failures (check all inputs).
c          -5 means repeated convergence failures (perhaps bad jacobian
c             supplied or wrong choice of tolerances).
c          -6 means error weight became zero during problem. (solution
c             component i vanished, and atol or atol(i) = 0.)
c          -7 cannot occur in casual use.
c          -8 means lsodi was unable to compute the initial dy/dt.
c             in casual use, this means a(t,y) is initially singular.
c             supply ydoti and use istate = 1 on the first call.
c
c  if lsodi returns istate = -1, -4, or -5, then the output of
c  lsodi also includes ydoti = array containing residual vector
c  r = g - a * dy/dt  evaluated at the current t, y, and dy/dt.
c
c h. to continue the integration after a successful return, simply
c reset tout and call lsodi again.  no other parameters need be reset.
c
c
c-----------------------------------------------------------------------
c example problem.
c
c the following is a simple example problem, with the coding
c needed for its solution by lsodi.  the problem is from chemical
c kinetics, and consists of the following three equations..
c     dy1/dt = -.04*y1 + 1.e4*y2*y3
c     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
c       0.   = y1 + y2 + y3 - 1.
c on the interval from t = 0.0 to t = 4.e10, with initial conditions
c y1 = 1.0, y2 = y3 = 0.
c
c the following coding solves this problem with lsodi, using mf = 21
c and printing results at t = .4, 4., ..., 4.e10.  it uses
c itol = 2 and atol much smaller for y2 than y1 or y3 because
c y2 has much smaller values.  dy/dt is supplied in ydoti. we had
c obtained the initial value of dy3/dt by differentiating the
c third equation and evaluating the first two at t=0.
c at the end of the run, statistical quantities of interest are
c printed (see optional outputs in the full description below).
c
c     external resid, aplusp, dgbydy
c     double precision atol, rtol, rwork, t, tout, y, ydoti
c     dimension y(3), ydoti(3), atol(3), rwork(58), iwork(23)
c     neq = 3
c     y(1) = 1.d0
c     y(2) = 0.d0
c     y(3) = 0.d0
c     ydoti(1) = -.04d0
c     ydoti(2) =  .04d0
c     ydoti(3) =  0.d0
c     t = 0.d0
c     tout = .4d0
c     itol = 2
c     rtol = 1.d-4
c     atol(1) = 1.d-6
c     atol(2) = 1.d-10
c     atol(3) = 1.d-6
c     itask = 1
c     istate = 1
c     iopt = 0
c     lrw = 58
c     liw = 23
c     mf = 21
c     do 40  iout = 1,12
c        call lsodi(resid, aplusp, dgbydy, neq, y, ydoti, t, tout, itol,
c    1      rtol, atol, itask, istate, iopt, rwork, lrw, iwork, liw, mf)
c        write (6,20)  t, y(1), y(2), y(3)
c  20    format(' at t =',e12.4,'   y =',3e14.6)
c        if (istate .lt. 0 )  go to 80
c  40    tout = tout*10.d0
c     write (6,60)  iwork(11), iwork(12), iwork(13)
c  60 format(/' no. steps =',i4,'  no. r-s =',i4,
c    1         '  no. j-s =',i4)
c     stop
c  80 write (6,90)  istate
c  90 format(///' error halt.. istate =',i3)
c     stop
c     end
c
c     subroutine resid(neq, t, y, s, r, ires)
c     double precision r, s, t, y
c     dimension y(3), s(3), r(3)
c     r(1) = -.04d0*y(1) + 1.d4*y(2)*y(3) - s(1)
c     r(2) = .04d0*y(1) - 1.d4*y(2)*y(3) - 3.d7*y(2)*y(2) - s(2)
c     r(3) = y(1) + y(2) + y(3) - 1.d0
c     return
c     end
c
c     subroutine aplusp(neq, t, y, ml, mu, p, nrowp)
c     double precision p, t, y
c     dimension y(3), p(nrowp,3)
c     p(1,1) = p(1,1) + 1.d0
c     p(2,2) = p(2,2) + 1.d0
c     return
c     end
c
c     subroutine dgbydy(neq, t, y, s, ml, mu, p, nrowp)
c     double precision s, t, p, y
c     dimension y(3), s(3), p(nrowp,3)
c     p(1,1) = -.04d0
c     p(1,2) = 1.d4*y(3)
c     p(1,3) = 1.d4*y(2)
c     p(2,1) = .04d0
c     p(2,2) = -1.d4*y(3) - 6.d7*y(2)
c     p(2,3) = -1.d4*y(2)
c     p(3,1) = 1.d0
c     p(3,2) = 1.d0
c     p(3,3) = 1.d0
c     return
c     end
c
c the output of this program (on a cdc-7600 in single precision)
c is as follows..
c
c   at t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02
c   at t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02
c   at t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01
c   at t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01
c   at t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01
c   at t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01
c   at t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01
c   at t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01
c   at t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01
c   at t =  4.0000e+08   y =  5.494532e-06  2.197826e-11  9.999945e-01
c   at t =  4.0000e+09   y =  5.129457e-07  2.051784e-12  9.999995e-01
c   at t =  4.0000e+10   y = -7.170472e-08 -2.868188e-13  1.000000e+00
c
c   no. steps = 330  no. r-s = 404  no. j-s =  69
c-----------------------------------------------------------------------
c full description of user interface to lsodi.
c
c the user interface to lsodi consists of the following parts.
c
c i.   the call sequence to subroutine lsodi, which is a driver
c      routine for the solver.  this includes descriptions of both
c      the call sequence arguments and of user-supplied routines.
c      following these descriptions is a description of
c      optional inputs available through the call sequence, and then
c      a description of optional outputs (in the work arrays).
c
c ii.  descriptions of other routines in the lsodi package that may be
c      (optionally) called by the user.  these provide the ability to
c      alter error message handling, save and restore the internal
c      common, and obtain specified derivatives of the solution y(t).
c
c iii. descriptions of common blocks to be declared in overlay
c      or similar environments, or to be saved when doing an interrupt
c      of the problem and continued solution later.
c
c iv.  description of two routines in the lsodi package, either of
c      which the user may replace with his own version, if desired.
c      these relate to the measurement of errors.
c
c-----------------------------------------------------------------------
c part i.  call sequence.
c
c the call sequence parameters used for input only are
c     res, adda, jac, neq, tout, itol, rtol, atol, itask,
c     iopt, lrw, liw, mf,
c and those used for both input and output are
c     y, t, istate, ydoti.
c the work arrays rwork and iwork are also used for conditional and
c optional inputs and optional outputs.  (the term output here refers
c to the return from subroutine lsodi to the user-s calling program.)
c
c the legality of input parameters will be thoroughly checked on the
c initial call for the problem, but not checked thereafter unless a
c change in input parameters is flagged by istate = 3 on input.
c
c the descriptions of the call arguments are as follows.
c
c res    = the name of the user-supplied subroutine which supplies
c          the residual vector for the ode system, defined by
c            r = g(t,y) - a(t,y) * s
c          as a function of the scalar t and the vectors
c          s and y ( s approximates dy/dt ). this
c          subroutine is to have the form
c              subroutine res ( neq, t, y, s, r, ires )
c              dimension y(1), s(1), r(1)
c          where neq, t, y, s, and ires are input, and r and
c          ires are output. y, s, and r are arrays of length neq.
c          in dimension statements such as that above, 1 is a
c          dummy dimension. it can be replaced by any value.
c             on input, ires indicates how lsodi will use the
c          returned array r, as follows..
c             ires = 1  means that lsodi needs the full residual,
c                       r = g - a*s, exactly.
c             ires = -1 means that lsodi is using r only to compute
c                       the jacobian dr/dy by difference quotients.
c          the res routine can ignore ires, or it can omit some terms
c          if ires = -1.  if a does not depend on y, then res can
c          just return r = g when ires = -1.  if g - a*s contains other
c          additive terms that are independent of y, these can also be
c          dropped, if done consistently, when ires = -1.
c             the subroutine should set the flag ires if it
c          encounters a halt condition or illegal input.
c          otherwise, it should not reset ires.  on output,
c             ires = 1 or -1 represents a normal return, and
c          lsodi continues integrating the ode.  leave ires
c          unchanged from its input value.
c             ires = 2 tells lsodi to immediately return control
c          to the calling program, with istate = 3.  this lets
c          the calling program change parameters of the prob-
c          lem if necessary.
c             ires = 3 represents an error condition (for example, an
c          illegal value of y). lsodi tries to integrate the ode without
c          getting ires = 3 from res.  if it cannot, lsodi returns
c          with istate = -7 or -1.
c             on an lsodi return with istate = 3, -1, or -7, the values
c          of t and y returned correspond to the last point reached
c          successfully without getting the flag ires = 2 or 3.
c             the flag values ires = 2 and 3 should not be used to
c          handle switches or root-stop conditions.  this is better
c          done by calling lsodi in a one-step mode and checking the
c          stopping function for a sign change at each step.
c             if quantities computed in the res routine are needed
c          externally to lsodi, an extra call to res should be made
c          for this purpose, for consistent and accurate results.
c          to get the current dy/dt for the s argument, use intdy.
c             res must be declared external in the calling
c          program.  see note below for more about res.
c
c adda   = the name of the user-supplied subroutine which adds
c          the matrix a = a(t,y) to another matrix stored in the same
c          form as a. the storage form is determined by miter (see
c          mf).  this subroutine is to have the form
c               subroutine adda ( neq, t, y, ml, mu, p, nrowp )
c               dimension y(1), p(nrowp,1)
c          where neq, t, y, ml, mu, and nrowp are input and p is
c          output. y is an array of length neq, and the matrix p is
c          stored in an nrowp by neq array.
c             in the full matrix case ( miter =  1 or 2 ) adda should
c          add  a    to p(i,j). ml and mu are ignored.
c                i,j
c             in the band matrix case ( miter = 4 or 5 ) adda should
c          add  a    to  p(i-j+mu+1,j).
c                i,j
c          see jac for details on this band storage form.
c             adda must be declared external in the calling program.
c          see note below for more information about adda.
c
c jac    = the name of the user-supplied subroutine which supplies
c          the jacobian matrix, dr/dy, where r = g-a*s. the form of the
c          jacobian matrix is determined by miter. jac is required
c          if miter = 1 or 4 -- otherwise a dummy name can be
c          passed. this subroutine is to have the form
c               subroutine jac ( neq, t, y, s, ml, mu, p, nrowp )
c               dimension y(1), s(1), p(nrowp,1)
c          where neq, t, y, s, ml, mu, and nrowp are input and p
c          is output. y and s are arrays of length neq, and the
c          matrix p is stored in an nrowp by neq array.
c          p is to be loaded with partial derivatives ( elements
c          of the jacobian matrix ) on output.
c             in the full matrix case ( miter = 1 ), ml and mu
c          are ignored and the jacobian is to be loaded into p
c          by columns- i.e., dr(i)/dy(j) is loaded into p(i,j).
c             in the band matrix case ( miter = 4 ), the ele-
c          ments within the band are to be loaded into p by
c          by columns, with diagonal lines of dr/dy loaded into
c          the rows of p.  thus dr(i)/dy(j) is to be loaded
c          into p(i-j+mu+1,j). the locations in p in the two
c          triangular areas which correspond to nonexistent matrix
c          elements can be ignored or loaded arbitrarily, as they
c          they are overwritten by lsodi. ml and mu are the half-
c          bandwidth parameters ( see iwork ).
c               in either case, p is preset to zero by the solver,
c          so that only the nonzero elements need be loaded by jac.
c          each call to jac is preceded by a call to res with the same
c          arguments neq, t, y, and s.  thus to gain some efficiency,
c          intermediate quantities shared by both calculations may be
c          saved in a user common block by res and not recomputed by jac
c          if desired.  also, jac may alter the y array, if desired.
c               jac need not provide dr/dy exactly.  a crude
c          approximation (possibly with a smaller bandwidth) will do.
c               jac must be declared external in the calling program.
c               see note below for more about jac.
c
c       note on res, adda, and jac--  these
c          subroutines may access user-defined quantities in
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array
c          (dimensioned in the subroutines) and/or y has length
c          exceeding neq(1).  however, these routines should not alter
c          neq(1), y(1),...,y(neq) or any other input variables.
c          see the descriptions of neq and y below.
c
c neq    = the size of the system (number of first order ordinary
c          differential equations or scalar algebraic equations).
c          used only for input.
c          neq may be decreased, but not increased, during the problem.
c          if neq is decreased (with istate = 3 on input), the
c          remaining components of y should be left undisturbed, if
c          these are to be accessed in res, adda, or jac.
c
c          normally, neq is a scalar, and it is generally referred to
c          as a scalar in this user interface description.  however,
c          neq may be an array, with neq(1) set to the system size.
c          (the lsodi package accesses only neq(1).)  in either case,
c          this parameter is passed as the neq argument in all calls
c          to res, adda, and jac.  hence, if it is an array,
c          locations neq(2),... may be used to store other integer data
c          and pass it to  res, adda, or jac.  each such subroutine
c          must include neq in a dimension statement in that case.
c
c y      = a real array for the vector of dependent variables, of
c          length neq or more.  used for both input and output on the
c          first call (istate = 0 or 1), and only for output on other
c          calls.  on the first call, y must contain the vector of
c          initial values.  on output, y contains the computed solution
c          vector, evaluated at t.  if desired, the y array may be used
c          for other purposes between calls to the solver.
c
c          this array is passed as the y argument in all calls to res,
c          adda, and jac.  hence its length may exceed neq,
c          and locations y(neq+1),... may be used to store other real
c          data and pass it to res, adda, or jac.  (the lsodi
c          package accesses only y(1),...,y(neq). )
c
c ydoti  = a real array for the initial value of the vector
c          dy/dt and for work space, of dimension at least neq.
c
c          on input...
c            if istate = 0 then lsodi will compute the initial value
c          of dy/dt, if a is nonsingular.  thus ydoti will
c          serve only as work space and may have any value.
c            if istate = 1 then ydoti must contain the initial value
c          of dy/dt.
c            if istate = 2 or 3 (continuation calls) then ydoti
c          may have any value.
c            n.b.- if the initial value of a is singular, then
c          lsodi cannot compute the initial value of dy/dt, so
c          it must be provided in ydoti, with istate=1.
c
c          on output, when lsodi terminates abnormally with istate =
c          -1, -4, or -5, ydoti will contain the residual
c          r = g(t,y) - a(t,y)*(dy/dt).  if r is large, t is near
c          its initial value, and ydoti is supplied with istate=1,
c          there may have been an incorrect input value of
c          ydoti = dy/dt or the problem ( as given to lsodi )
c          may not have a solution.
c
c          if desired, the ydoti array may be used for other
c          purposes between calls to the solver.
c
c t      = the independent variable.  on input, t is used only on the
c          first call, as the initial point of the integration.
c          on output, after each call, t is the value at which a
c          computed solution y is evaluated (usually the same as tout).
c          on an error return, t is the farthest point reached.
c
c tout   = the next value of t at which a computed solution is desired.
c          used only for input.
c
c          when starting the problem (istate = 0 or 1), tout may be
c          equal to t for one call, then should .ne. t for the next
c          call.  for the initial t, an input value of tout .ne. t is
c          used in order to determine the direction of the integration
c          (i.e. the algebraic sign of the step sizes) and the rough
c          scale of the problem.  integration in either direction
c          (forward or backward in t) is permitted.
c
c          if itask = 2 or 5 (one-step modes), tout is ignored after
c          the first call (i.e. the first call with tout .ne. t).
c          otherwise, tout is required on every call.
c
c          if itask = 1, 3, or 4, the values of tout need not be
c          monotone, but a value of tout which backs up is limited
c          to the current internal t interval, whose endpoints are
c          tcur - hu and tcur (see optional outputs, below, for
c          tcur and hu).
c
c itol   = an indicator for the type of error control.  see
c          description below under atol.  used only for input.
c
c rtol   = a relative error tolerance parameter, either a scalar or
c          an array of length neq.  see description below under atol.
c          input only.
c
c atol   = an absolute error tolerance parameter, either a scalar or
c          an array of length neq.  input only.
c
c             the input parameters itol, rtol, and atol determine
c          the error control performed by the solver.  the solver will
c          control the vector e = (e(i)) of estimated local errors
c          in y, according to an inequality of the form
c                      rms-norm of ( e(i)/ewt(i) )   .le.   1,
c          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),
c          and the rms-norm (root-mean-square norm) here is
c          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))
c          is a vector of weights which must always be positive, and
c          the values of rtol and atol should all be non-negative.
c          the following table gives the types (scalar/array) of
c          rtol and atol, and the corresponding form of ewt(i).
c
c             itol    rtol       atol          ewt(i)
c              1     scalar     scalar     rtol*abs(y(i)) + atol
c              2     scalar     array      rtol*abs(y(i)) + atol(i)
c              3     array      scalar     rtol(i)*abs(y(i)) + atol
c              4     array      scalar     rtol(i)*abs(y(i)) + atol(i)
c
c          when either of these parameters is a scalar, it need not
c          be dimensioned in the user-s calling program.
c
c          if none of the above choices (with itol, rtol, and atol
c          fixed throughout the problem) is suitable, more general
c          error controls can be obtained by substituting
c          user-supplied routines for the setting of ewt and/or for
c          the norm calculation.  see part iv below.
c
c          if global errors are to be estimated by making a repeated
c          run on the same problem with smaller tolerances, then all
c          components of rtol and atol (i.e. of ewt) should be scaled
c          down uniformly
c
c itask  = an index specifying the task to be performed.
c          input only.  itask has the following values and meanings.
c          1  means normal computation of output values of y(t) at
c             t = tout (by overshooting and interpolating).
c          2  means take one step only and return.
c          3  means stop at the first internal mesh point at or
c             beyond t = tout and return.
c          4  means normal computation of output values of y(t) at
c             t = tout but without overshooting t = tcrit.
c             tcrit must be input as rwork(1).  tcrit may be equal to
c             or beyond tout, but not behind it in the direction of
c             integration.  this option is useful if the problem
c             has a singularity at or beyond t = tcrit.
c          5  means take one step, without passing tcrit, and return.
c             tcrit must be input as rwork(1).
c
c          note..  if itask = 4 or 5 and the solver reaches tcrit
c          (within roundoff), it will return t = tcrit (exactly) to
c          indicate this (unless itask = 4 and tout comes before tcrit,
c          in which case answers at t = tout are returned first).
c
c istate = an index used for input and output to specify the
c          state of the calculation.
c
c          on input, the values of istate are as follows.
c          0  means this is the first call for the problem, and
c             lsodi is to compute the initial value of dy/dt
c             (while doing other initializations).  see note below.
c          1  means this is the first call for the problem, and
c             the initial value of dy/dt has been supplied in
c             ydoti (lsodi will do other initializations). see note
c             below.
c          2  means this is not the first call, and the calculation
c             is to continue normally, with no change in any input
c             parameters except possibly tout and itask.
c             (if itol, rtol, and/or atol are changed between calls
c             with istate = 2, the new values will be used but not
c             tested for legality.)
c          3  means this is not the first call, and the
c             calculation is to continue normally, but with
c             a change in input parameters other than
c             tout and itask.  changes are allowed in
c             neq, itol, rtol, atol, iopt, lrw, liw, mf, ml, mu,
c             and any of the optional inputs except h0.
c             (see iwork description for ml and mu.)
c          note..  a preliminary call with tout = t is not counted
c          as a first call here, as no initialization or checking of
c          input is done.  (such a call is sometimes useful for the
c          purpose of outputting the initial conditions.)
c          thus the first call for which tout .ne. t requires
c          istate = 0 or 1 on input.
c
c          on output, istate has the following values and meanings.
c           0 or 1  means nothing was done, as tout was equal to t with
c              istate = 0 or 1 on input.  (however, an internal counter
c              was set to detect and prevent repeated calls of this
c              type. )
c           2  means that the integration was performed successfully.
c           3  means that the user-supplied subroutine res signalled
c              lsodi to halt the integration and return (ires=2).
c              integration as far as t was achieved with no occurrence
c              of ires=2, but this flag was set on attempting the next
c              step.
c          -1  means an excessive amount of work (more than mxstep
c              steps) was done on this call, before completing the
c              requested task, but the integration was otherwise
c              successful as far as t.  (mxstep is an optional input
c              and is normally 500.)  to continue, the user may
c              simply reset istate to a value .gt. 1 and call again
c              (the excess work step counter will be reset to 0).
c              in addition, the user may increase mxstep to avoid
c              this error return (see below on optional inputs).
c          -2  means too much accuracy was requested for the precision
c              of the machine being used.  this was detected before
c              completing the requested task, but the integration
c              was successful as far as t.  to continue, the tolerance
c              parameters must be reset, and istate must be set
c              to 3.  the optional output tolsf may be used for this
c              purpose.  (note.. if this condition is detected before
c              taking any steps, then an illegal input return
c              (istate = -3) occurs instead.)
c          -3  means illegal input was detected, before taking any
c              integration steps.  see written message for details.
c              note..  if the solver detects an infinite loop of calls
c              to the solver with illegal input, it will cause
c              the run to stop.
c          -4  means there were repeated error test failures on
c              one attempted step, before completing the requested
c              task, but the integration was successful as far as t.
c              the problem may have a singularity, or the input
c              may be inappropriate.
c          -5  means there were repeated convergence test failures on
c              one attempted step, before completing the requested
c              task, but the integration was successful as far as t.
c              this may be caused by an inaccurate jacobian matrix.
c          -6  means ewt(i) became zero for some i during the
c              integration.  pure relative error control (atol(i)=0.0)
c              was requested on a variable which has now vanished.
c              the integration was successful as far as t.
c          -7  means that the user-supplied subroutine res set
c              its error flag (ires = 3) despite repeated tries by
c              lsodi to avoid that condition.
c          -8  means that istate was 0 on input but lsodi was unable
c              to compute the initial value of dy/dt.  see the
c              printed message for details.
c
c          note..  since the normal output value of istate is 2,
c          it does not need to be reset for normal continuation.
c          similarly, istate need not be reset if res told lsodi
c          to return because the calling program must change
c          the parameters of the problem.
c          also, since a negative input value of istate will be
c          regarded as illegal, a negative output value requires the
c          user to change it, and possibly other inputs, before
c          calling the solver again.
c
c iopt   = an integer flag to specify whether or not any optional
c          inputs are being used on this call.  input only.
c          the optional inputs are listed separately below.
c          iopt = 0 means no optional inputs are being used.
c                   default values will be used in all cases.
c          iopt = 1 means one or more optional inputs are being used.
c
c rwork  = a real working array (double precision).
c          the length of rwork must be at least
c             20 + nyh*(maxord + 1) + 3*neq + lenwm    where
c          nyh    = the initial value of neq,
c          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
c                   smaller value is given as an optional input),
c          lenwm   = neq**2 + 2    if miter is 1 or 2, and
c          lenwm   = (2*ml+mu+1)*neq + 2 if miter is 4 or 5.
c          (see mf description for the definition of meth and miter.)
c          thus if maxord has its default value and neq is constant,
c          this length is
c             22 + 16*neq + neq**2         for mf = 11 or 12,
c             22 + 17*neq + (2*ml+mu)*neq  for mf = 14 or 15,
c             22 +  9*neq + neq**2         for mf = 21 or 22,
c             22 + 10*neq + (2*ml+mu)*neq  for mf = 24 or 25.
c          the first 20 words of rwork are reserved for conditional
c          and optional inputs and optional outputs.
c
c          the following word in rwork is a conditional input..
c            rwork(1) = tcrit = critical value of t which the solver
c                       is not to overshoot.  required if itask is
c                       4 or 5, and ignored otherwise.  (see itask.)
c
c lrw    = the length of the array rwork, as declared by the user.
c          (this will be checked by the solver.)
c
c iwork  = an integer work array.  the length of iwork must be at least
c          20 + neq .  the first few words of iwork are used for
c          conditional and optional inputs and optional outputs.
c
c          the following 2 words in iwork are conditional inputs..
c            iwork(1) = ml     these are the lower and upper
c            iwork(2) = mu     half-bandwidths, respectively, of the
c                       matrices in the problem-- the jacobian dr/dy
c                       and the left-hand side matrix a. these half-
c                       bandwidths exclude the main diagonal, so
c                       the total bandwidth is ml + mu + 1 .
c                       the band is defined by the matrix locations
c                       (i,j) with i-ml .le. j .le. i+mu.  ml and mu
c                       must satisfy  0 .le.  ml,mu  .le. neq-1.
c                       these are required if miter is 4 or 5, and
c                       ignored otherwise.
c                       ml and mu may in fact be the band parameters for
c                       matrices to which dr/dy and a are only
c                       approximately equal.
c
c liw    = the length of the array iwork, as declared by the user.
c          (this will be checked by the solver.)
c
c note..  the work arrays must not be altered between calls to lsodi
c for the same problem, except possibly for the conditional and
c optional inputs, and except for the last 3*neq words of rwork.
c the latter space is used for internal scratch space, and so is
c available for use by the user outside lsodi between calls, if
c desired (but not for use by res, adda, or jac).
c
c mf     = the method flag.  used only for input.  the legal values of
c          mf are 11, 12, 14, 15, 21, 22, 24, and 25.
c          mf has decimal digits meth and miter.. mf = 10*meth + miter.
c            meth indicates the basic linear multistep method..
c              meth = 1 means the implicit adams method.
c              meth = 2 means the method based on backward
c                       differentiation formulas (bdf-s).
c                the bdf method is strongly preferred for stiff prob-
c              lems, while the adams method is preferred when the prob-
c              lem is not stiff. if the matrix a(t,y) is nonsingular,
c              stiffness here can be taken to mean that of the explicit
c              ode system dy/dt = a**(-1) * g.  if a is singular, the
c              concept of stiffness is not well defined.
c                if you do not know whether the problem is stiff, we
c              recommend using meth = 2.  if it is stiff, the advan-
c              tage of meth = 2 over 1 will be great, while if it is
c              not stiff, the advantage of meth = 1 will be slight.
c              if maximum efficiency is important, some experimentation
c              with meth may be necessary.
c            miter indicates the corrector iteration method..
c              miter = 1 means chord iteration with a user-supplied
c                        full (neq by neq) jacobian.
c              miter = 2 means chord iteration with an internally
c                        generated (difference quotient) full jacobian.
c                        this uses neq+1 extra calls to res per dr/dy
c                        evaluation.
c              miter = 4 means chord iteration with a user-supplied
c                        banded jacobian.
c              miter = 5 means chord iteration with an internally
c                        generated banded jacobian (using ml+mu+2
c                        extra calls to res per dr/dy evaluation).
c              if miter = 1 or 4, the user must supply a subroutine jac
c              (the name is arbitrary) as described above under jac.
c              for other values of miter, a dummy argument can be used.
c-----------------------------------------------------------------------
c optional inputs.
c
c the following is a list of the optional inputs provided for in the
c call sequence.  (see also part ii.)  for each such input variable,
c this table lists its name as used in this documentation, its
c location in the call sequence, its meaning, and the default value.
c the use of any of these inputs requires iopt = 1, and in that
c case all of these inputs are examined.  a value of zero for any
c of these optional inputs will cause the default value to be used.
c thus to use a subset of the optional inputs, simply preload
c locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
c then set those of interest to nonzero values.
c
c name    location      meaning and default value
c
c h0      rwork(5)  the step size to be attempted on the first step.
c                   the default value is determined by the solver.
c
c hmax    rwork(6)  the maximum absolute step size allowed.
c                   the default value is infinite.
c
c hmin    rwork(7)  the minimum absolute step size allowed.
c                   the default value is 0.  (this lower bound is not
c                   enforced on the final step before reaching tcrit
c                   when itask = 4 or 5.)
c
c maxord  iwork(5)  the maximum order to be allowed.  the default
c                   value is 12 if meth = 1, and 5 if meth = 2.
c                   if maxord exceeds the default value, it will
c                   be reduced to the default value.
c                   if maxord is changed during the problem, it may
c                   cause the current order to be reduced.
c
c mxstep  iwork(6)  maximum number of (internally defined) steps
c                   allowed during one call to the solver.
c                   the default value is 500.
c
c mxhnil  iwork(7)  maximum number of messages printed (per problem)
c                   warning that t + h = t on a step (h = step size).
c                   this must be positive to result in a non-default
c                   value.  the default value is 10.
c-----------------------------------------------------------------------
c optional outputs.
c
c as optional additional output from lsodi, the variables listed
c below are quantities related to the performance of lsodi
c which are available to the user.  these are communicated by way of
c the work arrays, but also have internal mnemonic names as shown.
c except where stated otherwise, all of these outputs are defined
c on any successful return from lsodi, and on any return with
c istate = -1, -2, -4, -5, -6, or -7. on a return with -3 (illegal
c input) or -8, they will be unchanged from their existing values
c (if any), except possibly for tolsf, lenrw, and leniw.
c on any error return, outputs relevant to the error will be defined,
c as noted below.
c
c name    location      meaning
c
c hu      rwork(11) the step size in t last used (successfully).
c
c hcur    rwork(12) the step size to be attempted on the next step.
c
c tcur    rwork(13) the current value of the independent variable
c                   which the solver has actually reached, i.e. the
c                   current internal mesh point in t.  on output, tcur
c                   will always be at least as far as the argument
c                   t, but may be farther (if interpolation was done).
c
c tolsf   rwork(14) a tolerance scale factor, greater than 1.0,
c                   computed when a request for too much accuracy was
c                   detected (istate = -3 if detected at the start of
c                   the problem, istate = -2 otherwise).  if itol is
c                   left unaltered but rtol and atol are uniformly
c                   scaled up by a factor of tolsf for the next call,
c                   then the solver is deemed likely to succeed.
c                   (the user may also ignore tolsf and alter the
c                   tolerance parameters in any other way appropriate.)
c
c nst     iwork(11) the number of steps taken for the problem so far.
c
c nre     iwork(12) the number of residual evaluations (res calls)
c                   for the problem so far.
c
c nje     iwork(13) the number of jacobian evaluations (each involving
c                   an evaluation of a and dr/dy) for the problem so
c                   far.  this equals the number of calls to adda and
c                   (if miter = 1 or 4) jac, and the number of matrix
c                   l-u decompositions.
c
c nqu     iwork(14) the method order last used (successfully).
c
c nqcur   iwork(15) the order to be attempted on the next step.
c
c imxer   iwork(16) the index of the component of largest magnitude in
c                   the weighted local error vector ( e(i)/ewt(i) ),
c                   on an error return with istate = -4 or -5.
c
c lenrw   iwork(17) the length of rwork actually required.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c leniw   iwork(18) the length of iwork actually required.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c
c the following two arrays are segments of the rwork array which
c may also be of interest to the user as optional outputs.
c for each array, the table below gives its internal name,
c its base address in rwork, and its description.
c
c name    base address      description
c
c yh      21             the nordsieck history array, of size nyh by
c                        (nqcur + 1), where nyh is the initial value
c                        of neq.  for j = 0,1,...,nqcur, column j+1
c                        of yh contains hcur**j/factorial(j) times
c                        the j-th derivative of the interpolating
c                        polynomial currently representing the solution,
c                        evaluated at t = tcur.
c
c acor     lenrw-neq+1   array of size neq used for the accumulated
c                        corrections on each step, scaled on output to
c                        represent the estimated local error in y on the
c                        last step. this is the vector e in the descrip-
c                        tion of the error control.  it is defined only
c                        on a return from lsodi with istate = 2.
c
c-----------------------------------------------------------------------
c part ii.  other routines callable.
c
c the following are optional calls which the user may make to
c gain additional capabilities in conjunction with lsodi.
c (the routines xsetun and xsetf are designed to conform to the
c slatec error handling package.)
c
c     form of call                  function
c   call xsetun(lun)          set the logical unit number, lun, for
c                             output of messages from lsodi, if
c                             the default is not desired.
c                             the default value of lun is 6.
c
c   call xsetf(mflag)         set a flag to control the printing of
c                             messages by lsodi.
c                             mflag = 0 means do not print. (danger..
c                             this risks losing valuable information.)
c                             mflag = 1 means print (the default).
c
c                             either of the above calls may be made at
c                             any time and will take effect immediately.
c
c   call srcom(rsav,isav,job) saves and restores the contents of
c                             the internal common blocks used by
c                             lsodi (see part iii below).
c                             rsav must be a real array of length 218
c                             or more, and isav must be an integer
c                             array of length 41 or more.
c                             job=1 means save common into rsav/isav.
c                             job=2 means restore common from rsav/isav.
c                                srcom is useful if one is
c                             interrupting a run and restarting
c                             later, or alternating between two or
c                             more problems solved with lsodi.
c
c   call intdy(,,,,,)         provide derivatives of y, of various
c        (see below)          orders, at a specified point t, if
c                             desired.  it may be called only after
c                             a successful return from lsodi.
c
c the detailed instructions for using intdy are as follows.
c the form of the call is..
c
c   call intdy (t, k, rwork(21), nyh, dky, iflag)
c
c the input parameters are..
c
c t         = value of independent variable where answers are desired
c             (normally the same as the t last returned by lsodi).
c             for valid results, t must lie between tcur - hu and tcur.
c             (see optional outputs for tcur and hu.)
c k         = integer order of the derivative desired.  k must satisfy
c             0 .le. k .le. nqcur, where nqcur is the current order
c             (see optional outputs).  the capability corresponding
c             to k = 0, i.e. computing y(t), is already provided
c             by lsodi directly.  since nqcur .ge. 1, the first
c             derivative dy/dt is always available with intdy.
c rwork(21) = the base address of the history array yh.
c nyh       = column length of yh, equal to the initial value of neq.
c
c the output parameters are..
c
c dky       = a real array of length neq containing the computed value
c             of the k-th derivative of y(t).
c iflag     = integer flag, returned as 0 if k and t were legal,
c             -1 if k was illegal, and -2 if t was illegal.
c             on an error return, a message is also written.
c-----------------------------------------------------------------------
c part iii.  common blocks.
c
c if lsodi is to be used in an overlay situation, the user
c must declare, in the primary overlay, the variables in..
c   (1) the call sequence to lsodi,
c   (2) the two internal common blocks
c         /ls0001/  of length  257  (218 double precision words
c                         followed by 39 integer words),
c         /eh0001/  of length  2 (integer words).
c
c if lsodi is used on a system in which the contents of internal
c common blocks are not preserved between calls, the user should
c declare the above two common blocks in his main program to insure
c that their contents are preserved.
c
c if the solution of a given problem by lsodi is to be interrupted
c and then later continued, such as when restarting an interrupted run
c or alternating between two or more problems, the user should save,
c following the return from the last lsodi call prior to the
c interruption, the contents of the call sequence variables and the
c internal common blocks, and later restore these values before the
c next lsodi call for that problem.  to save and restore the common
c blocks, use subroutine srcom (see part ii above).
c
c-----------------------------------------------------------------------
c part iv.  optionally replaceable solver routines.
c
c below are descriptions of two routines in the lsodi package which
c relate to the measurement of errors.  either routine can be
c replaced by a user-supplied version, if desired.  however, since such
c a replacement may have a major impact on performance, it should be
c done only when absolutely necessary, and only with great caution.
c (note.. the means by which the package version of a routine is
c superseded by the user-s version may be system-dependent.)
c
c (a) ewset.
c the following subroutine is called just before each internal
c integration step, and sets the array of error weights, ewt, as
c described under itol/rtol/atol above..
c     subroutine ewset (neq, itol, rtol, atol, ycur, ewt)
c where neq, itol, rtol, and atol are as in the lsodi call sequence,
c ycur contains the current dependent variable vector, and
c ewt is the array of weights set by ewset.
c
c if the user supplies this subroutine, it must return in ewt(i)
c (i = 1,...,neq) a positive quantity suitable for comparing errors
c in y(i) to.  the ewt array returned by ewset is passed to the
c vnorm routine (see below), and also used by lsodi in the computation
c of the optional output imxer, the diagonal jacobian approximation,
c and the increments for difference quotient jacobians.
c
c in the user-supplied version of ewset, it may be desirable to use
c the current values of derivatives of y.  derivatives up to order nq
c are available from the history array yh, described above under
c optional outputs.  in ewset, yh is identical to the ycur array,
c extended to nq + 1 columns with a column length of nyh and scale
c factors of h**j/factorial(j).  on the first call for the problem,
c given by nst = 0, nq is 1 and h is temporarily set to 1.0.
c the quantities nq, nyh, h, and nst can be obtained by including
c in ewset the statements..
c     double precision h, rls
c     common /ls0001/ rls(218),ils(39)
c     nq = ils(35)
c     nyh = ils(14)
c     nst = ils(36)
c     h = rls(212)
c thus, for example, the current value of dy/dt can be obtained as
c ycur(nyh+i)/h  (i=1,...,neq)  (and the division by h is
c unnecessary when nst = 0).
c
c (b) vnorm.
c the following is a real function routine which computes the weighted
c root-mean-square norm of a vector v..
c     d = vnorm (n, v, w)
c where..
c   n = the length of the vector,
c   v = real array of length n containing the vector,
c   w = real array of length n containing weights,
c   d = sqrt( (1/n) * sum(v(i)*w(i))**2 ).
c vnorm is called with n = neq and with w(i) = 1.0/ewt(i), where
c ewt is as set by subroutine ewset.
c
c if the user supplies this function, it should return a non-negative
c value of vnorm suitable for use in the error control in lsodi.
c none of the arguments should be altered by vnorm.
c for example, a user-supplied vnorm routine might..
c   -substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
c   -ignore some components of v in the norm, with the effect of
c    suppressing the error control on those components of y.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c other routines in the lsodi package.
c
c in addition to subroutine lsodi, the lsodi package includes the
c following subroutines and function routines..
c  ainvg    computes the initial value of the vector
c             dy/dt = inverse(a) * g
c  intdy    computes an interpolated value of the y vector at t = tout.
c  stodi    is the core integrator, which does one step of the
c           integration and the associated error control.
c  cfode    sets all method coefficients and test constants.
c  prepji   computes and preprocesses the jacobian matrix
c           and the newton iteration matrix p.
c  solsy    manages solution of linear system in chord iteration.
c  ewset    sets the error weight vector ewt before each step.
c  vnorm    computes the weighted r.m.s. norm of a vector.
c  srcom    is a user-callable routine to save and restore
c           the contents of the internal common blocks.
c  dgefa and dgesl   are routines from linpack for solving full
c           systems of linear algebraic equations.
c  dgbfa and dgbsl   are routines from linpack for solving banded
c           linear systems.
c  daxpy, dscal, idamax, and ddot   are basic linear algebra modules
c           (blas) used by the above linpack routines.
c  d1mach   computes the unit roundoff in a machine-independent manner.
c  xerrwv, xsetun, and xsetf   handle the printing of all error
c           messages and warnings.  xerrwv is machine-dependent.
c note..  vnorm, idamax, ddot, and d1mach are function routines.
c all the others are subroutines.
c
c the intrinsic and external routines used by lsodi are..  dabs,
c dmax1, dmin1, dfloat, iabs, max0, min0, mod, dsign, dsqrt, and write.
c
c a block data subprogram is also included with the package,
c for loading some of the variables in internal common.
c
c-----------------------------------------------------------------------
c the following card is for optimized compilation on llnl compilers.
clll. optimize
c-----------------------------------------------------------------------
      external prepji, solsy
      integer illin, init, lyh, lewt, lacor, lsavr, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nre, nje, nqu
      integer i, i1, i2, ier, iflag, imxer, ires, kgo,
     1   leniw, lenrw, lenwm, lp, lyd0, ml, mord, mu, mxhnl0, mxstp0
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0,
     2   d1mach, vnorm
      dimension mord(2)
      logical ihit
c-----------------------------------------------------------------------
c the following internal common block contains
c (a) variables which are local to any subroutine but whose values must
c     be preserved between calls to the routine (own variables), and
c (b) variables which are communicated between subroutines.
c common block ls0001 is shared by the lsodi and lsode packages.
c the structure of ls0001 is as follows..  all real variables are
c listed first, followed by all integers.  within each type, the
c variables are grouped with those local to subroutine lsodi first,
c then those local to subroutine stodi, and finally those used
c for communication.  the block is declared in subroutines
c lsodi, intdy, stodi, prepji, and solsy.  groups of variables are
c replaced by dummy arrays in the common declarations in routines
c where those variables are not used.
c-----------------------------------------------------------------------
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavr, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nre, nje, nqu
c
      data  mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
c-----------------------------------------------------------------------
c block a.
c this code block is executed on every call.
c it tests istate and itask for legality and branches appropriately.
c if istate .gt. 1 but the flag init shows that initialization has
c not yet been done, an error return occurs.
c if istate = 0 or 1 and tout = t, jump to block g and return
c immediately.
c-----------------------------------------------------------------------
      if (istate .lt. 0 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .le. 1) go to 10
      if (init .eq. 0) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   init = 0
      if (tout .eq. t) go to 430
 20   ntrep = 0
c-----------------------------------------------------------------------
c block b.
c the next code block is executed for the initial call (istate = 0 or 1)
c or for a continuation call with parameter changes (istate = 3).
c it contains checking of all inputs and various initializations.
c
c first check legality of the non-optional inputs neq, itol, iopt,
c mf, ml, and mu.
c-----------------------------------------------------------------------
      if (neq(1) .le. 0) go to 604
      if (istate .le. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      meth = mf/10
      miter = mf - 10*meth
      if (meth .lt. 1 .or. meth .gt. 2) go to 608
      if (miter .le. 0 .or. miter .gt. 5) go to 608
      if (miter .eq. 3)  go to 608
      if (miter .lt. 3) go to 30
      ml = iwork(1)
      mu = iwork(2)
      if (ml .lt. 0 .or. ml .ge. n) go to 609
      if (mu .lt. 0 .or. mu .ge. n) go to 610
 30   continue
c next process and check the optional inputs. --------------------------
      if (iopt .eq. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      if (istate .le. 1) h0 = 0.0d0
      hmxi = 0.0d0
      hmin = 0.0d0
      go to 60
 40   maxord = iwork(5)
      if (maxord .lt. 0) go to 611
      if (maxord .eq. 0) maxord = 100
      maxord = min0(maxord,mord(meth))
      mxstep = iwork(6)
      if (mxstep .lt. 0) go to 612
      if (mxstep .eq. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      if (mxhnil .lt. 0) go to 613
      if (mxhnil .eq. 0) mxhnil = mxhnl0
      if (istate .gt. 1) go to 50
      h0 = rwork(5)
      if ((tout - t)*h0 .lt. 0.0d0) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. 0.0d0) go to 615
      hmxi = 0.0d0
      if (hmax .gt. 0.0d0) hmxi = 1.0d0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0d0) go to 616
c-----------------------------------------------------------------------
c set work array pointers and check lengths lrw and liw.
c pointers to segments of rwork and iwork are named by prefixing l to
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).
c segments of rwork (in order) are denoted yh, wm, ewt, savr, acor.
c-----------------------------------------------------------------------
 60   lyh = 21
      if (istate .le. 1) nyh = n
      lwm = lyh + (maxord + 1)*nyh
      if (miter .le. 2) lenwm = n*n + 2
      if (miter .ge. 4) lenwm = (2*ml + mu + 1)*n + 2
      lewt = lwm + lenwm
      lsavr = lewt + n
      lacor = lsavr + n
      lenrw = lacor + n - 1
      iwork(17) = lenrw
      liwm = 1
      leniw = 20 + n
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618
c check rtol and atol for legality. ------------------------------------
      rtoli = rtol(1)
      atoli = atol(1)
      do 70 i = 1,n
        if (itol .ge. 3) rtoli = rtol(i)
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        if (rtoli .lt. 0.0d0) go to 619
        if (atoli .lt. 0.0d0) go to 620
 70     continue
      if (istate .le. 1) go to 100
c if istate = 3, set flag to signal parameter changes to stodi. --------
      jstart = -1
      if (nq .le. maxord) go to 90
c maxord was reduced below nq.  copy yh(*,maxord+2) into ydoti.---------
      do 80 i = 1,n
 80     ydoti(i) = rwork(i+lwm-1)
c reload wm(1) = rwork(lwm), since lwm may have changed. ---------------
 90   rwork(lwm) = dsqrt(uround)
      if (n .eq. nyh) go to 200
c neq was reduced.  zero part of yh to avoid undefined references. -----
      i1 = lyh + l*nyh
      i2 = lyh + (maxord + 1)*nyh - 1
      if (i1 .gt. i2) go to 200
      do 95 i = i1,i2
 95     rwork(i) = 0.0d0
      go to 200
c-----------------------------------------------------------------------
c block c.
c the next block is for the initial call only (istate = 0 or 1).
c it contains all remaining initializations, the call to ainvg
c (if istate = 1), and the calculation of the initial step size.
c the error weights in ewt are inverted after being loaded.
c-----------------------------------------------------------------------
 100  uround = d1mach(4)
      tn = t
      if (itask .ne. 4 .and. itask .ne. 5) go to 105
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. 0.0d0) go to 625
      if (h0 .ne. 0.0d0 .and. (t + h0 - tcrit)*h0 .gt. 0.0d0)
     1   h0 = tcrit - t
 105  jstart = 0
      rwork(lwm) = dsqrt(uround)
      nhnil = 0
      nst = 0
      nre = 0
      nje = 0
      nslast = 0
      hu = 0.0d0
      nqu = 0
      ccmax = 0.3d0
      maxcor = 3
      msbp = 20
      mxncf = 10
c compute initial dy/dt, if necessary, and load it and initial y into yh
      lyd0 = lyh + nyh
      lp = lwm + 1
      if (istate .eq. 1) go to 120
c lsodi must compute initial dy/dt (lyd0 points to yh(*,2)). -----------
         call ainvg( res, adda, neq, t, y, rwork(lyd0), miter,
     1               ml, mu, rwork(lp), iwork(21), ier )
         nre = nre + 1
         if (ier.lt.0) go to 560
         if (ier.eq.0) go to 110
         go to 565
 110     continue
         do 115 i = 1,n
 115        rwork(i+lyh-1) = y(i)
         go to 130
c initial dy/dt has been supplied. -------------------------------------
 120     do 125 i = 1,n
            rwork(i+lyh-1) = y(i)
 125        rwork(i+lyd0-1) = ydoti(i)
c load and invert the ewt array.  (h is temporarily set to 1.0.) -------
 130  continue
      nq = 1
      h = 1.0d0
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 135 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0d0) go to 621
 135    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
c-----------------------------------------------------------------------
c the coding below computes the step size, h0, to be attempted on the
c first step, unless the user has supplied a value for this.
c first check that tout - t differs significantly from zero.
c a scalar tolerance quantity tol is computed, as max(rtol(i))
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
c so as to be between 100*uround and 1.0e-3.
c then the computed value h0 is given by..
c                                      neq
c   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( ydot(i)/ywt(i) )**2  )
c                                       1
c where   w0      = max ( abs(t), abs(tout) ),
c         ydot(i) = i-th component of initial value of dy/dt,
c         ywt(i)  = ewt(i)/tol  (a weight for y(i)).
c the sign of h0 is inferred from the initial values of tout and t.
c-----------------------------------------------------------------------
      if (h0 .ne. 0.0d0) go to 180
      tdist = dabs(tout - t)
      w0 = dmax1(dabs(t),dabs(tout))
      if (tdist .lt. 2.0d0*uround*w0) go to 622
      tol = rtol(1)
      if (itol .le. 2) go to 145
      do 140 i = 1,n
 140    tol = dmax1(tol,rtol(i))
 145  if (tol .gt. 0.0d0) go to 160
      atoli = atol(1)
      do 150 i = 1,n
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        ayi = dabs(y(i))
        if (ayi .ne. 0.0d0) tol = dmax1(tol,atoli/ayi)
 150    continue
 160  tol = dmax1(tol,100.0d0*uround)
      tol = dmin1(tol,0.001d0)
      sum = vnorm (n, rwork(lyd0), rwork(lewt))
      sum = 1.0d0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0d0/dsqrt(sum)
      h0 = dmin1(h0,tdist)
      h0 = dsign(h0,tout-t)
c adjust h0 if necessary to meet hmax bound. ---------------------------
 180  rh = dabs(h0)*hmxi
      if (rh .gt. 1.0d0) h0 = h0/rh
c load h with h0 and scale yh(*,2) by h0. ------------------------------
      h = h0
      do 190 i = 1,n
 190    rwork(i+lyd0-1) = h0*rwork(i+lyd0-1)
      go to 270
c-----------------------------------------------------------------------
c block d.
c the next code block is for continuation calls only (istate = 2 or 3)
c and is to check stop conditions before taking a step.
c-----------------------------------------------------------------------
 200  nslast = nst
      go to (210, 250, 220, 230, 240), itask
 210  if ((tn - tout)*h .lt. 0.0d0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 220  tp = tn - hu*(1.0d0 + 100.0d0*uround)
      if ((tp - tout)*h .gt. 0.0d0) go to 623
      if ((tn - tout)*h .lt. 0.0d0) go to 250
      go to 400
 230  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0d0) go to 624
      if ((tcrit - tout)*h .lt. 0.0d0) go to 625
      if ((tn - tout)*h .lt. 0.0d0) go to 245
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 240  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0d0) go to 624
 245  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0d0 + 4.0d0*uround)
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
      if (istate .eq. 2) jstart = -2
c-----------------------------------------------------------------------
c block e.
c the next block is normally executed for all calls and contains
c the call to the one-step core integrator stodi.
c
c this is a looping point for the integration steps.
c
c first check for too many steps being taken, update ewt (if not at
c start of problem), check for too much accuracy being requested, and
c check for h below the roundoff level in t.
c-----------------------------------------------------------------------
 250  continue
      if ((nst-nslast) .ge. mxstep) go to 500
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 260 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0d0) go to 510
 260    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
 270  tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))
      if (tolsf .le. 1.0d0) go to 280
      tolsf = tolsf*2.0d0
      if (nst .eq. 0) go to 626
      go to 520
 280  if ((tn + h) .ne. tn) go to 290
      nhnil = nhnil + 1
      if (nhnil .gt. mxhnil) go to 290
      call xerrwv('lsodi--  warning..internal t (=r1) and h (=r2) are',
     1   50, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      such that in the machine, t + h = t on the next step  ',
     1   60, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      (h = step size). solver will continue anyway',
     1   50, 101, 0, 0, 0, 0, 2, tn, h)
      if (nhnil .lt. mxhnil) go to 290
      call xerrwv('lsodi--  above warning has been issued i1 times.  ',
     1   50, 102, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      it will not be issued again for this problem',
     1   50, 102, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
 290  continue
c-----------------------------------------------------------------------
c     call stodi(neq,y,yh,nyh,yh1,ewt,savf,savr,acor,wm,iwm,res,
c                adda,jac,prepji,solsy)
c note... savf in stodi occupies the same space as ydoti in lsodi.
c-----------------------------------------------------------------------
      call stodi (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),
     1   ydoti, rwork(lsavr), rwork(lacor), rwork(lwm),
     2   iwork(liwm), res, adda, jac, prepji, solsy )
      kgo = 1 - kflag
      go to (300, 530, 540, 400, 550), kgo
c
c kgo = 1,success. 2,error test failure. 3,convergence failure.
c       4,res ordered return. 5,res returned error.
c-----------------------------------------------------------------------
c block f.
c the following block handles the case of a successful return from the
c core integrator (kflag = 0).  test for stop conditions.
c-----------------------------------------------------------------------
 300  init = 1
      go to (310, 400, 330, 340, 350), itask
c itask = 1.  if tout has been reached, interpolate. -------------------
 310  if ((tn - tout)*h .lt. 0.0d0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
c itask = 3.  jump to exit if tout was reached. ------------------------
 330  if ((tn - tout)*h .ge. 0.0d0) go to 400
      go to 250
c itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
 340  if ((tn - tout)*h .lt. 0.0d0) go to 345
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
 345  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0d0 + 4.0d0*uround)
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
      jstart = -2
      go to 250
c itask = 5.  see if tcrit was reached and jump to exit. ---------------
 350  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
c-----------------------------------------------------------------------
c block g.
c the following block handles all successful returns from lsodi.
c if itask .ne. 1, y is loaded from yh and t is set accordingly.
c istate is set to 2, the illegal input counter is zeroed, and the
c optional outputs are loaded into the work arrays before returning.  if
c istate = 0 or 1 and tout = t, there is a return with no action taken,
c except that if this has happened repeatedly, the run is terminated.
c-----------------------------------------------------------------------
 400  do 410 i = 1,n
 410    y(i) = rwork(i+lyh-1)
      t = tn
      if (itask .ne. 4 .and. itask .ne. 5) go to 420
      if (ihit) t = tcrit
 420  istate = 2
      if (kflag .eq. -3) istate = 3
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nre
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c
 430  ntrep = ntrep + 1
      if (ntrep .lt. 5) return
      call xerrwv(
     1  'lsodi--  repeated calls with istate= 0 or 1 and tout= t(=r1)',
     1   60, 301, 0, 0, 0, 0, 1, t, 0.0d0)
      go to 800
c-----------------------------------------------------------------------
c block h.
c the following block handles all unsuccessful returns other than
c those for illegal input.  first the error message routine is called.
c if there was an error test or convergence test failure, imxer is set.
c then y is loaded from yh, t is set to tn, and the illegal input
c counter illin is set to 0.  the optional outputs are loaded into
c the work arrays before returning.
c-----------------------------------------------------------------------
c the maximum number of steps was taken before reaching tout. ----------
 500  call xerrwv('lsodi--  at current t (=r1), mxstep (=i1) steps   ',
     1   50, 201, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      taken on this call before reaching tout     ',
     1   50, 201, 0, 1, mxstep, 0, 1, tn, 0.0d0)
      istate = -1
      go to 580
c ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  ewti = rwork(lewt+i-1)
      call xerrwv('lsodi--  at t (=r1), ewt(i1) has become r2 .le. 0.',
     1   50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 590
c too much accuracy requested for machine precision. -------------------
 520  call xerrwv('lsodi--  at t (=r1), too much accuracy requested  ',
     1   50, 203, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      for precision of machine..  see tolsf (=r2) ',
     1   50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 590
c kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
 530  call xerrwv('lsodi--  at t(=r1) and step size h(=r2), the error',
     1   50, 204, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      test failed repeatedly or with abs(h) = hmin',
     1   50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 570
c kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
 540  call xerrwv('lsodi--  at t (=r1) and step size h (=r2), the    ',
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      corrector convergence failed repeatedly     ',
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      or with abs(h) = hmin   ',
     1   30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
      go to 570
c ires = 3 returned by res, despite retries by stodi. ------------------
 550  call xerrwv('lsodi--  at t (=r1) residual routine returned     ',
     1   50, 206, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      error ires = 3 repeatedly         ',
     1   40, 206, 0, 0, 0, 0, 1, tn, 0.0d0)
      istate = -7
      go to 590
c ainvg failed because a-matrix was singular. --------------------------
 560  ier = -ier
      call xerrwv(
     1  'lsodi-- attempt to initialize dy/dt failed.. matrix a is    ',
     1   60, 207, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      singular.  sgefa or sgbfa returned info=(i1)',
     2   50, 207, 0, 1, ier, 0, 0, 0.0d0, 0.0d0)
      istate = -8
      return
c ainvg failed because res set ires to 2 or 3. -------------------------
 565  call xerrwv('lsodi--  attempt to initialize dy/dt failed       ',
     1   50, 208, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      because residual routine set its error flag ',
     1   50, 208, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      to ires = (i1)',
     1   20, 208, 0, 1, ier, 0, 0, 0.0d0, 0.0d0)
      istate = -8
      return
c compute imxer if relevant. -------------------------------------------
 570  big = 0.0d0
      imxer = 1
      do 575 i = 1,n
        size = dabs(rwork(i+lacor-1)*rwork(i+lewt-1))
        if (big .ge. size) go to 575
        big = size
        imxer = i
 575    continue
      iwork(16) = imxer
c compute residual if relevant. ----------------------------------------
 580  lyd0 = lyh + nyh
      do 585 i = 1,n
         rwork(i+lsavr-1) = rwork(i+lyd0-1)/h
 585     y(i) = rwork(i+lyh-1)
      ires = 1
      call res ( neq, tn, y, rwork(lsavr), ydoti, ires )
      nre = nre + 1
      if (ires .le. 1) go to 595
      call xerrwv('lsodi--  residual routine set its flag ires       ',
     1   50, 210, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      to (i1) when called for final output.       ',
     1   50, 210, 0, 1, ires, 0, 0, 0.0d0, 0.0d0)
      go to 595
c set y vector, t, illin, and optional outputs. ------------------------
 590  do 592 i = 1,n
 592    y(i) = rwork(i+lyh-1)
 595  t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nre
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c-----------------------------------------------------------------------
c block i.
c the following block handles all error returns due to illegal input
c (istate = -3), as detected before calling the core integrator.
c first the error message routine is called.  then if there have been
c 5 consecutive such returns just before this call to the solver,
c the run is halted.
c-----------------------------------------------------------------------
 601  call xerrwv('lsodi--  istate (=i1) illegal ',
     1   30, 1, 0, 1, istate, 0, 0, 0.0d0, 0.0d0)
      go to 700
 602  call xerrwv('lsodi--  itask (=i1) illegal  ',
     1   30, 2, 0, 1, itask, 0, 0, 0.0d0, 0.0d0)
      go to 700
 603  call xerrwv('lsodi--  istate .gt. 1 but lsodi not initialized  ',
     1   50, 3, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      go to 700
 604  call xerrwv('lsodi--  neq (=i1) .lt. 1     ',
     1   30, 4, 0, 1, neq(1), 0, 0, 0.0d0, 0.0d0)
      go to 700
 605  call xerrwv('lsodi--  istate = 3 and neq increased (i1 to i2)  ',
     1   50, 5, 0, 2, n, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 606  call xerrwv('lsodi--  itol (=i1) illegal   ',
     1   30, 6, 0, 1, itol, 0, 0, 0.0d0, 0.0d0)
      go to 700
 607  call xerrwv('lsodi--  iopt (=i1) illegal   ',
     1   30, 7, 0, 1, iopt, 0, 0, 0.0d0, 0.0d0)
      go to 700
 608  call xerrwv('lsodi--  mf (=i1) illegal     ',
     1   30, 8, 0, 1, mf, 0, 0, 0.0d0, 0.0d0)
      go to 700
 609  call xerrwv('lsodi--  ml(=i1) illegal.. .lt. 0 or .ge. neq(=i2)',
     1   50, 9, 0, 2, ml, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 610  call xerrwv('lsodi--  mu(=i1) illegal.. .lt. 0 or .ge. neq(=i2)',
     1   50, 10, 0, 2, mu, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 611  call xerrwv('lsodi--  maxord (=i1) .lt. 0  ',
     1   30, 11, 0, 1, maxord, 0, 0, 0.0d0, 0.0d0)
      go to 700
 612  call xerrwv('lsodi--  mxstep (=i1) .lt. 0  ',
     1   30, 12, 0, 1, mxstep, 0, 0, 0.0d0, 0.0d0)
      go to 700
 613  call xerrwv('lsodi--  mxhnil (=i1) .lt. 0  ',
     1   30, 13, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
      go to 700
 614  call xerrwv('lsodi--  tout (=r1) behind t (=r2)      ',
     1   40, 14, 0, 0, 0, 0, 2, tout, t)
      call xerrwv('      integration direction is given by h0 (=r1)  ',
     1   50, 14, 0, 0, 0, 0, 1, h0, 0.0d0)
      go to 700
 615  call xerrwv('lsodi--  hmax (=r1) .lt. 0.0  ',
     1   30, 15, 0, 0, 0, 0, 1, hmax, 0.0d0)
      go to 700
 616  call xerrwv('lsodi--  hmin (=r1) .lt. 0.0  ',
     1   30, 16, 0, 0, 0, 0, 1, hmin, 0.0d0)
      go to 700
 617  call xerrwv(
     1  'lsodi--  rwork length needed, lenrw (=i1), exceeds lrw (=i2)',
     1   60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 618  call xerrwv(
     1  'lsodi--  iwork length needed, leniw (=i1), exceeds liw (=i2)',
     1   60, 18, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
      go to 700
 619  call xerrwv('lsodi--  rtol(=i1) is r1 .lt. 0.0       ',
     1   40, 19, 0, 1, i, 0, 1, rtoli, 0.0d0)
      go to 700
 620  call xerrwv('lsodi--  atol(=i1) is r1 .lt. 0.0       ',
     1   40, 20, 0, 1, i, 0, 1, atoli, 0.0d0)
      go to 700
 621  ewti = rwork(lewt+i-1)
      call xerrwv('lsodi--  ewt(=i1) is r1 .le. 0.0        ',
     1   40, 21, 0, 1, i, 0, 1, ewti, 0.0d0)
      go to 700
 622  call xerrwv(
     1  'lsodi--  tout (=r1) too close to t(=r2) to start integration',
     1   60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  call xerrwv(
     1  'lsodi--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  ',
     1   60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624  call xerrwv(
     1  'lsodi--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ',
     1   60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  call xerrwv(
     1  'lsodi--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ',
     1   60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  call xerrwv('lsodi--  at start of problem, too much accuracy   ',
     1   50, 26, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      requested for precision of machine..  see tolsf (=r1) ',
     1   60, 26, 0, 0, 0, 0, 1, tolsf, 0.0d0)
      rwork(14) = tolsf
      go to 700
 627  call xerrwv('lsodi--  trouble from intdy. itask = i1, tout = r1',
     1   50, 27, 0, 1, itask, 0, 1, tout, 0.0d0)
c
 700  if (illin .eq. 5) go to 710
      illin = illin + 1
      istate = -3
      return
 710  call xerrwv('lsodi--  repeated occurrences of illegal input    ',
     1   50, 302, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
c
 800  call xerrwv('lsodi--  run aborted.. apparent infinite loop     ',
     1   50, 303, 2, 0, 0, 0, 0, 0.0d0, 0.0d0)
      return
c----------------------- end of subroutine lsodi -----------------------
      end
