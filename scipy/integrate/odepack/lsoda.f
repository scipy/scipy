      subroutine lsoda (f, neq, y, t, tout, itol, rtol, atol, itask,
     1            istate, iopt, rwork, lrw, iwork, liw, jac, jt)
      external f, jac
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, jt, isav
      double precision y, t, tout, rtol, atol, rwork, rsav
      dimension neq(1), y(1), rtol(1), atol(1), rwork(lrw), iwork(liw)
      dimension rsav(240), isav(50)
c-----------------------------------------------------------------------
c this is the 24 feb 1997 version of
c lsoda.. livermore solver for ordinary differential equations, with
c         automatic method switching for stiff and nonstiff problems.
c
c this version is in double precision.
c
c lsoda solves the initial value problem for stiff or nonstiff
c systems of first order ode-s,
c     dy/dt = f(t,y) ,  or, in component form,
c     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
c
c this a variant version of the lsode package.
c it switches automatically between stiff and nonstiff methods.
c this means that the user does not have to determine whether the
c problem is stiff or not, and the solver will automatically choose the
c appropriate method.  it always starts with the nonstiff method.
c
c authors..
c                linda r. petzold  and  alan c. hindmarsh,
c                computing and mathematics research division, l-316
c                lawrence livermore national laboratory
c                livermore, ca 94550.
c
c references..
c 1.  alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c 2.  linda r. petzold, automatic selection of methods for solving
c     stiff and nonstiff systems of ordinary differential equations,
c     siam j. sci. stat. comput. 4 (1983), pp. 136-148.
c-----------------------------------------------------------------------
c summary of usage.
c
c communication between the user and the lsoda package, for normal
c situations, is summarized here.  this summary describes only a subset
c of the full set of options available.  see the full description for
c details, including alternative treatment of the jacobian matrix,
c optional inputs and outputs, nonstandard options, and
c instructions for special situations.  see also the example
c problem (with program and output) following this summary.
c
c a. first provide a subroutine of the form..
c               subroutine f (neq, t, y, ydot)
c               dimension y(neq), ydot(neq)
c which supplies the vector function f by loading ydot(i) with f(i).
c
c b. write a main program which calls subroutine lsoda once for
c each point at which answers are desired.  this should also provide
c for possible use of logical unit 6 for output of error messages
c by lsoda.  on the first call to lsoda, supply arguments as follows..
c f      = name of subroutine for right-hand side vector f.
c          this name must be declared external in calling program.
c neq    = number of first order ode-s.
c y      = array of initial values, of length neq.
c t      = the initial value of the independent variable.
c tout   = first point where output is desired (.ne. t).
c itol   = 1 or 2 according as atol (below) is a scalar or array.
c rtol   = relative tolerance parameter (scalar).
c atol   = absolute tolerance parameter (scalar or array).
c          the estimated local error in y(i) will be controlled so as
c          to be less than
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
c istate = integer flag (input and output).  set istate = 1.
c iopt   = 0 to indicate no optional inputs used.
c rwork  = real work array of length at least..
c             22 + neq * max(16, neq + 9).
c          see also paragraph e below.
c lrw    = declared length of rwork (in user-s dimension).
c iwork  = integer work array of length at least  20 + neq.
c liw    = declared length of iwork (in user-s dimension).
c jac    = name of subroutine for jacobian matrix.
c          use a dummy name.  see also paragraph e below.
c jt     = jacobian type indicator.  set jt = 2.
c          see also paragraph e below.
c note that the main program must declare arrays y, rwork, iwork,
c and possibly atol.
c
c c. the output from the first call (or any call) is..
c      y = array of computed values of y(t) vector.
c      t = corresponding value of independent variable (normally tout).
c istate = 2  if lsoda was successful, negative otherwise.
c          -1 means excess work done on this call (perhaps wrong jt).
c          -2 means excess accuracy requested (tolerances too small).
c          -3 means illegal input detected (see printed message).
c          -4 means repeated error test failures (check all inputs).
c          -5 means repeated convergence failures (perhaps bad jacobian
c             supplied or wrong choice of jt or tolerances).
c          -6 means error weight became zero during problem. (solution
c             component i vanished, and atol or atol(i) = 0.)
c          -7 means work space insufficient to finish (see messages).
c
c d. to continue the integration after a successful return, simply
c reset tout and call lsoda again.  no other parameters need be reset.
c
c e. note.. if and when lsoda regards the problem as stiff, and
c switches methods accordingly, it must make use of the neq by neq
c jacobian matrix, j = df/dy.  for the sake of simplicity, the
c inputs to lsoda recommended in paragraph b above cause lsoda to
c treat j as a full matrix, and to approximate it internally by
c difference quotients.  alternatively, j can be treated as a band
c matrix (with great potential reduction in the size of the rwork
c array).  also, in either the full or banded case, the user can supply
c j in closed form, with a routine whose name is passed as the jac
c argument.  these alternatives are described in the paragraphs on
c rwork, jac, and jt in the full description of the call sequence below.
c
c-----------------------------------------------------------------------
c example problem.
c
c the following is a simple example problem, with the coding
c needed for its solution by lsoda.  the problem is from chemical
c kinetics, and consists of the following three rate equations..
c     dy1/dt = -.04*y1 + 1.e4*y2*y3
c     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
c     dy3/dt = 3.e7*y2**2
c on the interval from t = 0.0 to t = 4.e10, with initial conditions
c y1 = 1.0, y2 = y3 = 0.  the problem is stiff.
c
c the following coding solves this problem with lsoda,
c printing results at t = .4, 4., ..., 4.e10.  it uses
c itol = 2 and atol much smaller for y2 than y1 or y3 because
c y2 has much smaller values.
c at the end of the run, statistical quantities of interest are
c printed (see optional outputs in the full description below).
c
c     external fex
c     double precision atol, rtol, rwork, t, tout, y
c     dimension y(3), atol(3), rwork(70), iwork(23)
c     neq = 3
c     y(1) = 1.0d0
c     y(2) = 0.0d0
c     y(3) = 0.0d0
c     t = 0.0d0
c     tout = 0.4d0
c     itol = 2
c     rtol = 1.0d-4
c     atol(1) = 1.0d-6
c     atol(2) = 1.0d-10
c     atol(3) = 1.0d-6
c     itask = 1
c     istate = 1
c     iopt = 0
c     lrw = 70
c     liw = 23
c     jt = 2
c     do 40 iout = 1,12
c       call lsoda(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,
c    1     iopt,rwork,lrw,iwork,liw,jdum,jt)
c       write(6,20)t,y(1),y(2),y(3)
c 20    format(' at t =',e12.4,'   y =',3e14.6)
c       if (istate .lt. 0) go to 80
c 40    tout = tout*10.0d0
c     write(6,60)iwork(11),iwork(12),iwork(13),iwork(19),rwork(15)
c 60  format(/' no. steps =',i4,'  no. f-s =',i4,'  no. j-s =',i4/
c    1   ' method last used =',i2,'   last switch was at t =',e12.4)
c     stop
c 80  write(6,90)istate
c 90  format(///' error halt.. istate =',i3)
c     stop
c     end
c
c     subroutine fex (neq, t, y, ydot)
c     double precision t, y, ydot
c     dimension y(3), ydot(3)
c     ydot(1) = -.04d0*y(1) + 1.0d4*y(2)*y(3)
c     ydot(3) = 3.0d7*y(2)*y(2)
c     ydot(2) = -ydot(1) - ydot(3)
c     return
c     end
c
c the output of this program (on a cdc-7600 in single precision)
c is as follows..
c
c   at t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02
c   at t =  4.0000e+00   y =  9.055333e-01  2.240655e-05  9.444430e-02
c   at t =  4.0000e+01   y =  7.158403e-01  9.186334e-06  2.841505e-01
c   at t =  4.0000e+02   y =  4.505250e-01  3.222964e-06  5.494717e-01
c   at t =  4.0000e+03   y =  1.831975e-01  8.941774e-07  8.168016e-01
c   at t =  4.0000e+04   y =  3.898730e-02  1.621940e-07  9.610125e-01
c   at t =  4.0000e+05   y =  4.936363e-03  1.984221e-08  9.950636e-01
c   at t =  4.0000e+06   y =  5.161831e-04  2.065786e-09  9.994838e-01
c   at t =  4.0000e+07   y =  5.179817e-05  2.072032e-10  9.999482e-01
c   at t =  4.0000e+08   y =  5.283401e-06  2.113371e-11  9.999947e-01
c   at t =  4.0000e+09   y =  4.659031e-07  1.863613e-12  9.999995e-01
c   at t =  4.0000e+10   y =  1.404280e-08  5.617126e-14  1.000000e+00
c
c   no. steps = 361  no. f-s = 693  no. j-s =  64
c   method last used = 2   last switch was at t =  6.0092e-03
c-----------------------------------------------------------------------
c full description of user interface to lsoda.
c
c the user interface to lsoda consists of the following parts.
c
c i.   the call sequence to subroutine lsoda, which is a driver
c      routine for the solver.  this includes descriptions of both
c      the call sequence arguments and of user-supplied routines.
c      following these descriptions is a description of
c      optional inputs available through the call sequence, and then
c      a description of optional outputs (in the work arrays).
c
c ii.  descriptions of other routines in the lsoda package that may be
c      (optionally) called by the user.  these provide the ability to
c      alter error message handling, save and restore the internal
c      common, and obtain specified derivatives of the solution y(t).
c
c iii. descriptions of common blocks to be declared in overlay
c      or similar environments, or to be saved when doing an interrupt
c      of the problem and continued solution later.
c
c iv.  description of a subroutine in the lsoda package,
c      which the user may replace with his own version, if desired.
c      this relates to the measurement of errors.
c
c-----------------------------------------------------------------------
c part i.  call sequence.
c
c the call sequence parameters used for input only are
c     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, jt,
c and those used for both input and output are
c     y, t, istate.
c the work arrays rwork and iwork are also used for conditional and
c optional inputs and optional outputs.  (the term output here refers
c to the return from subroutine lsoda to the user-s calling program.)
c
c the legality of input parameters will be thoroughly checked on the
c initial call for the problem, but not checked thereafter unless a
c change in input parameters is flagged by istate = 3 on input.
c
c the descriptions of the call arguments are as follows.
c
c f      = the name of the user-supplied subroutine defining the
c          ode system.  the system must be put in the first-order
c          form dy/dt = f(t,y), where f is a vector-valued function
c          of the scalar t and the vector y.  subroutine f is to
c          compute the function f.  it is to have the form
c               subroutine f (neq, t, y, ydot)
c               dimension y(1), ydot(1)
c          where neq, t, and y are input, and the array ydot = f(t,y)
c          is output.  y and ydot are arrays of length neq.
c          (in the dimension statement above, 1 is a dummy
c          dimension.. it can be replaced by any value.)
c          subroutine f should not alter y(1),...,y(neq).
c          f must be declared external in the calling program.
c
c          subroutine f may access user-defined quantities in
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array
c          (dimensioned in f) and/or y has length exceeding neq(1).
c          see the descriptions of neq and y below.
c
c          if quantities computed in the f routine are needed
c          externally to lsoda, an extra call to f should be made
c          for this purpose, for consistent and accurate results.
c          if only the derivative dy/dt is needed, use intdy instead.
c
c neq    = the size of the ode system (number of first order
c          ordinary differential equations).  used only for input.
c          neq may be decreased, but not increased, during the problem.
c          if neq is decreased (with istate = 3 on input), the
c          remaining components of y should be left undisturbed, if
c          these are to be accessed in f and/or jac.
c
c          normally, neq is a scalar, and it is generally referred to
c          as a scalar in this user interface description.  however,
c          neq may be an array, with neq(1) set to the system size.
c          (the lsoda package accesses only neq(1).)  in either case,
c          this parameter is passed as the neq argument in all calls
c          to f and jac.  hence, if it is an array, locations
c          neq(2),... may be used to store other integer data and pass
c          it to f and/or jac.  subroutines f and/or jac must include
c          neq in a dimension statement in that case.
c
c y      = a real array for the vector of dependent variables, of
c          length neq or more.  used for both input and output on the
c          first call (istate = 1), and only for output on other calls.
c          on the first call, y must contain the vector of initial
c          values.  on output, y contains the computed solution vector,
c          evaluated at t.  if desired, the y array may be used
c          for other purposes between calls to the solver.
c
c          this array is passed as the y argument in all calls to
c          f and jac.  hence its length may exceed neq, and locations
c          y(neq+1),... may be used to store other real data and
c          pass it to f and/or jac.  (the lsoda package accesses only
c          y(1),...,y(neq).)
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
c          when starting the problem (istate = 1), tout may be equal
c          to t for one call, then should .ne. t for the next call.
c          for the initial t, an input value of tout .ne. t is used
c          in order to determine the direction of the integration
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
c                      max-norm of ( e(i)/ewt(i) )   .le.   1,
c          where ewt = (ewt(i)) is a vector of positive error weights.
c          the values of rtol and atol should all be non-negative.
c          the following table gives the types (scalar/array) of
c          rtol and atol, and the corresponding form of ewt(i).
c
c             itol    rtol       atol          ewt(i)
c              1     scalar     scalar     rtol*abs(y(i)) + atol
c              2     scalar     array      rtol*abs(y(i)) + atol(i)
c              3     array      scalar     rtol(i)*abs(y(i)) + atol
c              4     array      array      rtol(i)*abs(y(i)) + atol(i)
c
c          when either of these parameters is a scalar, it need not
c          be dimensioned in the user-s calling program.
c
c          if none of the above choices (with itol, rtol, and atol
c          fixed throughout the problem) is suitable, more general
c          error controls can be obtained by substituting a
c          user-supplied routine for the setting of ewt.
c          see part iv below.
c
c          if global errors are to be estimated by making a repeated
c          run on the same problem with smaller tolerances, then all
c          components of rtol and atol (i.e. of ewt) should be scaled
c          down uniformly.
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
c          the state of the calculation.
c
c          on input, the values of istate are as follows.
c          1  means this is the first call for the problem
c             (initializations will be done).  see note below.
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
c             neq, itol, rtol, atol, iopt, lrw, liw, jt, ml, mu,
c             and any optional inputs except h0, mxordn, and mxords.
c             (see iwork description for ml and mu.)
c          note..  a preliminary call with tout = t is not counted
c          as a first call here, as no initialization or checking of
c          input is done.  (such a call is sometimes useful for the
c          purpose of outputting the initial conditions.)
c          thus the first call for which tout .ne. t requires
c          istate = 1 on input.
c
c          on output, istate has the following values and meanings.
c           1  means nothing was done, as tout was equal to t with
c              istate = 1 on input.  (however, an internal counter was
c              set to detect and prevent repeated calls of this type.)
c           2  means the integration was performed successfully.
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
c              this may be caused by an inaccurate jacobian matrix,
c              if one is being used.
c          -6  means ewt(i) became zero for some i during the
c              integration.  pure relative error control (atol(i)=0.0)
c              was requested on a variable which has now vanished.
c              the integration was successful as far as t.
c          -7  means the length of rwork and/or iwork was too small to
c              proceed, but the integration was successful as far as t.
c              this happens when lsoda chooses to switch methods
c              but lrw and/or liw is too small for the new method.
c
c          note..  since the normal output value of istate is 2,
c          it does not need to be reset for normal continuation.
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
c rwork  = a real array (double precision) for work space, and (in the
c          first 20 words) for conditional and optional inputs and
c          optional outputs.
c          as lsoda switches automatically between stiff and nonstiff
c          methods, the required length of rwork can change during the
c          problem.  thus the rwork array passed to lsoda can either
c          have a static (fixed) length large enough for both methods,
c          or have a dynamic (changing) length altered by the calling
c          program in response to output from lsoda.
c
c                       --- fixed length case ---
c          if the rwork length is to be fixed, it should be at least
c               max (lrn, lrs),
c          where lrn and lrs are the rwork lengths required when the
c          current method is nonstiff or stiff, respectively.
c
c          the separate rwork length requirements lrn and lrs are
c          as follows..
c          if neq is constant and the maximum method orders have
c          their default values, then
c             lrn = 20 + 16*neq,
c             lrs = 22 + 9*neq + neq**2           if jt = 1 or 2,
c             lrs = 22 + 10*neq + (2*ml+mu)*neq   if jt = 4 or 5.
c          under any other conditions, lrn and lrs are given by..
c             lrn = 20 + nyh*(mxordn+1) + 3*neq,
c             lrs = 20 + nyh*(mxords+1) + 3*neq + lmat,
c          where
c             nyh    = the initial value of neq,
c             mxordn = 12, unless a smaller  value is given as an
c                      optional input,
c             mxords = 5, unless a smaller value is given as an
c                      optional input,
c             lmat   = length of matrix work space..
c             lmat   = neq**2 + 2              if jt = 1 or 2,
c             lmat   = (2*ml + mu + 1)*neq + 2 if jt = 4 or 5.
c
c                       --- dynamic length case ---
c          if the length of rwork is to be dynamic, then it should
c          be at least lrn or lrs, as defined above, depending on the
c          current method.  initially, it must be at least lrn (since
c          lsoda starts with the nonstiff method).  on any return
c          from lsoda, the optional output mcur indicates the current
c          method.  if mcur differs from the value it had on the
c          previous return, or if there has only been one call to
c          lsoda and mcur is now 2, then lsoda has switched
c          methods during the last call, and the length of rwork
c          should be reset (to lrn if mcur = 1, or to lrs if
c          mcur = 2).  (an increase in the rwork length is required
c          if lsoda returned istate = -7, but not otherwise.)
c          after resetting the length, call lsoda with istate = 3
c          to signal that change.
c
c lrw    = the length of the array rwork, as declared by the user.
c          (this will be checked by the solver.)
c
c iwork  = an integer array for work space.
c          as lsoda switches automatically between stiff and nonstiff
c          methods, the required length of iwork can change during
c          problem, between
c             lis = 20 + neq   and   lin = 20,
c          respectively.  thus the iwork array passed to lsoda can
c          either have a fixed length of at least 20 + neq, or have a
c          dynamic length of at least lin or lis, depending on the
c          current method.  the comments on dynamic length under
c          rwork above apply here.  initially, this length need
c          only be at least lin = 20.
c
c          the first few words of iwork are used for conditional and
c          optional inputs and optional outputs.
c 
c          the following 2 words in iwork are conditional inputs..
c            iwork(1) = ml     these are the lower and upper
c            iwork(2) = mu     half-bandwidths, respectively, of the
c                       banded jacobian, excluding the main diagonal.
c                       the band is defined by the matrix locations
c                       (i,j) with i-ml .le. j .le. i+mu.  ml and mu
c                       must satisfy  0 .le.  ml,mu  .le. neq-1.
c                       these are required if jt is 4 or 5, and
c                       ignored otherwise.  ml and mu may in fact be
c                       the band parameters for a matrix to which
c                       df/dy is only approximately equal.
c
c liw    = the length of the array iwork, as declared by the user.
c          (this will be checked by the solver.)
c
c note.. the base addresses of the work arrays must not be
c altered between calls to lsoda for the same problem.
c the contents of the work arrays must not be altered
c between calls, except possibly for the conditional and
c optional inputs, and except for the last 3*neq words of rwork.
c the latter space is used for internal scratch space, and so is
c available for use by the user outside lsoda between calls, if
c desired (but not for use by f or jac).
c
c jac    = the name of the user-supplied routine to compute the
c          jacobian matrix, df/dy, if jt = 1 or 4.  the jac routine
c          is optional, but if the problem is expected to be stiff much
c          of the time, you are encouraged to supply jac, for the sake
c          of efficiency.  (alternatively, set jt = 2 or 5 to have
c          lsoda compute df/dy internally by difference quotients.)
c          if and when lsoda uses df/dy, if treats this neq by neq
c          matrix either as full (jt = 1 or 2), or as banded (jt =
c          4 or 5) with half-bandwidths ml and mu (discussed under
c          iwork above).  in either case, if jt = 1 or 4, the jac
c          routine must compute df/dy as a function of the scalar t
c          and the vector y.  it is to have the form
c               subroutine jac (neq, t, y, ml, mu, pd, nrowpd)
c               dimension y(1), pd(nrowpd,1)
c          where neq, t, y, ml, mu, and nrowpd are input and the array
c          pd is to be loaded with partial derivatives (elements of
c          the jacobian matrix) on output.  pd must be given a first
c          dimension of nrowpd.  t and y have the same meaning as in
c          subroutine f.  (in the dimension statement above, 1 is a
c          dummy dimension.. it can be replaced by any value.)
c               in the full matrix case (jt = 1), ml and mu are
c          ignored, and the jacobian is to be loaded into pd in
c          columnwise manner, with df(i)/dy(j) loaded into pd(i,j).
c               in the band matrix case (jt = 4), the elements
c          within the band are to be loaded into pd in columnwise
c          manner, with diagonal lines of df/dy loaded into the rows
c          of pd.  thus df(i)/dy(j) is to be loaded into pd(i-j+mu+1,j).
c          ml and mu are the half-bandwidth parameters (see iwork).
c          the locations in pd in the two triangular areas which
c          correspond to nonexistent matrix elements can be ignored
c          or loaded arbitrarily, as they are overwritten by lsoda.
c               jac need not provide df/dy exactly.  a crude
c          approximation (possibly with a smaller bandwidth) will do.
c               in either case, pd is preset to zero by the solver,
c          so that only the nonzero elements need be loaded by jac.
c          each call to jac is preceded by a call to f with the same
c          arguments neq, t, and y.  thus to gain some efficiency,
c          intermediate quantities shared by both calculations may be
c          saved in a user common block by f and not recomputed by jac,
c          if desired.  also, jac may alter the y array, if desired.
c          jac must be declared external in the calling program.
c               subroutine jac may access user-defined quantities in
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array
c          (dimensioned in jac) and/or y has length exceeding neq(1).
c          see the descriptions of neq and y above.
c
c jt     = jacobian type indicator.  used only for input.
c          jt specifies how the jacobian matrix df/dy will be
c          treated, if and when lsoda requires this matrix.
c          jt has the following values and meanings..
c           1 means a user-supplied full (neq by neq) jacobian.
c           2 means an internally generated (difference quotient) full
c             jacobian (using neq extra calls to f per df/dy value).
c           4 means a user-supplied banded jacobian.
c           5 means an internally generated banded jacobian (using
c             ml+mu+1 extra calls to f per df/dy evaluation).
c          if jt = 1 or 4, the user must supply a subroutine jac
c          (the name is arbitrary) as described above under jac.
c          if jt = 2 or 5, a dummy argument can be used.
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
c ixpr    iwork(5)  flag to generate extra printing at method switches.
c                   ixpr = 0 means no extra printing (the default).
c                   ixpr = 1 means print data on each switch.
c                   t, h, and nst will be printed on the same logical
c                   unit as used for error messages.
c
c mxstep  iwork(6)  maximum number of (internally defined) steps
c                   allowed during one call to the solver.
c                   the default value is 500.
c
c mxhnil  iwork(7)  maximum number of messages printed (per problem)
c                   warning that t + h = t on a step (h = step size).
c                   this must be positive to result in a non-default
c                   value.  the default value is 10.
c
c mxordn  iwork(8)  the maximum order to be allowed for the nonstiff
c                   (adams) method.  the default value is 12.
c                   if mxordn exceeds the default value, it will
c                   be reduced to the default value.
c                   mxordn is held constant during the problem.
c
c mxords  iwork(9)  the maximum order to be allowed for the stiff
c                   (bdf) method.  the default value is 5.
c                   if mxords exceeds the default value, it will
c                   be reduced to the default value.
c                   mxords is held constant during the problem.
c-----------------------------------------------------------------------
c optional outputs.
c
c as optional additional output from lsoda, the variables listed
c below are quantities related to the performance of lsoda
c which are available to the user.  these are communicated by way of
c the work arrays, but also have internal mnemonic names as shown.
c except where stated otherwise, all of these outputs are defined
c on any successful return from lsoda, and on any return with
c istate = -1, -2, -4, -5, or -6.  on an illegal input return
c (istate = -3), they will be unchanged from their existing values
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
c tsw     rwork(15) the value of t at the time of the last method
c                   switch, if any.
c
c nst     iwork(11) the number of steps taken for the problem so far.
c
c nfe     iwork(12) the number of f evaluations for the problem so far.
c
c nje     iwork(13) the number of jacobian evaluations (and of matrix
c                   lu decompositions) for the problem so far.
c
c nqu     iwork(14) the method order last used (successfully).
c
c nqcur   iwork(15) the order to be attempted on the next step.
c
c imxer   iwork(16) the index of the component of largest magnitude in
c                   the weighted local error vector ( e(i)/ewt(i) ),
c                   on an error return with istate = -4 or -5.
c
c lenrw   iwork(17) the length of rwork actually required, assuming
c                   that the length of rwork is to be fixed for the
c                   rest of the problem, and that switching may occur.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c leniw   iwork(18) the length of iwork actually required, assuming
c                   that the length of iwork is to be fixed for the
c                   rest of the problem, and that switching may occur.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c mused   iwork(19) the method indicator for the last successful step..
c                   1 means adams (nonstiff), 2 means bdf (stiff).
c
c mcur    iwork(20) the current method indicator..
c                   1 means adams (nonstiff), 2 means bdf (stiff).
c                   this is the method to be attempted
c                   on the next step.  thus it differs from mused
c                   only if a method switch has just been made.
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
c acor     lacor         array of size neq used for the accumulated
c         (from common   corrections on each step, scaled on output
c           as noted)    to represent the estimated local error in y
c                        on the last step.  this is the vector e in
c                        the description of the error control.  it is
c                        defined only on a successful return from lsoda.
c                        the base address lacor is obtained by
c                        including in the user-s program the
c                        following 3 lines..
c                           double precision rls
c                           common /ls0001/ rls(218), ils(39)
c                           lacor = ils(5)
c
c-----------------------------------------------------------------------
c part ii.  other routines callable.
c
c the following are optional calls which the user may make to
c gain additional capabilities in conjunction with lsoda.
c (the routines xsetun and xsetf are designed to conform to the
c slatec error handling package.)
c
c     form of call                  function
c   call xsetun(lun)          set the logical unit number, lun, for
c                             output of messages from lsoda, if
c                             the default is not desired.
c                             the default value of lun is 6.
c
c   call xsetf(mflag)         set a flag to control the printing of
c                             messages by lsoda.
c                             mflag = 0 means do not print. (danger..
c                             this risks losing valuable information.)
c                             mflag = 1 means print (the default).
c
c                             either of the above calls may be made at
c                             any time and will take effect immediately.
c
c   call srcma(rsav,isav,job) saves and restores the contents of
c                             the internal common blocks used by
c                             lsoda (see part iii below).
c                             rsav must be a real array of length 240
c                             or more, and isav must be an integer
c                             array of length 50 or more.
c                             job=1 means save common into rsav/isav.
c                             job=2 means restore common from rsav/isav.
c                                srcma is useful if one is
c                             interrupting a run and restarting
c                             later, or alternating between two or
c                             more problems solved with lsoda.
c
c   call intdy(,,,,,)         provide derivatives of y, of various
c        (see below)          orders, at a specified point t, if
c                             desired.  it may be called only after
c                             a successful return from lsoda.
c
c the detailed instructions for using intdy are as follows.
c the form of the call is..
c
c   call intdy (t, k, rwork(21), nyh, dky, iflag)
c
c the input parameters are..
c
c t         = value of independent variable where answers are desired
c             (normally the same as the t last returned by lsoda).
c             for valid results, t must lie between tcur - hu and tcur.
c             (see optional outputs for tcur and hu.)
c k         = integer order of the derivative desired.  k must satisfy
c             0 .le. k .le. nqcur, where nqcur is the current order
c             (see optional outputs).  the capability corresponding
c             to k = 0, i.e. computing y(t), is already provided
c             by lsoda directly.  since nqcur .ge. 1, the first
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
c if lsoda is to be used in an overlay situation, the user
c must declare, in the primary overlay, the variables in..
c   (1) the call sequence to lsoda,
c   (2) the three internal common blocks
c         /ls0001/  of length  257  (218 double precision words
c                         followed by 39 integer words),
c         /lsa001/  of length  31    (22 double precision words
c                         followed by  9 integer words),
c         /eh0001/  of length  2 (integer words).
c
c if lsoda is used on a system in which the contents of internal
c common blocks are not preserved between calls, the user should
c declare the above common blocks in his main program to insure
c that their contents are preserved.
c
c if the solution of a given problem by lsoda is to be interrupted
c and then later continued, such as when restarting an interrupted run
c or alternating between two or more problems, the user should save,
c following the return from the last lsoda call prior to the
c interruption, the contents of the call sequence variables and the
c internal common blocks, and later restore these values before the
c next lsoda call for that problem.  to save and restore the common
c blocks, use subroutine srcma (see part ii above).
c
c-----------------------------------------------------------------------
c part iv.  optionally replaceable solver routines.
c
c below is a description of a routine in the lsoda package which
c relates to the measurement of errors, and can be
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
c where neq, itol, rtol, and atol are as in the lsoda call sequence,
c ycur contains the current dependent variable vector, and
c ewt is the array of weights set by ewset.
c
c if the user supplies this subroutine, it must return in ewt(i)
c (i = 1,...,neq) a positive quantity suitable for comparing errors
c in y(i) to.  the ewt array returned by ewset is passed to the
c vmnorm routine, and also used by lsoda in the computation
c of the optional output imxer, and the increments for difference
c quotient jacobians.
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
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c other routines in the lsoda package.
c
c in addition to subroutine lsoda, the lsoda package includes the
c following subroutines and function routines..
c  intdy    computes an interpolated value of the y vector at t = tout.
c  stoda    is the core integrator, which does one step of the
c           integration and the associated error control.
c  cfode    sets all method coefficients and test constants.
c  prja     computes and preprocesses the jacobian matrix j = df/dy
c           and the newton iteration matrix p = i - h*l0*j.
c  solsy    manages solution of linear system in chord iteration.
c  ewset    sets the error weight vector ewt before each step.
c  vmnorm   computes the weighted max-norm of a vector.
c  fnorm    computes the norm of a full matrix consistent with the
c           weighted max-norm on vectors.
c  bnorm    computes the norm of a band matrix consistent with the
c           weighted max-norm on vectors.
c  srcma    is a user-callable routine to save and restore
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
c note..  vmnorm, fnorm, bnorm, idamax, ddot, and d1mach are function
c routines.  all the others are subroutines.
c
c the intrinsic and external routines used by lsoda are..
c dabs, dmax1, dmin1, dfloat, max0, min0, mod, dsign, dsqrt, and write.
c
c a block data subprogram is also included with the package,
c for loading some of the variables in internal common.
c
c-----------------------------------------------------------------------
c the following card is for optimized compilation on lll compilers.
clll. optimize
c-----------------------------------------------------------------------
      external prja, solsy
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer insufr, insufi, ixpr, iowns2, jtyp, mused, mxordn, mxords
      integer i, i1, i2, iflag, imxer, kgo, lf0,
     1   leniw, lenrw, lenwm, ml, mord, mu, mxhnl0, mxstp0
      integer len1, len1c, len1n, len1s, len2, leniwc,
     1   lenrwc, lenrwn, lenrws
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision tsw, rowns2, pdnorm
      double precision atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0,
     2   d1mach, vmnorm
      dimension mord(2)
      logical ihit
c-----------------------------------------------------------------------
c the following two internal common blocks contain
c (a) variables which are local to any subroutine but whose values must
c     be preserved between calls to the routine (own variables), and
c (b) variables which are communicated between subroutines.
c the structure of each block is as follows..  all real variables are
c listed first, followed by all integers.  within each type, the
c variables are grouped with those local to subroutine lsoda first,
c then those local to subroutine stoda, and finally those used
c for communication.  the block ls0001 is declared in subroutines
c lsoda, intdy, stoda, prja, and solsy.  the block lsa001 is declared
c in subroutines lsoda, stoda, and prja.  groups of variables are
c replaced by dummy arrays in the common declarations in routines
c where those variables are not used.
c-----------------------------------------------------------------------
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lsa001/ tsw, rowns2(20), pdnorm,
     1   insufr, insufi, ixpr, iowns2(2), jtyp, mused, mxordn, mxords
c
      data mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
c-----------------------------------------------------------------------
c block a.
c this code block is executed on every call.
c it tests istate and itask for legality and branches appropriately.
c if istate .gt. 1 but the flag init shows that initialization has
c not yet been done, an error return occurs.
c if istate = 1 and tout = t, jump to block g and return immediately.
c-----------------------------------------------------------------------
      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .eq. 1) go to 10
      if (init .eq. 0) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   init = 0
      if (tout .eq. t) go to 430
 20   ntrep = 0
c-----------------------------------------------------------------------
c block b.
c the next code block is executed for the initial call (istate = 1),
c or for a continuation call with parameter changes (istate = 3).
c it contains checking of all inputs and various initializations.
c
c first check legality of the non-optional inputs neq, itol, iopt,
c jt, ml, and mu.
c-----------------------------------------------------------------------
      if (neq(1) .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      if (jt .eq. 3 .or. jt .lt. 1 .or. jt .gt. 5) go to 608
      jtyp = jt
      if (jt .le. 2) go to 30
      ml = iwork(1)
      mu = iwork(2)
      if (ml .lt. 0 .or. ml .ge. n) go to 609
      if (mu .lt. 0 .or. mu .ge. n) go to 610
 30   continue
c next process and check the optional inputs. --------------------------
      if (iopt .eq. 1) go to 40
      ixpr = 0
      mxstep = mxstp0
      mxhnil = mxhnl0
      hmxi = 0.0d0
      hmin = 0.0d0
      if (istate .ne. 1) go to 60
      h0 = 0.0d0
      mxordn = mord(1)
      mxords = mord(2)
      go to 60
 40   ixpr = iwork(5)
      if (ixpr .lt. 0 .or. ixpr .gt. 1) go to 611
      mxstep = iwork(6)
      if (mxstep .lt. 0) go to 612
      if (mxstep .eq. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      if (mxhnil .lt. 0) go to 613
      if (mxhnil .eq. 0) mxhnil = mxhnl0
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      mxordn = iwork(8)
      if (mxordn .lt. 0) go to 628
      if (mxordn .eq. 0) mxordn = 100
      mxordn = min0(mxordn,mord(1))
      mxords = iwork(9)
      if (mxords .lt. 0) go to 629
      if (mxords .eq. 0) mxords = 100
      mxords = min0(mxords,mord(2))
      if ((tout - t)*h0 .lt. 0.0d0) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. 0.0d0) go to 615
      hmxi = 0.0d0
      if (hmax .gt. 0.0d0) hmxi = 1.0d0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0d0) go to 616
c-----------------------------------------------------------------------
c set work array pointers and check lengths lrw and liw.
c if istate = 1, meth is initialized to 1 here to facilitate the
c checking of work space lengths.
c pointers to segments of rwork and iwork are named by prefixing l to
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).
c segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
c if the lengths provided are insufficient for the current method,
c an error return occurs.  this is treated as illegal input on the
c first call, but as a problem interruption with istate = -7 on a
c continuation call.  if the lengths are sufficient for the current
c method but not for both methods, a warning message is sent.
c-----------------------------------------------------------------------
 60   if (istate .eq. 1) meth = 1
      if (istate .eq. 1) nyh = n
      lyh = 21
      len1n = 20 + (mxordn + 1)*nyh
      len1s = 20 + (mxords + 1)*nyh
      lwm = len1s + 1
      if (jt .le. 2) lenwm = n*n + 2
      if (jt .ge. 4) lenwm = (2*ml + mu + 1)*n + 2
      len1s = len1s + lenwm
      len1c = len1n
      if (meth .eq. 2) len1c = len1s
      len1 = max0(len1n,len1s)
      len2 = 3*n
      lenrw = len1 + len2
      lenrwn = len1n + len2
      lenrws = len1s + len2
      lenrwc = len1c + len2
      iwork(17) = lenrw
      liwm = 1
      leniw = 20 + n
      leniwc = 20
      if (meth .eq. 2) leniwc = leniw
      iwork(18) = leniw
      if (istate .eq. 1 .and. lrw .lt. lenrwc) go to 617
      if (istate .eq. 1 .and. liw .lt. leniwc) go to 618
      if (istate .eq. 3 .and. lrw .lt. lenrwc) go to 550
      if (istate .eq. 3 .and. liw .lt. leniwc) go to 555
      lewt = len1 + 1
      insufr = 0
      if (lrw .ge. lenrw) go to 65
      insufr = 2
      lewt = len1c + 1
      call xerrwv(
     1  'lsoda--  warning.. rwork length is sufficient for now, but  ',
     1   60, 103, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      may not be later.  integration will proceed anyway.   ',
     1   60, 103, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      length needed is lenrw = i1, while lrw = i2.',
     1   50, 103, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
 65   lsavf = lewt + n
      lacor = lsavf + n
      insufi = 0
      if (liw .ge. leniw) go to 70
      insufi = 2
      call xerrwv(
     1  'lsoda--  warning.. iwork length is sufficient for now, but  ',
     1   60, 104, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      may not be later.  integration will proceed anyway.   ',
     1   60, 104, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      length needed is leniw = i1, while liw = i2.',
     1   50, 104, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
 70   continue
c check rtol and atol for legality. ------------------------------------
      rtoli = rtol(1)
      atoli = atol(1)
      do 75 i = 1,n
        if (itol .ge. 3) rtoli = rtol(i)
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        if (rtoli .lt. 0.0d0) go to 619
        if (atoli .lt. 0.0d0) go to 620
 75     continue
      if (istate .eq. 1) go to 100
c if istate = 3, set flag to signal parameter changes to stoda. --------
      jstart = -1
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
c the next block is for the initial call only (istate = 1).
c it contains all remaining initializations, the initial call to f,
c and the calculation of the initial step size.
c the error weights in ewt are inverted after being loaded.
c-----------------------------------------------------------------------
 100  uround = d1mach(4)
      tn = t
      tsw = t
      maxord = mxordn
      if (itask .ne. 4 .and. itask .ne. 5) go to 110
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. 0.0d0) go to 625
      if (h0 .ne. 0.0d0 .and. (t + h0 - tcrit)*h0 .gt. 0.0d0)
     1   h0 = tcrit - t
 110  jstart = 0
      nhnil = 0
      nst = 0
      nje = 0
      nslast = 0
      hu = 0.0d0
      nqu = 0
      mused = 0
      miter = 0
      ccmax = 0.3d0
      maxcor = 3
      msbp = 20
      mxncf = 10
c initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      call srcma(rsav, isav, 1)
      call f (neq, t, y, rwork(lf0))
      call srcma(rsav, isav, 2)
      nfe = 1
c load the initial value vector in yh. ---------------------------------
      do 115 i = 1,n
 115    rwork(i+lyh-1) = y(i)
c load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      nq = 1
      h = 1.0d0
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 120 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0d0) go to 621
 120    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
c-----------------------------------------------------------------------
c the coding below computes the step size, h0, to be attempted on the
c first step, unless the user has supplied a value for this.
c first check that tout - t differs significantly from zero.
c a scalar tolerance quantity tol is computed, as max(rtol(i))
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
c so as to be between 100*uround and 1.0e-3.
c then the computed value h0 is given by..
c
c   h0**(-2)  =  1./(tol * w0**2)  +  tol * (norm(f))**2
c
c where   w0     = max ( abs(t), abs(tout) ),
c         f      = the initial value of the vector f(t,y), and
c         norm() = the weighted vector norm used throughout, given by
c                  the vmnorm function routine, and weighted by the
c                  tolerances initially loaded into the ewt array.
c the sign of h0 is inferred from the initial values of tout and t.
c abs(h0) is made .le. abs(tout-t) in any case.
c-----------------------------------------------------------------------
      if (h0 .ne. 0.0d0) go to 180
      tdist = dabs(tout - t)
      w0 = dmax1(dabs(t),dabs(tout))
      if (tdist .lt. 2.0d0*uround*w0) go to 622
      tol = rtol(1)
      if (itol .le. 2) go to 140
      do 130 i = 1,n
 130    tol = dmax1(tol,rtol(i))
 140  if (tol .gt. 0.0d0) go to 160
      atoli = atol(1)
      do 150 i = 1,n
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        ayi = dabs(y(i))
        if (ayi .ne. 0.0d0) tol = dmax1(tol,atoli/ayi)
 150    continue
 160  tol = dmax1(tol,100.0d0*uround)
      tol = dmin1(tol,0.001d0)
      sum = vmnorm (n, rwork(lf0), rwork(lewt))
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
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)
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
      t = tn
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
      if (ihit) t = tcrit
      if (ihit) go to 400
      tnext = tn + h*(1.0d0 + 4.0d0*uround)
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
      if (istate .eq. 2 .and. jstart .ge. 0) jstart = -2
c-----------------------------------------------------------------------
c block e.
c the next block is normally executed for all calls and contains
c the call to the one-step core integrator stoda.
c
c this is a looping point for the integration steps.
c
c first check for too many steps being taken, update ewt (if not at
c start of problem), check for too much accuracy being requested, and
c check for h below the roundoff level in t.
c-----------------------------------------------------------------------
 250  continue
      if (meth .eq. mused) go to 255
      if (insufr .eq. 1) go to 550
      if (insufi .eq. 1) go to 555
 255  if ((nst-nslast) .ge. mxstep) go to 500
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 260 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0d0) go to 510
 260    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
 270  tolsf = uround*vmnorm (n, rwork(lyh), rwork(lewt))
      if (tolsf .le. 0.01d0) go to 280
      tolsf = tolsf*200.0d0
      if (nst .eq. 0) go to 626
      go to 520
 280  if ((tn + h) .ne. tn) go to 290
      nhnil = nhnil + 1
      if (nhnil .gt. mxhnil) go to 290
      call xerrwv('lsoda--  warning..internal t (=r1) and h (=r2) are',
     1   50, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      such that in the machine, t + h = t on the next step  ',
     1   60, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      (h = step size). solver will continue anyway',
     1   50, 101, 0, 0, 0, 0, 2, tn, h)
      if (nhnil .lt. mxhnil) go to 290
      call xerrwv('lsoda--  above warning has been issued i1 times.  ',
     1   50, 102, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      it will not be issued again for this problem',
     1   50, 102, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
 290  continue
c-----------------------------------------------------------------------
c     call stoda(neq,y,yh,nyh,yh,ewt,savf,acor,wm,iwm,f,jac,prja,solsy)
c-----------------------------------------------------------------------
      call stoda (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),
     1   rwork(lsavf), rwork(lacor), rwork(lwm), iwork(liwm),
     2   f, jac, prja, solsy)
      kgo = 1 - kflag
      go to (300, 530, 540), kgo
c-----------------------------------------------------------------------
c block f.
c the following block handles the case of a successful return from the
c core integrator (kflag = 0).
c if a method switch was just made, record tsw, reset maxord,
c set jstart to -1 to signal stoda to complete the switch,
c and do extra printing of data if ixpr = 1.
c then, in any case, check for stop conditions.
c-----------------------------------------------------------------------
 300  init = 1
      if (meth .eq. mused) go to 310
      tsw = tn
      maxord = mxordn
      if (meth .eq. 2) maxord = mxords
      if (meth .eq. 2) rwork(lwm) = dsqrt(uround)
      insufr = min0(insufr,1)
      insufi = min0(insufi,1)
      jstart = -1
      if (ixpr .eq. 0) go to 310
      if (meth .eq. 2) call xerrwv(
     1  'lsoda-- a switch to the bdf (stiff) method has occurred     ',
     1   60, 105, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      if (meth .eq. 1) call xerrwv(
     1  'lsoda-- a switch to the adams (nonstiff) method has occurred',
     1   60, 106, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '     at t = r1,  tentative step size h = r2,  step nst = i1 ',
     1   60, 107, 0, 1, nst, 0, 2, tn, h)
 310  go to (320, 400, 330, 340, 350), itask
c itask = 1.  if tout has been reached, interpolate. -------------------
 320  if ((tn - tout)*h .lt. 0.0d0) go to 250
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
      if (jstart .ge. 0) jstart = -2
      go to 250
c itask = 5.  see if tcrit was reached and jump to exit. ---------------
 350  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
c-----------------------------------------------------------------------
c block g.
c the following block handles all successful returns from lsoda.
c if itask .ne. 1, y is loaded from yh and t is set accordingly.
c istate is set to 2, the illegal input counter is zeroed, and the
c optional outputs are loaded into the work arrays before returning.
c if istate = 1 and tout = t, there is a return with no action taken,
c except that if this has happened repeatedly, the run is terminated.
c-----------------------------------------------------------------------
 400  do 410 i = 1,n
 410    y(i) = rwork(i+lyh-1)
      t = tn
      if (itask .ne. 4 .and. itask .ne. 5) go to 420
      if (ihit) t = tcrit
 420  istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      rwork(15) = tsw
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = mused
      iwork(20) = meth
      return
c
 430  ntrep = ntrep + 1
      if (ntrep .lt. 5) return
      call xerrwv(
     1  'lsoda--  repeated calls with istate = 1 and tout = t (=r1)  ',
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
 500  call xerrwv('lsoda--  at current t (=r1), mxstep (=i1) steps   ',
     1   50, 201, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      taken on this call before reaching tout     ',
     1   50, 201, 0, 1, mxstep, 0, 1, tn, 0.0d0)
      istate = -1
      go to 580
c ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  ewti = rwork(lewt+i-1)
      call xerrwv('lsoda--  at t (=r1), ewt(i1) has become r2 .le. 0.',
     1   50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 580
c too much accuracy requested for machine precision. -------------------
 520  call xerrwv('lsoda--  at t (=r1), too much accuracy requested  ',
     1   50, 203, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      for precision of machine..  see tolsf (=r2) ',
     1   50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580
c kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
 530  call xerrwv('lsoda--  at t(=r1) and step size h(=r2), the error',
     1   50, 204, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      test failed repeatedly or with abs(h) = hmin',
     1   50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 560
c kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
 540  call xerrwv('lsoda--  at t (=r1) and step size h (=r2), the    ',
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      corrector convergence failed repeatedly     ',
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      or with abs(h) = hmin   ',
     1   30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
      go to 560
c rwork length too small to proceed. -----------------------------------
 550  call xerrwv('lsoda--  at current t(=r1), rwork length too small',
     1   50, 206, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      to proceed.  the integration was otherwise successful.',
     1   60, 206, 0, 0, 0, 0, 1, tn, 0.0d0)
      istate = -7
      go to 580
c iwork length too small to proceed. -----------------------------------
 555  call xerrwv('lsoda--  at current t(=r1), iwork length too small',
     1   50, 207, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      to proceed.  the integration was otherwise successful.',
     1   60, 207, 0, 0, 0, 0, 1, tn, 0.0d0)
      istate = -7
      go to 580
c compute imxer if relevant. -------------------------------------------
 560  big = 0.0d0
      imxer = 1
      do 570 i = 1,n
        size = dabs(rwork(i+lacor-1)*rwork(i+lewt-1))
        if (big .ge. size) go to 570
        big = size
        imxer = i
 570    continue
      iwork(16) = imxer
c set y vector, t, illin, and optional outputs. ------------------------
 580  do 590 i = 1,n
 590    y(i) = rwork(i+lyh-1)
      t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      rwork(15) = tsw
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = mused
      iwork(20) = meth
      return
c-----------------------------------------------------------------------
c block i.
c the following block handles all error returns due to illegal input
c (istate = -3), as detected before calling the core integrator.
c first the error message routine is called.  then if there have been
c 5 consecutive such returns just before this call to the solver,
c the run is halted.
c-----------------------------------------------------------------------
 601  call xerrwv('lsoda--  istate (=i1) illegal ',
     1   30, 1, 0, 1, istate, 0, 0, 0.0d0, 0.0d0)
      go to 700
 602  call xerrwv('lsoda--  itask (=i1) illegal  ',
     1   30, 2, 0, 1, itask, 0, 0, 0.0d0, 0.0d0)
      go to 700
 603  call xerrwv('lsoda--  istate .gt. 1 but lsoda not initialized  ',
     1   50, 3, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      go to 700
 604  call xerrwv('lsoda--  neq (=i1) .lt. 1     ',
     1   30, 4, 0, 1, neq(1), 0, 0, 0.0d0, 0.0d0)
      go to 700
 605  call xerrwv('lsoda--  istate = 3 and neq increased (i1 to i2)  ',
     1   50, 5, 0, 2, n, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 606  call xerrwv('lsoda--  itol (=i1) illegal   ',
     1   30, 6, 0, 1, itol, 0, 0, 0.0d0, 0.0d0)
      go to 700
 607  call xerrwv('lsoda--  iopt (=i1) illegal   ',
     1   30, 7, 0, 1, iopt, 0, 0, 0.0d0, 0.0d0)
      go to 700
 608  call xerrwv('lsoda--  jt (=i1) illegal     ',
     1   30, 8, 0, 1, jt, 0, 0, 0.0d0, 0.0d0)
      go to 700
 609  call xerrwv('lsoda--  ml (=i1) illegal.. .lt.0 or .ge.neq (=i2)',
     1   50, 9, 0, 2, ml, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 610  call xerrwv('lsoda--  mu (=i1) illegal.. .lt.0 or .ge.neq (=i2)',
     1   50, 10, 0, 2, mu, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 611  call xerrwv('lsoda--  ixpr (=i1) illegal   ',
     1   30, 11, 0, 1, ixpr, 0, 0, 0.0d0, 0.0d0)
      go to 700
 612  call xerrwv('lsoda--  mxstep (=i1) .lt. 0  ',
     1   30, 12, 0, 1, mxstep, 0, 0, 0.0d0, 0.0d0)
      go to 700
 613  call xerrwv('lsoda--  mxhnil (=i1) .lt. 0  ',
     1   30, 13, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
      go to 700
 614  call xerrwv('lsoda--  tout (=r1) behind t (=r2)      ',
     1   40, 14, 0, 0, 0, 0, 2, tout, t)
      call xerrwv('      integration direction is given by h0 (=r1)  ',
     1   50, 14, 0, 0, 0, 0, 1, h0, 0.0d0)
      go to 700
 615  call xerrwv('lsoda--  hmax (=r1) .lt. 0.0  ',
     1   30, 15, 0, 0, 0, 0, 1, hmax, 0.0d0)
      go to 700
 616  call xerrwv('lsoda--  hmin (=r1) .lt. 0.0  ',
     1   30, 16, 0, 0, 0, 0, 1, hmin, 0.0d0)
      go to 700
 617  call xerrwv(
     1  'lsoda--  rwork length needed, lenrw (=i1), exceeds lrw (=i2)',
     1   60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 618  call xerrwv(
     1  'lsoda--  iwork length needed, leniw (=i1), exceeds liw (=i2)',
     1   60, 18, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
      go to 700
 619  call xerrwv('lsoda--  rtol(i1) is r1 .lt. 0.0        ',
     1   40, 19, 0, 1, i, 0, 1, rtoli, 0.0d0)
      go to 700
 620  call xerrwv('lsoda--  atol(i1) is r1 .lt. 0.0        ',
     1   40, 20, 0, 1, i, 0, 1, atoli, 0.0d0)
      go to 700
 621  ewti = rwork(lewt+i-1)
      call xerrwv('lsoda--  ewt(i1) is r1 .le. 0.0         ',
     1   40, 21, 0, 1, i, 0, 1, ewti, 0.0d0)
      go to 700
 622  call xerrwv(
     1  'lsoda--  tout (=r1) too close to t(=r2) to start integration',
     1   60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  call xerrwv(
     1  'lsoda--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  ',
     1   60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624  call xerrwv(
     1  'lsoda--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ',
     1   60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  call xerrwv(
     1  'lsoda--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ',
     1   60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  call xerrwv('lsoda--  at start of problem, too much accuracy   ',
     1   50, 26, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      requested for precision of machine..  see tolsf (=r1) ',
     1   60, 26, 0, 0, 0, 0, 1, tolsf, 0.0d0)
      rwork(14) = tolsf
      go to 700
 627  call xerrwv('lsoda--  trouble from intdy. itask = i1, tout = r1',
     1   50, 27, 0, 1, itask, 0, 1, tout, 0.0d0)
      go to 700
 628  call xerrwv('lsoda--  mxordn (=i1) .lt. 0  ',
     1   30, 28, 0, 1, mxordn, 0, 0, 0.0d0, 0.0d0)
      go to 700
 629  call xerrwv('lsoda--  mxords (=i1) .lt. 0  ',
     1   30, 29, 0, 1, mxords, 0, 0, 0.0d0, 0.0d0)
c
 700  if (illin .eq. 5) go to 710
      illin = illin + 1
      istate = -3
      return
 710  call xerrwv('lsoda--  repeated occurrences of illegal input    ',
     1   50, 302, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
c
 800  call xerrwv('lsoda--  run aborted.. apparent infinite loop     ',
     1   50, 303, 2, 0, 0, 0, 0, 0.0d0, 0.0d0)
      return
c----------------------- end of subroutine lsoda -----------------------
      end
