      subroutine lsodes (f, neq, y, t, tout, itol, rtol, atol, itask,
     1            istate, iopt, rwork, lrw, iwork, liw, jac, mf)
      external f, jac
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
      double precision y, t, tout, rtol, atol, rwork
      dimension neq(1), y(1), rtol(1), atol(1), rwork(lrw), iwork(liw)
c-----------------------------------------------------------------------
c this is the march 30, 1987 version of
c lsodes.. livermore solver for ordinary differential equations
c          with general sparse jacobian matrices.
c this version is in double precision.
c
c lsodes solves the initial value problem for stiff or nonstiff
c systems of first order ode-s,
c     dy/dt = f(t,y) ,  or, in component form,
c     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
c lsodes is a variant of the lsode package, and is intended for
c problems in which the jacobian matrix df/dy has an arbitrary
c sparse structure (when the problem is stiff).
c
c authors..      alan c. hindmarsh,
c                computing and mathematics research division, l-316
c                lawrence livermore national laboratory
c                livermore, ca 94550.
c
c and            andrew h. sherman
c                j. s. nolen and associates
c                houston, tx 77084
c-----------------------------------------------------------------------
c references..
c 1.  alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c
c 2.  s. c. eisenstat, m. c. gursky, m. h. schultz, and a. h. sherman,
c     yale sparse matrix package.. i. the symmetric codes,
c     int. j. num. meth. eng., 18 (1982), pp. 1145-1151.
c
c 3.  s. c. eisenstat, m. c. gursky, m. h. schultz, and a. h. sherman,
c     yale sparse matrix package.. ii. the nonsymmetric codes,
c     research report no. 114, dept. of computer sciences, yale
c     university, 1977.
c-----------------------------------------------------------------------
c summary of usage.
c
c communication between the user and the lsodes package, for normal
c situations, is summarized here.  this summary describes only a subset
c of the full set of options available.  see the full description for
c details, including optional communication, nonstandard options,
c and instructions for special situations.  see also the example
c problem (with program and output) following this summary.
c
c a. first provide a subroutine of the form..
c               subroutine f (neq, t, y, ydot)
c               dimension y(neq), ydot(neq)
c which supplies the vector function f by loading ydot(i) with f(i).
c
c b. next determine (or guess) whether or not the problem is stiff.
c stiffness occurs when the jacobian matrix df/dy has an eigenvalue
c whose real part is negative and large in magnitude, compared to the
c reciprocal of the t span of interest.  if the problem is nonstiff,
c use a method flag mf = 10.  if it is stiff, there are two standard
c for the method flag, mf = 121 and mf = 222.  in both cases, lsodes
c requires the jacobian matrix in some form, and it treats this matrix
c in general sparse form, with sparsity structure determined internally.
c (for options where the user supplies the sparsity structure, see
c the full description of mf below.)
c
c c. if the problem is stiff, you are encouraged to supply the jacobian
c directly (mf = 121), but if this is not feasible, lsodes will
c compute it internally by difference quotients (mf = 222).
c if you are supplying the jacobian, provide a subroutine of the form..
c               subroutine jac (neq, t, y, j, ian, jan, pdj)
c               dimension y(1), ian(1), jan(1), pdj(1)
c here neq, t, y, and j are input arguments, and the jac routine is to
c load the array pdj (of length neq) with the j-th column of df/dy.
c i.e., load pdj(i) with df(i)/dy(j) for all relevant values of i.
c the arguments ian and jan should be ignored for normal situations.
c lsodes will call the jac routine with j = 1,2,...,neq.
c only nonzero elements need be loaded.  usually, a crude approximation
c to df/dy, possibly with fewer nonzero elements, will suffice.
c
c d. write a main program which calls subroutine lsodes once for
c each point at which answers are desired.  this should also provide
c for possible use of logical unit 6 for output of error messages
c by lsodes.  on the first call to lsodes, supply arguments as follows..
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
c istate = integer flag (input and output).  set istate = 1.
c iopt   = 0 to indicate no optional inputs used.
c rwork  = real work array of length at least..
c             20 + 16*neq            for mf = 10,
c             20 + (2 + 1./lenrat)*nnz + (11 + 9./lenrat)*neq
c                                    for mf = 121 or 222,
c          where..
c          nnz    = the number of nonzero elements in the sparse
c                   jacobian (if this is unknown, use an estimate), and
c          lenrat = the real to integer wordlength ratio (usually 1 in
c                   single precision and 2 in double precision).
c          in any case, the required size of rwork cannot generally
c          be predicted in advance if mf = 121 or 222, and the value
c          above is a rough estimate of a crude lower bound.  some
c          experimentation with this size may be necessary.
c          (when known, the correct required length is an optional
c          output, available in iwork(17).)
c lrw    = declared length of rwork (in user-s dimension).
c iwork  = integer work array of length at least 30.
c liw    = declared length of iwork (in user-s dimension).
c jac    = name of subroutine for jacobian matrix (mf = 121).
c          if used, this name must be declared external in calling
c          program.  if not used, pass a dummy name.
c mf     = method flag.  standard values are..
c          10  for nonstiff (adams) method, no jacobian used.
c          121 for stiff (bdf) method, user-supplied sparse jacobian.
c          222 for stiff method, internally generated sparse jacobian.
c note that the main program must declare arrays y, rwork, iwork,
c and possibly atol.
c
c e. the output from the first call (or any call) is..
c      y = array of computed values of y(t) vector.
c      t = corresponding value of independent variable (normally tout).
c istate = 2  if lsodes was successful, negative otherwise.
c          -1 means excess work done on this call (perhaps wrong mf).
c          -2 means excess accuracy requested (tolerances too small).
c          -3 means illegal input detected (see printed message).
c          -4 means repeated error test failures (check all inputs).
c          -5 means repeated convergence failures (perhaps bad jacobian
c             supplied or wrong choice of mf or tolerances).
c          -6 means error weight became zero during problem. (solution
c             component i vanished, and atol or atol(i) = 0.)
c          -7 means a fatal error return flag came from the sparse
c             solver cdrv by way of prjs or slss.  should never happen.
c          a return with istate = -1, -4, or -5 may result from using
c          an inappropriate sparsity structure, one that is quite
c          different from the initial structure.  consider calling
c          lsodes again with istate = 3 to force the structure to be
c          reevaluated.  see the full description of istate below.
c
c f. to continue the integration after a successful return, simply
c reset tout and call lsodes again.  no other parameters need be reset.
c
c-----------------------------------------------------------------------
c example problem.
c
c the following is a simple example problem, with the coding
c needed for its solution by lsodes.  the problem is from chemical
c kinetics, and consists of the following 12 rate equations..
c    dy1/dt  = -rk1*y1
c    dy2/dt  = rk1*y1 + rk11*rk14*y4 + rk19*rk14*y5
c                - rk3*y2*y3 - rk15*y2*y12 - rk2*y2
c    dy3/dt  = rk2*y2 - rk5*y3 - rk3*y2*y3 - rk7*y10*y3
c                + rk11*rk14*y4 + rk12*rk14*y6
c    dy4/dt  = rk3*y2*y3 - rk11*rk14*y4 - rk4*y4
c    dy5/dt  = rk15*y2*y12 - rk19*rk14*y5 - rk16*y5
c    dy6/dt  = rk7*y10*y3 - rk12*rk14*y6 - rk8*y6
c    dy7/dt  = rk17*y10*y12 - rk20*rk14*y7 - rk18*y7
c    dy8/dt  = rk9*y10 - rk13*rk14*y8 - rk10*y8
c    dy9/dt  = rk4*y4 + rk16*y5 + rk8*y6 + rk18*y7
c    dy10/dt = rk5*y3 + rk12*rk14*y6 + rk20*rk14*y7
c                + rk13*rk14*y8 - rk7*y10*y3 - rk17*y10*y12
c                - rk6*y10 - rk9*y10
c    dy11/dt = rk10*y8
c    dy12/dt = rk6*y10 + rk19*rk14*y5 + rk20*rk14*y7
c                - rk15*y2*y12 - rk17*y10*y12
c
c with rk1 = rk5 = 0.1,  rk4 = rk8 = rk16 = rk18 = 2.5,
c      rk10 = 5.0,  rk2 = rk6 = 10.0,  rk14 = 30.0,
c      rk3 = rk7 = rk9 = rk11 = rk12 = rk13 = rk19 = rk20 = 50.0,
c      rk15 = rk17 = 100.0.
c
c the t interval is from 0 to 1000, and the initial conditions
c are y1 = 1, y2 = y3 = ... = y12 = 0.  the problem is stiff.
c
c the following coding solves this problem with lsodes, using mf = 121
c and printing results at t = .1, 1., 10., 100., 1000.  it uses
c itol = 1 and mixed relative/absolute tolerance controls.
c during the run and at the end, statistical quantities of interest
c are printed (see optional outputs in the full description below).
c
c     external fex, jex
c     double precision atol, rtol, rwork, t, tout, y
c     dimension y(12), rwork(500), iwork(30)
c     data lrw/500/, liw/30/
c     neq = 12
c     do 10 i = 1,neq
c 10    y(i) = 0.0d0
c     y(1) = 1.0d0
c     t = 0.0d0
c     tout = 0.1d0
c     itol = 1
c     rtol = 1.0d-4
c     atol = 1.0d-6
c     itask = 1
c     istate = 1
c     iopt = 0
c     mf = 121
c     do 40 iout = 1,5
c       call lsodes (fex, neq, y, t, tout, itol, rtol, atol,
c    1     itask, istate, iopt, rwork, lrw, iwork, liw, jex, mf)
c       write(6,30)t,iwork(11),rwork(11),(y(i),i=1,neq)
c 30    format(//7h at t =,e11.3,4x,
c    1    12h no. steps =,i5,4x,12h last step =,e11.3/
c    2    13h  y array =  ,4e14.5/13x,4e14.5/13x,4e14.5)
c       if (istate .lt. 0) go to 80
c       tout = tout*10.0d0
c 40    continue
c     lenrw = iwork(17)
c     leniw = iwork(18)
c     nst = iwork(11)
c     nfe = iwork(12)
c     nje = iwork(13)
c     nlu = iwork(21)
c     nnz = iwork(19)
c     nnzlu = iwork(25) + iwork(26) + neq
c     write (6,70) lenrw,leniw,nst,nfe,nje,nlu,nnz,nnzlu
c 70  format(//22h required rwork size =,i4,15h   iwork size =,i4/
c    1   12h no. steps =,i4,12h   no. f-s =,i4,12h   no. j-s =,i4,
c    2   13h   no. lu-s =,i4/23h no. of nonzeros in j =,i5,
c    3   26h   no. of nonzeros in lu =,i5)
c     stop
c 80  write(6,90)istate
c 90  format(///22h error halt.. istate =,i3)
c     stop
c     end
c
c     subroutine fex (neq, t, y, ydot)
c     double precision t, y, ydot
c     double precision rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9,
c    1   rk10, rk11, rk12, rk13, rk14, rk15, rk16, rk17
c     dimension y(12), ydot(12)
c     data rk1/0.1d0/, rk2/10.0d0/, rk3/50.0d0/, rk4/2.5d0/, rk5/0.1d0/,
c    1   rk6/10.0d0/, rk7/50.0d0/, rk8/2.5d0/, rk9/50.0d0/, rk10/5.0d0/,
c    2   rk11/50.0d0/, rk12/50.0d0/, rk13/50.0d0/, rk14/30.0d0/,
c    3   rk15/100.0d0/, rk16/2.5d0/, rk17/100.0d0/, rk18/2.5d0/,
c    4   rk19/50.0d0/, rk20/50.0d0/
c     ydot(1)  = -rk1*y(1)
c     ydot(2)  = rk1*y(1) + rk11*rk14*y(4) + rk19*rk14*y(5)
c    1           - rk3*y(2)*y(3) - rk15*y(2)*y(12) - rk2*y(2)
c     ydot(3)  = rk2*y(2) - rk5*y(3) - rk3*y(2)*y(3) - rk7*y(10)*y(3)
c    1           + rk11*rk14*y(4) + rk12*rk14*y(6)
c     ydot(4)  = rk3*y(2)*y(3) - rk11*rk14*y(4) - rk4*y(4)
c     ydot(5)  = rk15*y(2)*y(12) - rk19*rk14*y(5) - rk16*y(5)
c     ydot(6)  = rk7*y(10)*y(3) - rk12*rk14*y(6) - rk8*y(6)
c     ydot(7)  = rk17*y(10)*y(12) - rk20*rk14*y(7) - rk18*y(7)
c     ydot(8)  = rk9*y(10) - rk13*rk14*y(8) - rk10*y(8)
c     ydot(9)  = rk4*y(4) + rk16*y(5) + rk8*y(6) + rk18*y(7)
c     ydot(10) = rk5*y(3) + rk12*rk14*y(6) + rk20*rk14*y(7)
c    1           + rk13*rk14*y(8) - rk7*y(10)*y(3) - rk17*y(10)*y(12)
c    2           - rk6*y(10) - rk9*y(10)
c     ydot(11) = rk10*y(8)
c     ydot(12) = rk6*y(10) + rk19*rk14*y(5) + rk20*rk14*y(7)
c    1           - rk15*y(2)*y(12) - rk17*y(10)*y(12)
c     return
c     end
c
c     subroutine jex (neq, t, y, j, ia, ja, pdj)
c     double precision t, y, pdj
c     double precision rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9,
c    1   rk10, rk11, rk12, rk13, rk14, rk15, rk16, rk17
c     dimension y(1), ia(1), ja(1), pdj(1)
c     data rk1/0.1d0/, rk2/10.0d0/, rk3/50.0d0/, rk4/2.5d0/, rk5/0.1d0/,
c    1   rk6/10.0d0/, rk7/50.0d0/, rk8/2.5d0/, rk9/50.0d0/, rk10/5.0d0/,
c    2   rk11/50.0d0/, rk12/50.0d0/, rk13/50.0d0/, rk14/30.0d0/,
c    3   rk15/100.0d0/, rk16/2.5d0/, rk17/100.0d0/, rk18/2.5d0/,
c    4   rk19/50.0d0/, rk20/50.0d0/
c     go to (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), j
c 1   pdj(1) = -rk1
c     pdj(2) = rk1
c     return
c 2   pdj(2) = -rk3*y(3) - rk15*y(12) - rk2
c     pdj(3) = rk2 - rk3*y(3)
c     pdj(4) = rk3*y(3)
c     pdj(5) = rk15*y(12)
c     pdj(12) = -rk15*y(12)
c     return
c 3   pdj(2) = -rk3*y(2)
c     pdj(3) = -rk5 - rk3*y(2) - rk7*y(10)
c     pdj(4) = rk3*y(2)
c     pdj(6) = rk7*y(10)
c     pdj(10) = rk5 - rk7*y(10)
c     return
c 4   pdj(2) = rk11*rk14
c     pdj(3) = rk11*rk14
c     pdj(4) = -rk11*rk14 - rk4
c     pdj(9) = rk4
c     return
c 5   pdj(2) = rk19*rk14
c     pdj(5) = -rk19*rk14 - rk16
c     pdj(9) = rk16
c     pdj(12) = rk19*rk14
c     return
c 6   pdj(3) = rk12*rk14
c     pdj(6) = -rk12*rk14 - rk8
c     pdj(9) = rk8
c     pdj(10) = rk12*rk14
c     return
c 7   pdj(7) = -rk20*rk14 - rk18
c     pdj(9) = rk18
c     pdj(10) = rk20*rk14
c     pdj(12) = rk20*rk14
c     return
c 8   pdj(8) = -rk13*rk14 - rk10
c     pdj(10) = rk13*rk14
c     pdj(11) = rk10
c 9   return
c 10  pdj(3) = -rk7*y(3)
c     pdj(6) = rk7*y(3)
c     pdj(7) = rk17*y(12)
c     pdj(8) = rk9
c     pdj(10) = -rk7*y(3) - rk17*y(12) - rk6 - rk9
c     pdj(12) = rk6 - rk17*y(12)
c 11  return
c 12  pdj(2) = -rk15*y(2)
c     pdj(5) = rk15*y(2)
c     pdj(7) = rk17*y(10)
c     pdj(10) = -rk17*y(10)
c     pdj(12) = -rk15*y(2) - rk17*y(10)
c     return
c     end
c
c the output of this program (on a cray-1 in single precision)
c is as follows..
c
c
c at t =  1.000e-01     no. steps =   12     last step =  1.515e-02
c  y array =     9.90050e-01   6.28228e-03   3.65313e-03   7.51934e-07
c                1.12167e-09   1.18458e-09   1.77291e-12   3.26476e-07
c                5.46720e-08   9.99500e-06   4.48483e-08   2.76398e-06
c
c
c at t =  1.000e+00     no. steps =   33     last step =  7.880e-02
c  y array =     9.04837e-01   9.13105e-03   8.20622e-02   2.49177e-05
c                1.85055e-06   1.96797e-06   1.46157e-07   2.39557e-05
c                3.26306e-05   7.21621e-04   5.06433e-05   3.05010e-03
c
c
c at t =  1.000e+01     no. steps =   48     last step =  1.239e+00
c  y array =     3.67876e-01   3.68958e-03   3.65133e-01   4.48325e-05
c                6.10798e-05   4.33148e-05   5.90211e-05   1.18449e-04
c                3.15235e-03   3.56531e-03   4.15520e-03   2.48741e-01
c
c
c at t =  1.000e+02     no. steps =   91     last step =  3.764e+00
c  y array =     4.44981e-05   4.42666e-07   4.47273e-04  -3.53257e-11
c                2.81577e-08  -9.67741e-11   2.77615e-07   1.45322e-07
c                1.56230e-02   4.37394e-06   1.60104e-02   9.52246e-01
c
c
c at t =  1.000e+03     no. steps =  111     last step =  4.156e+02
c  y array =    -2.65492e-13   2.60539e-14  -8.59563e-12   6.29355e-14
c               -1.78066e-13   5.71471e-13  -1.47561e-12   4.58078e-15
c                1.56314e-02   1.37878e-13   1.60184e-02   9.52719e-01
c
c
c required rwork size = 442   iwork size =  30
c no. steps = 111   no. f-s = 142   no. j-s =   2   no. lu-s =  20
c no. of nonzeros in j =   44   no. of nonzeros in lu =   50
c-----------------------------------------------------------------------
c full description of user interface to lsodes.
c
c the user interface to lsodes consists of the following parts.
c
c i.   the call sequence to subroutine lsodes, which is a driver
c      routine for the solver.  this includes descriptions of both
c      the call sequence arguments and of user-supplied routines.
c      following these descriptions is a description of
c      optional inputs available through the call sequence, and then
c      a description of optional outputs (in the work arrays).
c
c ii.  descriptions of other routines in the lsodes package that may be
c      (optionally) called by the user.  these provide the ability to
c      alter error message handling, save and restore the internal
c      common, and obtain specified derivatives of the solution y(t).
c
c iii. descriptions of common blocks to be declared in overlay
c      or similar environments, or to be saved when doing an interrupt
c      of the problem and continued solution later.
c
c iv.  description of two routines in the lsodes package, either of
c      which the user may replace with his own version, if desired.
c      these relate to the measurement of errors.
c
c-----------------------------------------------------------------------
c part i.  call sequence.
c
c the call sequence parameters used for input only are
c     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, mf,
c and those used for both input and output are
c     y, t, istate.
c the work arrays rwork and iwork are also used for conditional and
c optional inputs and optional outputs.  (the term output here refers
c to the return from subroutine lsodes to the user-s calling program.)
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
c          externally to lsodes, an extra call to f should be made
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
c          (the lsodes package accesses only neq(1).)  in either case,
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
c          pass it to f and/or jac.  (the lsodes package accesses only
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
c              4     array      array      rtol(i)*abs(y(i)) + atol(i)
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
c             neq, itol, rtol, atol, iopt, lrw, liw, mf,
c             the conditional inputs ia and ja,
c             and any of the optional inputs except h0.
c             in particular, if miter = 1 or 2, a call with istate = 3
c             will cause the sparsity structure of the problem to be
c             recomputed (or reread from ia and ja if moss = 0).
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
c          -7  means a fatal error return flag came from the sparse
c              solver cdrv by way of prjs or slss (numerical
c              factorization or backsolve).  this should never happen.
c              the integration was successful as far as t.
c
c          note.. an error return with istate = -1, -4, or -5 and with
c          miter = 1 or 2 may mean that the sparsity structure of the
c          problem has changed significantly since it was last
c          determined (or input).  in that case, one can attempt to
c          complete the integration by setting istate = 3 on the next
c          call, so that a new structure determination is done.
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
c rwork  = a work array used for a mixture of real (double precision)
c          and integer work space.
c          the length of rwork (in real words) must be at least
c             20 + nyh*(maxord + 1) + 3*neq + lwm    where
c          nyh    = the initial value of neq,
c          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
c                   smaller value is given as an optional input),
c          lwm = 0                                    if miter = 0,
c          lwm = 2*nnz + 2*neq + (nnz+9*neq)/lenrat   if miter = 1,
c          lwm = 2*nnz + 2*neq + (nnz+10*neq)/lenrat  if miter = 2,
c          lwm = neq + 2                              if miter = 3.
c          in the above formulas,
c          nnz    = number of nonzero elements in the jacobian matrix.
c          lenrat = the real to integer wordlength ratio (usually 1 in
c                   single precision and 2 in double precision).
c          (see the mf description for meth and miter.)
c          thus if maxord has its default value and neq is constant,
c          the minimum length of rwork is..
c             20 + 16*neq        for mf = 10,
c             20 + 16*neq + lwm  for mf = 11, 111, 211, 12, 112, 212,
c             22 + 17*neq        for mf = 13,
c             20 +  9*neq        for mf = 20,
c             20 +  9*neq + lwm  for mf = 21, 121, 221, 22, 122, 222,
c             22 + 10*neq        for mf = 23.
c          if miter = 1 or 2, the above formula for lwm is only a
c          crude lower bound.  the required length of rwork cannot
c          be readily predicted in general, as it depends on the
c          sparsity structure of the problem.  some experimentation
c          may be necessary.
c
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
c             31 + neq + nnz   if moss = 0 and miter = 1 or 2, or
c             30               otherwise.
c          (nnz is the number of nonzero elements in df/dy.)
c
c          in lsodes, iwork is used only for conditional and
c          optional inputs and optional outputs.
c
c          the following two blocks of words in iwork are conditional
c          inputs, required if moss = 0 and miter = 1 or 2, but not
c          otherwise (see the description of mf for moss).
c            iwork(30+j) = ia(j)     (j=1,...,neq+1)
c            iwork(31+neq+k) = ja(k) (k=1,...,nnz)
c          the two arrays ia and ja describe the sparsity structure
c          to be assumed for the jacobian matrix.  ja contains the row
c          indices where nonzero elements occur, reading in columnwise
c          order, and ia contains the starting locations in ja of the
c          descriptions of columns 1,...,neq, in that order, with
c          ia(1) = 1.  thus, for each column index j = 1,...,neq, the
c          values of the row index i in column j where a nonzero
c          element may occur are given by
c            i = ja(k),  where   ia(j) .le. k .lt. ia(j+1).
c          if nnz is the total number of nonzero locations assumed,
c          then the length of the ja array is nnz, and ia(neq+1) must
c          be nnz + 1.  duplicate entries are not allowed.
c
c liw    = the length of the array iwork, as declared by the user.
c          (this will be checked by the solver.)
c
c note..  the work arrays must not be altered between calls to lsodes
c for the same problem, except possibly for the conditional and
c optional inputs, and except for the last 3*neq words of rwork.
c the latter space is used for internal scratch space, and so is
c available for use by the user outside lsodes between calls, if
c desired (but not for use by f or jac).
c
c jac    = name of user-supplied routine (miter = 1 or moss = 1) to
c          compute the jacobian matrix, df/dy, as a function of
c          the scalar t and the vector y.  it is to have the form
c               subroutine jac (neq, t, y, j, ian, jan, pdj)
c               dimension y(1), ian(1), jan(1), pdj(1)
c          where neq, t, y, j, ian, and jan are input, and the array
c          pdj, of length neq, is to be loaded with column j
c          of the jacobian on output.  thus df(i)/dy(j) is to be
c          loaded into pdj(i) for all relevant values of i.
c          here t and y have the same meaning as in subroutine f,
c          and j is a column index (1 to neq).  ian and jan are
c          undefined in calls to jac for structure determination
c          (moss = 1).  otherwise, ian and jan are structure
c          descriptors, as defined under optional outputs below, and
c          so can be used to determine the relevant row indices i, if
c          desired.  (in the dimension statement above, 1 is a
c          dummy dimension.. it can be replaced by any value.)
c               jac need not provide df/dy exactly.  a crude
c          approximation (possibly with greater sparsity) will do.
c               in any case, pdj is preset to zero by the solver,
c          so that only the nonzero elements need be loaded by jac.
c          calls to jac are made with j = 1,...,neq, in that order, and
c          each such set of calls is preceded by a call to f with the
c          same arguments neq, t, and y.  thus to gain some efficiency,
c          intermediate quantities shared by both calculations may be
c          saved in a user common block by f and not recomputed by jac,
c          if desired.  jac must not alter its input arguments.
c          jac must be declared external in the calling program.
c               subroutine jac may access user-defined quantities in
c          neq(2),... and y(neq(1)+1),... if neq is an array
c          (dimensioned in jac) and y has length exceeding neq(1).
c          see the descriptions of neq and y above.
c
c mf     = the method flag.  used only for input.
c          mf has three decimal digits-- moss, meth, miter--
c             mf = 100*moss + 10*meth + miter.
c          moss indicates the method to be used to obtain the sparsity
c          structure of the jacobian matrix if miter = 1 or 2..
c            moss = 0 means the user has supplied ia and ja
c                     (see descriptions under iwork above).
c            moss = 1 means the user has supplied jac (see below)
c                     and the structure will be obtained from neq
c                     initial calls to jac.
c            moss = 2 means the structure will be obtained from neq+1
c                     initial calls to f.
c          meth indicates the basic linear multistep method..
c            meth = 1 means the implicit adams method.
c            meth = 2 means the method based on backward
c                     differentiation formulas (bdf-s).
c          miter indicates the corrector iteration method..
c            miter = 0 means functional iteration (no jacobian matrix
c                      is involved).
c            miter = 1 means chord iteration with a user-supplied
c                      sparse jacobian, given by subroutine jac.
c            miter = 2 means chord iteration with an internally
c                      generated (difference quotient) sparse jacobian
c                      (using ngp extra calls to f per df/dy value,
c                      where ngp is an optional output described below.)
c            miter = 3 means chord iteration with an internally
c                      generated diagonal jacobian approximation.
c                      (using 1 extra call to f per df/dy evaluation).
c          if miter = 1 or moss = 1, the user must supply a subroutine
c          jac (the name is arbitrary) as described above under jac.
c          otherwise, a dummy argument can be used.
c
c          the standard choices for mf are..
c            mf = 10  for a nonstiff problem,
c            mf = 21 or 22 for a stiff problem with ia/ja supplied
c                     (21 if jac is supplied, 22 if not),
c            mf = 121 for a stiff problem with jac supplied,
c                     but not ia/ja,
c            mf = 222 for a stiff problem with neither ia/ja nor
c                     jac supplied.
c          the sparseness structure can be changed during the
c          problem by making a call to lsodes with istate = 3.
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
c seth    rwork(8)  the element threshhold for sparsity determination
c                   when moss = 1 or 2.  if the absolute value of
c                   an estimated jacobian element is .le. seth, it
c                   will be assumed to be absent in the structure.
c                   the default value of seth is 0.
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
c as optional additional output from lsodes, the variables listed
c below are quantities related to the performance of lsodes
c which are available to the user.  these are communicated by way of
c the work arrays, but also have internal mnemonic names as shown.
c except where stated otherwise, all of these outputs are defined
c on any successful return from lsodes, and on any return with
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
c nst     iwork(11) the number of steps taken for the problem so far.
c
c nfe     iwork(12) the number of f evaluations for the problem so far,
c                   excluding those for structure determination
c                   (moss = 2).
c
c nje     iwork(13) the number of jacobian evaluations for the problem
c                   so far, excluding those for structure determination
c                   (moss = 1).
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
c nnz     iwork(19) the number of nonzero elements in the jacobian
c                   matrix, including the diagonal (miter = 1 or 2).
c                   (this may differ from that given by ia(neq+1)-1
c                   if moss = 0, because of added diagonal entries.)
c
c ngp     iwork(20) the number of groups of column indices, used in
c                   difference quotient jacobian aproximations if
c                   miter = 2.  this is also the number of extra f
c                   evaluations needed for each jacobian evaluation.
c
c nlu     iwork(21) the number of sparse lu decompositions for the
c                   problem so far.
c
c lyh     iwork(22) the base address in rwork of the history array yh,
c                   described below in this list.
c
c ipian   iwork(23) the base address of the structure descriptor array
c                   ian, described below in this list.
c
c ipjan   iwork(24) the base address of the structure descriptor array
c                   jan, described below in this list.
c
c nzl     iwork(25) the number of nonzero elements in the strict lower
c                   triangle of the lu factorization used in the chord
c                   iteration (miter = 1 or 2).
c
c nzu     iwork(26) the number of nonzero elements in the strict upper
c                   triangle of the lu factorization used in the chord
c                   iteration (miter = 1 or 2).
c                   the total number of nonzeros in the factorization
c                   is therefore nzl + nzu + neq.
c
c the following four arrays are segments of the rwork array which
c may also be of interest to the user as optional outputs.
c for each array, the table below gives its internal name,
c its base address, and its description.
c for yh and acor, the base addresses are in rwork (a real array).
c the integer arrays ian and jan are to be obtained by declaring an
c integer array iwk and identifying iwk(1) with rwork(21), using either
c an equivalence statement or a subroutine call.  then the base
c addresses ipian (of ian) and ipjan (of jan) in iwk are to be obtained
c as optional outputs iwork(23) and iwork(24), respectively.
c thus ian(1) is iwk(ipian), etc.
c
c name    base address      description
c
c ian    ipian (in iwk)  structure descriptor array of size neq + 1.
c jan    ipjan (in iwk)  structure descriptor array of size nnz.
c         (see above)    ian and jan together describe the sparsity
c                        structure of the jacobian matrix, as used by
c                        lsodes when miter = 1 or 2.
c                        jan contains the row indices of the nonzero
c                        locations, reading in columnwise order, and
c                        ian contains the starting locations in jan of
c                        the descriptions of columns 1,...,neq, in
c                        that order, with ian(1) = 1.  thus for each
c                        j = 1,...,neq, the row indices i of the
c                        nonzero locations in column j are
c                        i = jan(k),  ian(j) .le. k .lt. ian(j+1).
c                        note that ian(neq+1) = nnz + 1.
c                        (if moss = 0, ian/jan may differ from the
c                        input ia/ja because of a different ordering
c                        in each column, and added diagonal entries.)
c
c yh      lyh            the nordsieck history array, of size nyh by
c          (optional     (nqcur + 1), where nyh is the initial value
c          output)       of neq.  for j = 0,1,...,nqcur, column j+1
c                        of yh contains hcur**j/factorial(j) times
c                        the j-th derivative of the interpolating
c                        polynomial currently representing the solution,
c                        evaluated at t = tcur.  the base address lyh
c                        is another optional output, listed above.
c
c acor     lenrw-neq+1   array of size neq used for the accumulated
c                        corrections on each step, scaled on output
c                        to represent the estimated local error in y
c                        on the last step.  this is the vector e in
c                        the description of the error control.  it is
c                        defined only on a successful return from
c                        lsodes.
c
c-----------------------------------------------------------------------
c part ii.  other routines callable.
c
c the following are optional calls which the user may make to
c gain additional capabilities in conjunction with lsodes.
c (the routines xsetun and xsetf are designed to conform to the
c slatec error handling package.)
c
c     form of call                  function
c   call xsetun(lun)          set the logical unit number, lun, for
c                             output of messages from lsodes, if
c                             the default is not desired.
c                             the default value of lun is 6.
c
c   call xsetf(mflag)         set a flag to control the printing of
c                             messages by lsodes.
c                             mflag = 0 means do not print. (danger..
c                             this risks losing valuable information.)
c                             mflag = 1 means print (the default).
c
c                             either of the above calls may be made at
c                             any time and will take effect immediately.
c
c   call srcms(rsav,isav,job) saves and restores the contents of
c                             the internal common blocks used by
c                             lsodes (see part iii below).
c                             rsav must be a real array of length 224
c                             or more, and isav must be an integer
c                             array of length 75 or more.
c                             job=1 means save common into rsav/isav.
c                             job=2 means restore common from rsav/isav.
c                                srcms is useful if one is
c                             interrupting a run and restarting
c                             later, or alternating between two or
c                             more problems solved with lsodes.
c
c   call intdy(,,,,,)         provide derivatives of y, of various
c        (see below)          orders, at a specified point t, if
c                             desired.  it may be called only after
c                             a successful return from lsodes.
c
c the detailed instructions for using intdy are as follows.
c the form of the call is..
c
c   lyh = iwork(22)
c   call intdy (t, k, rwork(lyh), nyh, dky, iflag)
c
c the input parameters are..
c
c t         = value of independent variable where answers are desired
c             (normally the same as the t last returned by lsodes).
c             for valid results, t must lie between tcur - hu and tcur.
c             (see optional outputs for tcur and hu.)
c k         = integer order of the derivative desired.  k must satisfy
c             0 .le. k .le. nqcur, where nqcur is the current order
c             (see optional outputs).  the capability corresponding
c             to k = 0, i.e. computing y(t), is already provided
c             by lsodes directly.  since nqcur .ge. 1, the first
c             derivative dy/dt is always available with intdy.
c lyh       = the base address of the history array yh, obtained
c             as an optional output as shown above.
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
c if lsodes is to be used in an overlay situation, the user
c must declare, in the primary overlay, the variables in..
c   (1) the call sequence to lsodes,
c   (2) the three internal common blocks
c         /ls0001/  of length  257  (218 double precision words
c                         followed by 39 integer words),
c         /lss001/  of length  40    ( 6 double precision words
c                         followed by 34 integer words),
c         /eh0001/  of length  2 (integer words).
c
c if lsodes is used on a system in which the contents of internal
c common blocks are not preserved between calls, the user should
c declare the above three common blocks in his main program to insure
c that their contents are preserved.
c
c if the solution of a given problem by lsodes is to be interrupted
c and then later continued, such as when restarting an interrupted run
c or alternating between two or more problems, the user should save,
c following the return from the last lsodes call prior to the
c interruption, the contents of the call sequence variables and the
c internal common blocks, and later restore these values before the
c next lsodes call for that problem.  to save and restore the common
c blocks, use subroutine srcms (see part ii above).
c
c-----------------------------------------------------------------------
c part iv.  optionally replaceable solver routines.
c
c below are descriptions of two routines in the lsodes package which
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
c where neq, itol, rtol, and atol are as in the lsodes call sequence,
c ycur contains the current dependent variable vector, and
c ewt is the array of weights set by ewset.
c
c if the user supplies this subroutine, it must return in ewt(i)
c (i = 1,...,neq) a positive quantity suitable for comparing errors
c in y(i) to.  the ewt array returned by ewset is passed to the
c vnorm routine (see below), and also used by lsodes in the computation
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
c value of vnorm suitable for use in the error control in lsodes.
c none of the arguments should be altered by vnorm.
c for example, a user-supplied vnorm routine might..
c   -substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
c   -ignore some components of v in the norm, with the effect of
c    suppressing the error control on those components of y.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c other routines in the lsodes package.
c
c in addition to subroutine lsodes, the lsodes package includes the
c following subroutines and function routines..
c  iprep    acts as an iterface between lsodes and prep, and also does
c           adjusting of work space pointers and work arrays.
c  prep     is called by iprep to compute sparsity and do sparse matrix
c           preprocessing if miter = 1 or 2.
c  jgroup   is called by prep to compute groups of jacobian column
c           indices for use when miter = 2.
c  adjlr    adjusts the length of required sparse matrix work space.
c           it is called by prep.
c  cntnzu   is called by prep and counts the nonzero elements in the
c           strict upper triangle of j + j-transpose, where j = df/dy.
c  intdy    computes an interpolated value of the y vector at t = tout.
c  stode    is the core integrator, which does one step of the
c           integration and the associated error control.
c  cfode    sets all method coefficients and test constants.
c  prjs     computes and preprocesses the jacobian matrix j = df/dy
c           and the newton iteration matrix p = i - h*l0*j.
c  slss     manages solution of linear system in chord iteration.
c  ewset    sets the error weight vector ewt before each step.
c  vnorm    computes the weighted r.m.s. norm of a vector.
c  srcms    is a user-callable routine to save and restore
c           the contents of the internal common blocks.
c  odrv     constructs a reordering of the rows and columns of
c           a matrix by the minimum degree algorithm.  odrv is a
c           driver routine which calls subroutines md, mdi, mdm,
c           mdp, mdu, and sro.  see ref. 2 for details.  (the odrv
c           module has been modified since ref. 2, however.)
c  cdrv     performs reordering, symbolic factorization, numerical
c           factorization, or linear system solution operations,
c           depending on a path argument ipath.  cdrv is a
c           driver routine which calls subroutines nroc, nsfc,
c           nnfc, nnsc, and nntc.  see ref. 3 for details.
c           lsodes uses cdrv to solve linear systems in which the
c           coefficient matrix is  p = i - con*j, where i is the
c           identity, con is a scalar, and j is an approximation to
c           the jacobian df/dy.  because cdrv deals with rowwise
c           sparsity descriptions, cdrv works with p-transpose, not p.
c  d1mach   computes the unit roundoff in a machine-independent manner.
c  xerrwv, xsetun, and xsetf   handle the printing of all error
c           messages and warnings.  xerrwv is machine-dependent.
c note..  vnorm and d1mach are function routines.
c all the others are subroutines.
c
c the intrinsic and external routines used by lsodes are..
c dabs, dmax1, dmin1, dfloat, max0, min0, mod, dsign, dsqrt, and write.
c
c a block data subprogram is also included with the package,
c for loading some of the variables in internal common.
c
c-----------------------------------------------------------------------
c the following card is for optimized compilation on lll compilers.
clll. optimize
c-----------------------------------------------------------------------
      external prjs, slss
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, i1, i2, iflag, imax, imul, imxer, ipflag, ipgo, irem,
     1   j, kgo, lenrat, lenyht, leniw, lenrw, lf0, lia, lja,
     2   lrtem, lwtem, lyhd, lyhn, mf1, mord, mxhnl0, mxstp0, ncolm
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision con0, conmin, ccmxj, psmall, rbig, seth
      double precision atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0,
     2   d1mach, vnorm
      dimension mord(2)
      logical ihit
c-----------------------------------------------------------------------
c the following two internal common blocks contain
c (a) variables which are local to any subroutine but whose values must
c     be preserved between calls to the routine (own variables), and
c (b) variables which are communicated between subroutines.
c the structure of each block is as follows..  all real variables are
c listed first, followed by all integers.  within each type, the
c variables are grouped with those local to subroutine lsodes first,
c then those local to subroutine stode or subroutine prjs
c (no other routines have own variables), and finally those used
c for communication.  the block ls0001 is declared in subroutines
c lsodes, iprep, prep, intdy, stode, prjs, and slss.  the block lss001
c is declared in subroutines lsodes, iprep, prep, prjs, and slss.
c groups of variables are replaced by dummy arrays in the common
c declarations in routines where those variables are not used.
c-----------------------------------------------------------------------
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c
      common /lss001/ con0, conmin, ccmxj, psmall, rbig, seth,
     1   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     4   nslj, ngp, nlu, nnz, nsp, nzl, nzu
c
      data mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
c-----------------------------------------------------------------------
c in the data statement below, set lenrat equal to the ratio of
c the wordlength for a real number to that for an integer.  usually,
c lenrat = 1 for single precision and 2 for double precision.  if the
c true ratio is not an integer, use the next smaller integer (.ge. 1).
c-----------------------------------------------------------------------
      data lenrat/2/
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
c if istate = 1, the final setting of work space pointers, the matrix
c preprocessing, and other initializations are done in block c.
c
c first check legality of the non-optional inputs neq, itol, iopt,
c mf, ml, and mu.
c-----------------------------------------------------------------------
      if (neq(1) .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      moss = mf/100
      mf1 = mf - 100*moss
      meth = mf1/10
      miter = mf1 - 10*meth
      if (moss .lt. 0 .or. moss .gt. 2) go to 608
      if (meth .lt. 1 .or. meth .gt. 2) go to 608
      if (miter .lt. 0 .or. miter .gt. 3) go to 608
      if (miter .eq. 0 .or. miter .eq. 3) moss = 0
c next process and check the optional inputs. --------------------------
      if (iopt .eq. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      if (istate .eq. 1) h0 = 0.0d0
      hmxi = 0.0d0
      hmin = 0.0d0
      seth = 0.0d0
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
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      if ((tout - t)*h0 .lt. 0.0d0) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. 0.0d0) go to 615
      hmxi = 0.0d0
      if (hmax .gt. 0.0d0) hmxi = 1.0d0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0d0) go to 616
      seth = rwork(8)
      if (seth .lt. 0.0d0) go to 609
c check rtol and atol for legality. ------------------------------------
 60   rtoli = rtol(1)
      atoli = atol(1)
      do 65 i = 1,n
        if (itol .ge. 3) rtoli = rtol(i)
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        if (rtoli .lt. 0.0d0) go to 619
        if (atoli .lt. 0.0d0) go to 620
 65     continue
c-----------------------------------------------------------------------
c compute required work array lengths, as far as possible, and test
c these against lrw and liw.  then set tentative pointers for work
c arrays.  pointers to rwork/iwork segments are named by prefixing l to
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).
c segments of rwork (in order) are denoted  wm, yh, savf, ewt, acor.
c if miter = 1 or 2, the required length of the matrix work space wm
c is not yet known, and so a crude minimum value is used for the
c initial tests of lrw and liw, and yh is temporarily stored as far
c to the right in rwork as possible, to leave the maximum amount
c of space for wm for matrix preprocessing.  thus if miter = 1 or 2
c and moss .ne. 2, some of the segments of rwork are temporarily
c omitted, as they are not needed in the preprocessing.  these
c omitted segments are.. acor if istate = 1, ewt and acor if istate = 3
c and moss = 1, and savf, ewt, and acor if istate = 3 and moss = 0.
c-----------------------------------------------------------------------
      lrat = lenrat
      if (istate .eq. 1) nyh = n
      lwmin = 0
      if (miter .eq. 1) lwmin = 4*n + 10*n/lrat
      if (miter .eq. 2) lwmin = 4*n + 11*n/lrat
      if (miter .eq. 3) lwmin = n + 2
      lenyh = (maxord+1)*nyh
      lrest = lenyh + 3*n
      lenrw = 20 + lwmin + lrest
      iwork(17) = lenrw
      leniw = 30
      if (moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3)
     1   leniw = leniw + n + 1
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618
      lia = 31
      if (moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3)
     1   leniw = leniw + iwork(lia+n) - 1
      iwork(18) = leniw
      if (leniw .gt. liw) go to 618
      lja = lia + n + 1
      lia = min0(lia,liw)
      lja = min0(lja,liw)
      lwm = 21
      if (istate .eq. 1) nq = 1
      ncolm = min0(nq+1,maxord+2)
      lenyhm = ncolm*nyh
      lenyht = lenyh
      if (miter .eq. 1 .or. miter .eq. 2) lenyht = lenyhm
      imul = 2
      if (istate .eq. 3) imul = moss
      if (moss .eq. 2) imul = 3
      lrtem = lenyht + imul*n
      lwtem = lwmin
      if (miter .eq. 1 .or. miter .eq. 2) lwtem = lrw - 20 - lrtem
      lenwk = lwtem
      lyhn = lwm + lwtem
      lsavf = lyhn + lenyht
      lewt = lsavf + n
      lacor = lewt + n
      istatc = istate
      if (istate .eq. 1) go to 100
c-----------------------------------------------------------------------
c istate = 3.  move yh to its new location.
c note that only the part of yh needed for the next step, namely
c min(nq+1,maxord+2) columns, is actually moved.
c a temporary error weight array ewt is loaded if moss = 2.
c sparse matrix processing is done in iprep/prep if miter = 1 or 2.
c if maxord was reduced below nq, then the pointers are finally set
c so that savf is identical to yh(*,maxord+2).
c-----------------------------------------------------------------------
      lyhd = lyh - lyhn
      imax = lyhn - 1 + lenyhm
c move yh.  branch for move right, no move, or move left. --------------
      if (lyhd) 70,80,74
 70   do 72 i = lyhn,imax
        j = imax + lyhn - i
 72     rwork(j) = rwork(j+lyhd)
      go to 80
 74   do 76 i = lyhn,imax
 76     rwork(i) = rwork(i+lyhd)
 80   lyh = lyhn
      iwork(22) = lyh
      if (miter .eq. 0 .or. miter .eq. 3) go to 92
      if (moss .ne. 2) go to 85
c temporarily load ewt if miter = 1 or 2 and moss = 2. -----------------
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 82 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0d0) go to 621
 82     rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
 85   continue
c iprep and prep do sparse matrix preprocessing if miter = 1 or 2. -----
      lsavf = min0(lsavf,lrw)
      lewt = min0(lewt,lrw)
      lacor = min0(lacor,lrw)
      call iprep (neq, y, rwork, iwork(lia), iwork(lja), ipflag, f, jac)
      lenrw = lwm - 1 + lenwk + lrest
      iwork(17) = lenrw
      if (ipflag .ne. -1) iwork(23) = ipian
      if (ipflag .ne. -1) iwork(24) = ipjan
      ipgo = -ipflag + 1
      go to (90, 628, 629, 630, 631, 632, 633), ipgo
 90   iwork(22) = lyh
      if (lenrw .gt. lrw) go to 617
c set flag to signal parameter changes to stode. -----------------------
 92   jstart = -1
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
c the sparse matrix preprocessing (miter = 1 or 2), and the
c calculation of the initial step size.
c the error weights in ewt are inverted after being loaded.
c-----------------------------------------------------------------------
 100  continue
      lyh = lyhn
      iwork(22) = lyh
      tn = t
      nst = 0
      h = 1.0d0
      nnz = 0
      ngp = 0
      nzl = 0
      nzu = 0
c load the initial value vector in yh. ---------------------------------
      do 105 i = 1,n
 105    rwork(i+lyh-1) = y(i)
c initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      call f (neq, t, y, rwork(lf0))
      nfe = 1
c load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 110 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0d0) go to 621
 110    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
      if (miter .eq. 0 .or. miter .eq. 3) go to 120
c iprep and prep do sparse matrix preprocessing if miter = 1 or 2. -----
      lacor = min0(lacor,lrw)
      call iprep (neq, y, rwork, iwork(lia), iwork(lja), ipflag, f, jac)
      lenrw = lwm - 1 + lenwk + lrest
      iwork(17) = lenrw
      if (ipflag .ne. -1) iwork(23) = ipian
      if (ipflag .ne. -1) iwork(24) = ipjan
      ipgo = -ipflag + 1
      go to (115, 628, 629, 630, 631, 632, 633), ipgo
 115  iwork(22) = lyh
      if (lenrw .gt. lrw) go to 617
c check tcrit for legality (itask = 4 or 5). ---------------------------
 120  continue
      if (itask .ne. 4 .and. itask .ne. 5) go to 125
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. 0.0d0) go to 625
      if (h0 .ne. 0.0d0 .and. (t + h0 - tcrit)*h0 .gt. 0.0d0)
     1   h0 = tcrit - t
c initialize all remaining parameters. ---------------------------------
 125  uround = d1mach(4)
      jstart = 0
      if (miter .ne. 0) rwork(lwm) = dsqrt(uround)
      msbj = 50
      nslj = 0
      ccmxj = 0.2d0
      psmall = 1000.0d0*uround
      rbig = 0.01d0/psmall
      nhnil = 0
      nje = 0
      nlu = 0
      nslast = 0
      hu = 0.0d0
      nqu = 0
      ccmax = 0.3d0
      maxcor = 3
      msbp = 20
      mxncf = 10
c-----------------------------------------------------------------------
c the coding below computes the step size, h0, to be attempted on the
c first step, unless the user has supplied a value for this.
c first check that tout - t differs significantly from zero.
c a scalar tolerance quantity tol is computed, as max(rtol(i))
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
c so as to be between 100*uround and 1.0e-3.
c then the computed value h0 is given by..
c                                      neq
c   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
c                                       1
c where   w0     = max ( abs(t), abs(tout) ),
c         f(i)   = i-th component of initial value of f,
c         ywt(i) = ewt(i)/tol  (a weight for y(i)).
c the sign of h0 is inferred from the initial values of tout and t.
c-----------------------------------------------------------------------
      lf0 = lyh + nyh
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
      sum = vnorm (n, rwork(lf0), rwork(lewt))
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
c the call to the one-step core integrator stode.
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
      call xerrwv(50hlsodes-- warning..internal t (=r1) and h (=r2) are,
     1   50, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  60h      such that in the machine, t + h = t on the next step  ,
     1   60, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      (h = step size). solver will continue anyway,
     1   50, 101, 0, 0, 0, 0, 2, tn, h)
      if (nhnil .lt. mxhnil) go to 290
      call xerrwv(50hlsodes-- above warning has been issued i1 times.  ,
     1   50, 102, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      it will not be issued again for this problem,
     1   50, 102, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
 290  continue
c-----------------------------------------------------------------------
c    call stode(neq,y,yh,nyh,yh,ewt,savf,acor,wm,wm,f,jac,prjs,slss)
c-----------------------------------------------------------------------
      call stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),
     1   rwork(lsavf), rwork(lacor), rwork(lwm), rwork(lwm),
     2   f, jac, prjs, slss)
      kgo = 1 - kflag
      go to (300, 530, 540, 550), kgo
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
c the following block handles all successful returns from lsodes.
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
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = nnz
      iwork(20) = ngp
      iwork(21) = nlu
      iwork(25) = nzl
      iwork(26) = nzu
      return
c
 430  ntrep = ntrep + 1
      if (ntrep .lt. 5) return
      call xerrwv(
     1  60hlsodes-- repeated calls with istate = 1 and tout = t (=r1)  ,
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
 500  call xerrwv(50hlsodes-- at current t (=r1), mxstep (=i1) steps   ,
     1   50, 201, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      taken on this call before reaching tout     ,
     1   50, 201, 0, 1, mxstep, 0, 1, tn, 0.0d0)
      istate = -1
      go to 580
c ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  ewti = rwork(lewt+i-1)
      call xerrwv(50hlsodes-- at t (=r1), ewt(i1) has become r2 .le. 0.,
     1   50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 580
c too much accuracy requested for machine precision. -------------------
 520  call xerrwv(50hlsodes-- at t (=r1), too much accuracy requested  ,
     1   50, 203, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      for precision of machine..  see tolsf (=r2) ,
     1   50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580
c kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
 530  call xerrwv(50hlsodes-- at t(=r1) and step size h(=r2), the error,
     1   50, 204, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      test failed repeatedly or with abs(h) = hmin,
     1   50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 560
c kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
 540  call xerrwv(50hlsodes-- at t (=r1) and step size h (=r2), the    ,
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      corrector convergence failed repeatedly     ,
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(30h      or with abs(h) = hmin   ,
     1   30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
      go to 560
c kflag = -3.  fatal error flag returned by prjs or slss (cdrv). -------
 550  call xerrwv(50hlsodes-- at t (=r1) and step size h (=r2), a fatal,
     1   50, 207, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      error flag was returned by cdrv (by way of  ,
     1   50, 207, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(30h      subroutine prjs or slss),
     1   30, 207, 0, 0, 0, 0, 2, tn, h)
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
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = nnz
      iwork(20) = ngp
      iwork(21) = nlu
      iwork(25) = nzl
      iwork(26) = nzu
      return
c-----------------------------------------------------------------------
c block i.
c the following block handles all error returns due to illegal input
c (istate = -3), as detected before calling the core integrator.
c first the error message routine is called.  then if there have been
c 5 consecutive such returns just before this call to the solver,
c the run is halted.
c-----------------------------------------------------------------------
 601  call xerrwv(30hlsodes-- istate (=i1) illegal ,
     1   30, 1, 0, 1, istate, 0, 0, 0.0d0, 0.0d0)
      go to 700
 602  call xerrwv(30hlsodes-- itask (=i1) illegal  ,
     1   30, 2, 0, 1, itask, 0, 0, 0.0d0, 0.0d0)
      go to 700
 603  call xerrwv(50hlsodes-- istate .gt. 1 but lsodes not initialized ,
     1   50, 3, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      go to 700
 604  call xerrwv(30hlsodes-- neq (=i1) .lt. 1     ,
     1   30, 4, 0, 1, neq(1), 0, 0, 0.0d0, 0.0d0)
      go to 700
 605  call xerrwv(50hlsodes-- istate = 3 and neq increased (i1 to i2)  ,
     1   50, 5, 0, 2, n, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 606  call xerrwv(30hlsodes-- itol (=i1) illegal   ,
     1   30, 6, 0, 1, itol, 0, 0, 0.0d0, 0.0d0)
      go to 700
 607  call xerrwv(30hlsodes-- iopt (=i1) illegal   ,
     1   30, 7, 0, 1, iopt, 0, 0, 0.0d0, 0.0d0)
      go to 700
 608  call xerrwv(30hlsodes-- mf (=i1) illegal     ,
     1   30, 8, 0, 1, mf, 0, 0, 0.0d0, 0.0d0)
      go to 700
 609  call xerrwv(30hlsodes-- seth (=r1) .lt. 0.0  ,
     1   30, 9, 0, 0, 0, 0, 1, seth, 0.0d0)
      go to 700
 611  call xerrwv(30hlsodes-- maxord (=i1) .lt. 0  ,
     1   30, 11, 0, 1, maxord, 0, 0, 0.0d0, 0.0d0)
      go to 700
 612  call xerrwv(30hlsodes-- mxstep (=i1) .lt. 0  ,
     1   30, 12, 0, 1, mxstep, 0, 0, 0.0d0, 0.0d0)
      go to 700
 613  call xerrwv(30hlsodes-- mxhnil (=i1) .lt. 0  ,
     1   30, 13, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
      go to 700
 614  call xerrwv(40hlsodes-- tout (=r1) behind t (=r2)      ,
     1   40, 14, 0, 0, 0, 0, 2, tout, t)
      call xerrwv(50h      integration direction is given by h0 (=r1)  ,
     1   50, 14, 0, 0, 0, 0, 1, h0, 0.0d0)
      go to 700
 615  call xerrwv(30hlsodes-- hmax (=r1) .lt. 0.0  ,
     1   30, 15, 0, 0, 0, 0, 1, hmax, 0.0d0)
      go to 700
 616  call xerrwv(30hlsodes-- hmin (=r1) .lt. 0.0  ,
     1   30, 16, 0, 0, 0, 0, 1, hmin, 0.0d0)
      go to 700
 617  call xerrwv(50hlsodes-- rwork length is insufficient to proceed. ,
     1   50, 17, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  60h        length needed is .ge. lenrw (=i1), exceeds lrw (=i2),
     1   60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 618  call xerrwv(50hlsodes-- iwork length is insufficient to proceed. ,
     1   50, 18, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  60h        length needed is .ge. leniw (=i1), exceeds liw (=i2),
     1   60, 18, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
      go to 700
 619  call xerrwv(40hlsodes-- rtol(i1) is r1 .lt. 0.0        ,
     1   40, 19, 0, 1, i, 0, 1, rtoli, 0.0d0)
      go to 700
 620  call xerrwv(40hlsodes-- atol(i1) is r1 .lt. 0.0        ,
     1   40, 20, 0, 1, i, 0, 1, atoli, 0.0d0)
      go to 700
 621  ewti = rwork(lewt+i-1)
      call xerrwv(40hlsodes-- ewt(i1) is r1 .le. 0.0         ,
     1   40, 21, 0, 1, i, 0, 1, ewti, 0.0d0)
      go to 700
 622  call xerrwv(
     1  60hlsodes-- tout (=r1) too close to t(=r2) to start integration,
     1   60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  call xerrwv(
     1  60hlsodes-- itask = i1 and tout (=r1) behind tcur - hu (= r2)  ,
     1   60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624  call xerrwv(
     1  60hlsodes-- itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ,
     1   60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  call xerrwv(
     1  60hlsodes-- itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ,
     1   60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  call xerrwv(50hlsodes-- at start of problem, too much accuracy   ,
     1   50, 26, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  60h      requested for precision of machine..  see tolsf (=r1) ,
     1   60, 26, 0, 0, 0, 0, 1, tolsf, 0.0d0)
      rwork(14) = tolsf
      go to 700
 627  call xerrwv(50hlsodes-- trouble from intdy. itask = i1, tout = r1,
     1   50, 27, 0, 1, itask, 0, 1, tout, 0.0d0)
      go to 700
 628  call xerrwv(
     1  60hlsodes-- rwork length insufficient (for subroutine prep).   ,
     1   60, 28, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  60h        length needed is .ge. lenrw (=i1), exceeds lrw (=i2),
     1   60, 28, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 629  call xerrwv(
     1  60hlsodes-- rwork length insufficient (for subroutine jgroup). ,
     1   60, 29, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  60h        length needed is .ge. lenrw (=i1), exceeds lrw (=i2),
     1   60, 29, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 630  call xerrwv(
     1  60hlsodes-- rwork length insufficient (for subroutine odrv).   ,
     1   60, 30, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  60h        length needed is .ge. lenrw (=i1), exceeds lrw (=i2),
     1   60, 30, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 631  call xerrwv(
     1  60hlsodes-- error from odrv in yale sparse matrix package      ,
     1   60, 31, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      imul = (iys - 1)/n
      irem = iys - imul*n
      call xerrwv(
     1  60h      at t (=r1), odrv returned error flag = i1*neq + i2.   ,
     1   60, 31, 0, 2, imul, irem, 1, tn, 0.0d0)
      go to 700
 632  call xerrwv(
     1  60hlsodes-- rwork length insufficient (for subroutine cdrv).   ,
     1   60, 32, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  60h        length needed is .ge. lenrw (=i1), exceeds lrw (=i2),
     1   60, 32, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 633  call xerrwv(
     1  60hlsodes-- error from cdrv in yale sparse matrix package      ,
     1   60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      imul = (iys - 1)/n
      irem = iys - imul*n
      call xerrwv(
     1  60h      at t (=r1), cdrv returned error flag = i1*neq + i2.   ,
     1   60, 33, 0, 2, imul, irem, 1, tn, 0.0d0)
      if (imul .eq. 2) call xerrwv(
     1  60h        duplicate entry in sparsity structure descriptors   ,
     1   60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      if (imul .eq. 3 .or. imul .eq. 6) call xerrwv(
     1  60h        insufficient storage for nsfc (called by cdrv)      ,
     1   60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
c
 700  if (illin .eq. 5) go to 710
      illin = illin + 1
      istate = -3
      return
 710  call xerrwv(50hlsodes-- repeated occurrences of illegal input    ,
     1   50, 302, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
c
 800  call xerrwv(50hlsodes-- run aborted.. apparent infinite loop     ,
     1   50, 303, 2, 0, 0, 0, 0, 0.0d0, 0.0d0)
      return
c----------------------- end of subroutine lsodes ----------------------
      end
