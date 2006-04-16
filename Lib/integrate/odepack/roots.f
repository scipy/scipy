      subroutine roots (ng, hmin, jflag, x0, x1, g0, g1, gx, x, jroot)
clll. optimize
      integer ng, jflag, jroot
      double precision hmin, x0, x1, g0, g1, gx, x
      dimension g0(ng), g1(ng), gx(ng), jroot(ng)
      integer iownd3, imax, last, idum3
      double precision alpha, x2, rdum3
      common /lsr001/ alpha, x2, rdum3(3),
     1   iownd3(3), imax, last, idum3(4)
c-----------------------------------------------------------------------
c this subroutine finds the leftmost root of a set of arbitrary
c functions gi(x) (i = 1,...,ng) in an interval (x0,x1).  only roots
c of odd multiplicity (i.e. changes of sign of the gi) are found.
c here the sign of x1 - x0 is arbitrary, but is constant for a given
c problem, and -leftmost- means nearest to x0.
c the values of the vector-valued function g(x) = (gi, i=1...ng)
c are communicated through the call sequence of roots.
c the method used is the illinois algorithm.
c
c reference..
c kathie l. hiebert and lawrence f. shampine, implicitly defined
c output points for solutions of ode-s, sandia report sand80-0180,
c february, 1980.
c
c description of parameters.
c
c ng     = number of functions gi, or the number of components of
c          the vector valued function g(x).  input only.
c
c hmin   = resolution parameter in x.  input only.  when a root is
c          found, it is located only to within an error of hmin in x.
c          typically, hmin should be set to something on the order of
c               100 * uround * max(abs(x0),abs(x1)),
c          where uround is the unit roundoff of the machine.
c
c jflag  = integer flag for input and output communication.
c
c          on input, set jflag = 0 on the first call for the problem,
c          and leave it unchanged until the problem is completed.
c          (the problem is completed when jflag .ge. 2 on return.)
c
c          on output, jflag has the following values and meanings..
c          jflag = 1 means roots needs a value of g(x).  set gx = g(x)
c                    and call roots again.
c          jflag = 2 means a root has been found.  the root is
c                    at x, and gx contains g(x).  (actually, x is the
c                    rightmost approximation to the root on an interval
c                    (x0,x1) of size hmin or less.)
c          jflag = 3 means x = x1 is a root, with one or more of the gi
c                    being zero at x1 and no sign changes in (x0,x1).
c                    gx contains g(x) on output.
c          jflag = 4 means no roots (of odd multiplicity) were
c                    found in (x0,x1) (no sign changes).
c
c x0,x1  = endpoints of the interval where roots are sought.
c          x1 and x0 are input when jflag = 0 (first call), and
c          must be left unchanged between calls until the problem is
c          completed.  x0 and x1 must be distinct, but x1 - x0 may be
c          of either sign.  however, the notion of -left- and -right-
c          will be used to mean nearer to x0 or x1, respectively.
c          when jflag .ge. 2 on return, x0 and x1 are output, and
c          are the endpoints of the relevant interval.
c
c g0,g1  = arrays of length ng containing the vectors g(x0) and g(x1),
c          respectively.  when jflag = 0, g0 and g1 are input and
c          none of the g0(i) should be be zero.
c          when jflag .ge. 2 on return, g0 and g1 are output.
c
c gx     = array of length ng containing g(x).  gx is input
c          when jflag = 1, and output when jflag .ge. 2.
c
c x      = independent variable value.  output only.
c          when jflag = 1 on output, x is the point at which g(x)
c          is to be evaluated and loaded into gx.
c          when jflag = 2 or 3, x is the root.
c          when jflag = 4, x is the right endpoint of the interval, x1.
c
c jroot  = integer array of length ng.  output only.
c          when jflag = 2 or 3, jroot indicates which components
c          of g(x) have a root at x.  jroot(i) is 1 if the i-th
c          component has a root, and jroot(i) = 0 otherwise.
c
c note.. this routine uses the common block /lsr001/ to save
c the values of certain variables between calls (own variables).
c-----------------------------------------------------------------------
      integer i, imxold, nxlast
      double precision t2, tmax, zero
      logical zroot, sgnchg, xroot
      data zero/0.0d0/
c
      if (jflag .eq. 1) go to 200
c jflag .ne. 1.  check for change in sign of g or zero at x1. ----------
      imax = 0
      tmax = zero
      zroot = .false.
      do 120 i = 1,ng
        if (dabs(g1(i)) .gt. zero) go to 110
        zroot = .true.
        go to 120
c at this point, g0(i) has been checked and cannot be zero. ------------
 110    if (dsign(1.0d0,g0(i)) .eq. dsign(1.0d0,g1(i))) go to 120
          t2 = dabs(g1(i)/(g1(i)-g0(i)))
          if (t2 .le. tmax) go to 120
            tmax = t2
            imax = i
 120    continue
      if (imax .gt. 0) go to 130
      sgnchg = .false.
      go to 140
 130  sgnchg = .true.
 140  if (.not. sgnchg) go to 400
c there is a sign change.  find the first root in the interval. --------
      xroot = .false.
      nxlast = 0
      last = 1
c
c repeat until the first root in the interval is found.  loop point. ---
 150  continue
      if (xroot) go to 300
      if (nxlast .eq. last) go to 160
      alpha = 1.0d0
      go to 180
 160  if (last .eq. 0) go to 170
      alpha = 0.5d0*alpha
      go to 180
 170  alpha = 2.0d0*alpha
 180  x2 = x1 - (x1-x0)*g1(imax)/(g1(imax) - alpha*g0(imax))
      if ((dabs(x2-x0) .lt. hmin) .and.
     1   (dabs(x1-x0) .gt. 10.0d0*hmin)) x2 = x0 + 0.1d0*(x1-x0)
      jflag = 1
      x = x2
c return to the calling routine to get a value of gx = g(x). -----------
      return
c check to see in which interval g changes sign. -----------------------
 200  imxold = imax
      imax = 0
      tmax = zero
      zroot = .false.
      do 220 i = 1,ng
        if (dabs(gx(i)) .gt. zero) go to 210
        zroot = .true.
        go to 220
c neither g0(i) nor gx(i) can be zero at this point. -------------------
 210    if (dsign(1.0d0,g0(i)) .eq. dsign(1.0d0,gx(i))) go to 220
          t2 = dabs(gx(i)/(gx(i) - g0(i)))
          if (t2 .le. tmax) go to 220
            tmax = t2
            imax = i
 220    continue
      if (imax .gt. 0) go to 230
      sgnchg = .false.
      imax = imxold
      go to 240
 230  sgnchg = .true.
 240  nxlast = last
      if (.not. sgnchg) go to 250
c sign change between x0 and x2, so replace x1 with x2. ----------------
      x1 = x2
      call dcopy (ng, gx, 1, g1, 1)
      last = 1
      xroot = .false.
      go to 270
 250  if (.not. zroot) go to 260
c zero value at x2 and no sign change in (x0,x2), so x2 is a root. -----
      x1 = x2
      call dcopy (ng, gx, 1, g1, 1)
      xroot = .true.
      go to 270
c no sign change between x0 and x2.  replace x0 with x2. ---------------
 260  continue
      call dcopy (ng, gx, 1, g0, 1)
      x0 = x2
      last = 0
      xroot = .false.
 270  if (dabs(x1-x0) .le. hmin) xroot = .true.
      go to 150
c
c return with x1 as the root.  set jroot.  set x = x1 and gx = g1. -----
 300  jflag = 2
      x = x1
      call dcopy (ng, g1, 1, gx, 1)
      do 320 i = 1,ng
        jroot(i) = 0
        if (dabs(g1(i)) .gt. zero) go to 310
          jroot(i) = 1
          go to 320
 310    if (dsign(1.0d0,g0(i)) .ne. dsign(1.0d0,g1(i))) jroot(i) = 1
 320    continue
      return
c
c no sign change in the interval.  check for zero at right endpoint. ---
 400  if (.not. zroot) go to 420
c
c zero value at x1 and no sign change in (x0,x1).  return jflag = 3. ---
      x = x1
      call dcopy (ng, g1, 1, gx, 1)
      do 410 i = 1,ng
        jroot(i) = 0
        if (dabs(g1(i)) .le. zero) jroot (i) = 1
 410  continue
      jflag = 3
      return
c
c no sign changes in this interval.  set x = x1, return jflag = 4. -----
 420  call dcopy (ng, g1, 1, gx, 1)
      x = x1
      jflag = 4
      return
c----------------------- end of subroutine roots -----------------------
      end
