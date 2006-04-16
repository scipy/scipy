      subroutine rchek (job, g, neq, y, yh, nyh, g0, g1, gx, jroot, irt)
clll. optimize
      external g
      integer job, neq, nyh, jroot, irt
      double precision y, yh, g0, g1, gx
      dimension neq(1), y(1), yh(nyh,*), g0(1), g1(1), gx(1), jroot(1)
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iownd3, iownr3, irfnd, itaskc, ngc, nge
      integer i, iflag, jflag
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision rownr3, t0, tlast, toutc
      double precision hming, t1, temp1, temp2, x
      logical zroot
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lsr001/ rownr3(2), t0, tlast, toutc,
     1   iownd3(3), iownr3(2), irfnd, itaskc, ngc, nge
c-----------------------------------------------------------------------
c this routine checks for the presence of a root in the
c vicinity of the current t, in a manner depending on the
c input flag job.  it calls subroutine roots to locate the root
c as precisely as possible.
c
c in addition to variables described previously, rchek
c uses the following for communication..
c job    = integer flag indicating type of call..
c          job = 1 means the problem is being initialized, and rchek
c                  is to look for a root at or very near the initial t.
c          job = 2 means a continuation call to the solver was just
c                  made, and rchek is to check for a root in the
c                  relevant part of the step last taken.
c          job = 3 means a successful step was just taken, and rchek
c                  is to look for a root in the interval of the step.
c g0     = array of length ng, containing the value of g at t = t0.
c          g0 is input for job .ge. 2 and on output in all cases.
c g1,gx  = arrays of length ng for work space.
c irt    = completion flag..
c          irt = 0  means no root was found.
c          irt = -1 means job = 1 and a root was found too near to t.
c          irt = 1  means a legitimate root was found (job = 2 or 3).
c                   on return, t0 is the root location, and y is the
c                   corresponding solution vector.
c t0     = value of t at one endpoint of interval of interest.  only
c          roots beyond t0 in the direction of integration are sought.
c          t0 is input if job .ge. 2, and output in all cases.
c          t0 is updated by rchek, whether a root is found or not.
c tlast  = last value of t returned by the solver (input only).
c toutc  = copy of tout (input only).
c irfnd  = input flag showing whether the last step taken had a root.
c          irfnd = 1 if it did, = 0 if not.
c itaskc = copy of itask (input only).
c ngc    = copy of ng (input only).
c-----------------------------------------------------------------------
c
      irt = 0
      do 10 i = 1,ngc
 10     jroot(i) = 0
      hming = (dabs(tn) + dabs(h))*uround*100.0d0
c
      go to (100, 200, 300), job
c
c evaluate g at initial t, and check for zero values. ------------------
 100  continue
      t0 = tn
      call g (neq, t0, y, ngc, g0)
      nge = 1
      zroot = .false.
      do 110 i = 1,ngc
 110    if (dabs(g0(i)) .le. 0.0d0) zroot = .true.
      if (.not. zroot) go to 190
c g has a zero at t.  look at g at t + (small increment). --------------
      temp1 = dsign(hming,h)
      t0 = t0 + temp1
      temp2 = temp1/h
      do 120 i = 1,n
 120    y(i) = y(i) + temp2*yh(i,2)
      call g (neq, t0, y, ngc, g0)
      nge = nge + 1
      zroot = .false.
      do 130 i = 1,ngc
 130    if (dabs(g0(i)) .le. 0.0d0) zroot = .true.
      if (.not. zroot) go to 190
c g has a zero at t and also close to t.  take error return. -----------
      irt = -1
      return
c
 190  continue
      return
c
c
 200  continue
      if (irfnd .eq. 0) go to 260
c if a root was found on the previous step, evaluate g0 = g(t0). -------
      call intdy (t0, 0, yh, nyh, y, iflag)
      call g (neq, t0, y, ngc, g0)
      nge = nge + 1
      zroot = .false.
      do 210 i = 1,ngc
 210    if (dabs(g0(i)) .le. 0.0d0) zroot = .true.
      if (.not. zroot) go to 260
c g has a zero at t0.  look at g at t + (small increment). -------------
      temp1 = dsign(hming,h)
      t0 = t0 + temp1
      if ((t0 - tn)*h .lt. 0.0d0) go to 230
      temp2 = temp1/h
      do 220 i = 1,n
 220    y(i) = y(i) + temp2*yh(i,2)
      go to 240
 230  call intdy (t0, 0, yh, nyh, y, iflag)
 240  call g (neq, t0, y, ngc, g0)
      nge = nge + 1
      zroot = .false.
      do 250 i = 1,ngc
        if (dabs(g0(i)) .gt. 0.0d0) go to 250
        jroot(i) = 1
        zroot = .true.
 250    continue
      if (.not. zroot) go to 260
c g has a zero at t0 and also close to t0.  return root. ---------------
      irt = 1
      return
c     here, g0 does not have a root
c g0 has no zero components.  proceed to check relevant interval. ------
 260  if (tn .eq. tlast) go to 390
c
 300  continue
c set t1 to tn or toutc, whichever comes first, and get g at t1. -------
      if (itaskc.eq.2 .or. itaskc.eq.3 .or. itaskc.eq.5) go to 310
      if ((toutc - tn)*h .ge. 0.0d0) go to 310
      t1 = toutc
      if ((t1 - t0)*h .le. 0.0d0) go to 390
      call intdy (t1, 0, yh, nyh, y, iflag)
      go to 330
 310  t1 = tn
      do 320 i = 1,n
 320    y(i) = yh(i,1)
 330  call g (neq, t1, y, ngc, g1)
      nge = nge + 1
c call roots to search for root in interval from t0 to t1. -------------
      jflag = 0
 350  continue
      call roots (ngc, hming, jflag, t0, t1, g0, g1, gx, x, jroot)
      if (jflag .gt. 1) go to 360
      call intdy (x, 0, yh, nyh, y, iflag)
      call g (neq, x, y, ngc, gx)
      nge = nge + 1
      go to 350
 360  t0 = x
      call dcopy (ngc, gx, 1, g0, 1)
      if (jflag .eq. 4) go to 390
c found a root.  interpolate to x and return. --------------------------
      call intdy (x, 0, yh, nyh, y, iflag)
      irt = 1
      return
c
 390  continue
      return
c----------------------- end of subroutine rchek -----------------------
      end
