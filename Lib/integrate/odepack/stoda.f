      subroutine stoda (neq, y, yh, nyh, yh1, ewt, savf, acor,
     1   wm, iwm, f, jac, pjac, slvs)
clll. optimize
      external f, jac, pjac, slvs
      integer neq, nyh, iwm
      integer iownd, ialth, ipup, lmax, meo, nqnyh, nslp,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iownd2, icount, irflag, jtyp, mused, mxordn, mxords
      integer i, i1, iredo, iret, j, jb, m, ncf, newq
      integer lm1, lm1p1, lm2, lm2p1, nqm1, nqm2, isav
      double precision y, yh, yh1, ewt, savf, acor, wm, rsav
      double precision conit, crate, el, elco, hold, rmax, tesco,
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision rownd2, pdest, pdlast, ratio, cm1, cm2,
     1   pdnorm
      double precision dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
     1   r, rh, rhdn, rhsm, rhup, told, vmnorm
      double precision alpha, dm1, dm2, exm1, exm2, pdh, pnorm, rate,
     1   rh1, rh1it, rh2, rm, sm1
      dimension neq(1), y(1), yh(nyh,*), yh1(1), ewt(1), savf(1),
     1   acor(1), wm(*), iwm(*), rsav(240), isav(50)
      dimension sm1(12)
      common /ls0001/ conit, crate, el(13), elco(13,12),
     1   hold, rmax, tesco(3,12),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14),
     3   ialth, ipup, lmax, meo, nqnyh, nslp,
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lsa001/ rownd2, pdest, pdlast, ratio, cm1(12), cm2(5),
     1   pdnorm,
     2   iownd2(3), icount, irflag, jtyp, mused, mxordn, mxords
      data sm1/0.5d0, 0.575d0, 0.55d0, 0.45d0, 0.35d0, 0.25d0,
     1   0.20d0, 0.15d0, 0.10d0, 0.075d0, 0.050d0, 0.025d0/
c-----------------------------------------------------------------------
c stoda performs one step of the integration of an initial value
c problem for a system of ordinary differential equations.
c note.. stoda is independent of the value of the iteration method
c indicator miter, when this is .ne. 0, and hence is independent
c of the type of chord method used, or the jacobian structure.
c communication with stoda is done with the following variables..
c
c y      = an array of length .ge. n used as the y argument in
c          all calls to f and jac.
c neq    = integer array containing problem size in neq(1), and
c          passed as the neq argument in all calls to f and jac.
c yh     = an nyh by lmax array containing the dependent variables
c          and their approximate scaled derivatives, where
c          lmax = maxord + 1.  yh(i,j+1) contains the approximate
c          j-th derivative of y(i), scaled by h**j/factorial(j)
c          (j = 0,1,...,nq).  on entry for the first step, the first
c          two columns of yh must be set from the initial values.
c nyh    = a constant integer .ge. n, the first dimension of yh.
c yh1    = a one-dimensional array occupying the same space as yh.
c ewt    = an array of length n containing multiplicative weights
c          for local error measurements.  local errors in y(i) are
c          compared to 1.0/ewt(i) in various error tests.
c savf   = an array of working storage, of length n.
c acor   = a work array of length n, used for the accumulated
c          corrections.  on a successful return, acor(i) contains
c          the estimated one-step local error in y(i).
c wm,iwm = real and integer work arrays associated with matrix
c          operations in chord iteration (miter .ne. 0).
c pjac   = name of routine to evaluate and preprocess jacobian matrix
c          and p = i - h*el0*jac, if a chord method is being used.
c          it also returns an estimate of norm(jac) in pdnorm.
c slvs   = name of routine to solve linear system in chord iteration.
c ccmax  = maximum relative change in h*el0 before pjac is called.
c h      = the step size to be attempted on the next step.
c          h is altered by the error control algorithm during the
c          problem.  h can be either positive or negative, but its
c          sign must remain constant throughout the problem.
c hmin   = the minimum absolute value of the step size h to be used.
c hmxi   = inverse of the maximum absolute value of h to be used.
c          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
c          hmin and hmxi may be changed at any time, but will not
c          take effect until the next change of h is considered.
c tn     = the independent variable. tn is updated on each step taken.
c jstart = an integer used for input only, with the following
c          values and meanings..
c               0  perform the first step.
c           .gt.0  take a new step continuing from the last.
c              -1  take the next step with a new value of h,
c                    n, meth, miter, and/or matrix parameters.
c              -2  take the next step with a new value of h,
c                    but with other inputs unchanged.
c          on return, jstart is set to 1 to facilitate continuation.
c kflag  = a completion code with the following meanings..
c               0  the step was succesful.
c              -1  the requested error could not be achieved.
c              -2  corrector convergence could not be achieved.
c              -3  fatal error in pjac or slvs.
c          a return with kflag = -1 or -2 means either
c          abs(h) = hmin or 10 consecutive failures occurred.
c          on a return with kflag negative, the values of tn and
c          the yh array are as of the beginning of the last
c          step, and h is the last step size attempted.
c maxord = the maximum order of integration method to be allowed.
c maxcor = the maximum number of corrector iterations allowed.
c msbp   = maximum number of steps between pjac calls (miter .gt. 0).
c mxncf  = maximum number of convergence failures allowed.
c meth   = current method.
c          meth = 1 means adams method (nonstiff)
c          meth = 2 means bdf method (stiff)
c          meth may be reset by stoda.
c miter  = corrector iteration method.
c          miter = 0 means functional iteration.
c          miter = jt .gt. 0 means a chord iteration corresponding
c          to jacobian type jt.  (the lsoda argument jt is
c          communicated here as jtyp, but is not used in stoda
c          except to load miter following a method switch.)
c          miter may be reset by stoda.
c n      = the number of first-order differential equations.
c-----------------------------------------------------------------------
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0d0
      if (jstart .gt. 0) go to 200
      if (jstart .eq. -1) go to 100
      if (jstart .eq. -2) go to 160
c-----------------------------------------------------------------------
c on the first call, the order is set to 1, and other variables are
c initialized.  rmax is the maximum ratio by which h can be increased
c in a single step.  it is initially 1.e4 to compensate for the small
c initial h, but then is normally equal to 10.  if a failure
c occurs (in corrector convergence or error test), rmax is set at 2
c for the next increase.
c cfode is called to get the needed coefficients for both methods.
c-----------------------------------------------------------------------
      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0d0
      rc = 0.0d0
      el0 = 1.0d0
      crate = 0.7d0
      hold = h
      nslp = 0
      ipup = miter
      iret = 3
c initialize switching parameters.  meth = 1 is assumed initially. -----
      icount = 20
      irflag = 0
      pdest = 0.0d0
      pdlast = 0.0d0
      ratio = 5.0d0
      call cfode (2, elco, tesco)
      do 10 i = 1,5
 10     cm2(i) = tesco(2,i)*elco(i+1,i)
      call cfode (1, elco, tesco)
      do 20 i = 1,12
 20     cm1(i) = tesco(2,i)*elco(i+1,i)
      go to 150
c-----------------------------------------------------------------------
c the following block handles preliminaries needed when jstart = -1.
c ipup is set to miter to force a matrix update.
c if an order increase is about to be considered (ialth = 1),
c ialth is reset to 2 to postpone consideration one more step.
c if the caller has changed meth, cfode is called to reset
c the coefficients of the method.
c if h is to be changed, yh must be rescaled.
c if h or meth is being changed, ialth is reset to l = nq + 1
c to prevent further changes in h for that many steps.
c-----------------------------------------------------------------------
 100  ipup = miter
      lmax = maxord + 1
      if (ialth .eq. 1) ialth = 2
      if (meth .eq. mused) go to 160
      call cfode (meth, elco, tesco)
      ialth = l
      iret = 1
c-----------------------------------------------------------------------
c the el vector and related constants are reset
c whenever the order nq is changed, or at the start of the problem.
c-----------------------------------------------------------------------
 150  do 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5d0/dfloat(nq+2)
      go to (160, 170, 200), iret
c-----------------------------------------------------------------------
c if h is being changed, the h ratio rh is checked against
c rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
c l = nq + 1 to prevent a change of h for that many steps, unless
c forced by a convergence or error test failure.
c-----------------------------------------------------------------------
 160  if (h .eq. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = dmax1(rh,hmin/dabs(h))
 175  rh = dmin1(rh,rmax)
      rh = rh/dmax1(1.0d0,dabs(h)*hmxi*rh)
c-----------------------------------------------------------------------
c if meth = 1, also restrict the new step size by the stability region.
c if this reduces h, set irflag to 1 so that if there are roundoff
c problems later, we can assume that is the cause of the trouble.
c-----------------------------------------------------------------------
      if (meth .eq. 2) go to 178
      irflag = 0
      pdh = dmax1(dabs(h)*pdlast,0.000001d0)
      if (rh*pdh*1.00001d0 .lt. sm1(nq)) go to 178
      rh = sm1(nq)/pdh
      irflag = 1
 178  continue
      r = 1.0d0
      do 180 j = 2,l
        r = r*rh
        do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 690
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the yh array by the pascal triangle matrix.
c rc is the ratio of new to old values of the coefficient  h*el(1).
c when rc differs from 1 by more than ccmax, ipup is set to miter
c to force pjac to be called, if a jacobian is involved.
c in any case, pjac is called at least every msbp steps.
c-----------------------------------------------------------------------
 200  if (dabs(rc-1.0d0) .gt. ccmax) ipup = miter
      if (nst .ge. nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do 215 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    continue
      pnorm = vmnorm (n, yh1, ewt)
c-----------------------------------------------------------------------
c up to maxcor corrector iterations are taken.  a convergence test is
c made on the r.m.s. norm of each correction, weighted by the error
c weight vector ewt.  the sum of the corrections is accumulated in the
c vector acor(i).  the yh array is not altered in the corrector loop.
c-----------------------------------------------------------------------
 220  m = 0
      rate = 0.0d0
      del = 0.0d0
      do 230 i = 1,n
 230    y(i) = yh(i,1)
      call srcma (rsav, isav, 1)
      call f (neq, tn, y, savf)
      call srcma (rsav, isav, 2)
      nfe = nfe + 1
      if (ipup .le. 0) go to 250
c-----------------------------------------------------------------------
c if indicated, the matrix p = i - h*el(1)*j is reevaluated and
c preprocessed before starting the corrector iteration.  ipup is set
c to 0 as an indicator that this has been done.
c-----------------------------------------------------------------------
      call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
      ipup = 0
      rc = 1.0d0
      nslp = nst
      crate = 0.7d0
      if (ierpj .ne. 0) go to 430
 250  do 260 i = 1,n
 260    acor(i) = 0.0d0
 270  if (miter .ne. 0) go to 350
c-----------------------------------------------------------------------
c in the case of functional iteration, update y directly from
c the result of the last function evaluation.
c-----------------------------------------------------------------------
      do 290 i = 1,n
        savf(i) = h*savf(i) - yh(i,2)
 290    y(i) = savf(i) - acor(i)
      del = vmnorm (n, y, ewt)
      do 300 i = 1,n
        y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400
c-----------------------------------------------------------------------
c in the case of the chord method, compute the corrector error,
c and solve the linear system with that as right-hand side and
c p as coefficient matrix.
c-----------------------------------------------------------------------
 350  do 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))
      call slvs (wm, iwm, y, savf)
      if (iersl .lt. 0) go to 430
      if (iersl .gt. 0) go to 410
      del = vmnorm (n, y, ewt)
      do 380 i = 1,n
        acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)
c-----------------------------------------------------------------------
c test for convergence.  if m.gt.0, an estimate of the convergence
c rate constant is stored in crate, and this is used in the test.
c
c we first check for a change of iterates that is the size of
c roundoff error.  if this occurs, the iteration has converged, and a
c new rate estimate is not formed.
c in all other cases, force at least two iterations to estimate a
c local lipschitz constant estimate for adams methods.
c on convergence, form pdest = local maximum lipschitz constant
c estimate.  pdlast is the most recent nonzero estimate.
c-----------------------------------------------------------------------
 400  continue
      if (del .le. 100.0d0*pnorm*uround) go to 450
      if (m .eq. 0 .and. meth .eq. 1) go to 405
      if (m .eq. 0) go to 402
      rm = 1024.0d0
      if (del .le. 1024.0d0*delp) rm = del/delp
      rate = dmax1(rate,rm)
      crate = dmax1(0.2d0*crate,rm)
 402  dcon = del*dmin1(1.0d0,1.5d0*crate)/(tesco(2,nq)*conit)
      if (dcon .gt. 1.0d0) go to 405
      pdest = dmax1(pdest,rate/dabs(h*el(1)))
      if (pdest .ne. 0.0d0) pdlast = pdest
      go to 450
 405  continue
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. 2.0d0*delp) go to 410
      delp = del
      call srcma (rsav, isav, 1)
      call f (neq, tn, y, savf)
      call srcma (rsav, isav, 2)
      nfe = nfe + 1
      go to 270
c-----------------------------------------------------------------------
c the corrector iteration failed to converge.
c if miter .ne. 0 and the jacobian is out of date, pjac is called for
c the next try.  otherwise the yh array is retracted to its values
c before prediction, and h is reduced, if possible.  if h cannot be
c reduced or mxncf failures have occurred, exit with kflag = -2.
c-----------------------------------------------------------------------
 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0d0
      tn = told
      i1 = nqnyh + 1
      do 445 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    continue
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
      if (dabs(h) .le. hmin*1.00001d0) go to 670
      if (ncf .eq. mxncf) go to 670
      rh = 0.25d0
      ipup = miter
      iredo = 1
      go to 170
c-----------------------------------------------------------------------
c the corrector has converged.  jcur is set to 0
c to signal that the jacobian involved may need updating later.
c the local error test is made and control passes to statement 500
c if it fails.
c-----------------------------------------------------------------------
 450  jcur = 0
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vmnorm (n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0d0) go to 500
c-----------------------------------------------------------------------
c after a successful step, update the yh array.
c decrease icount by 1, and if it is -1, consider switching methods.
c if a method switch is made, reset various parameters,
c rescale the yh array, and exit.  if there is no switch,
c consider changing h if ialth = 1.  otherwise decrease ialth by 1.
c if ialth is then 1 and nq .lt. maxord, then acor is saved for
c use in a possible order increase on the next step.
c if a change in h is considered, an increase or decrease in order
c by one is considered also.  a change in h is made only if it is by a
c factor of at least 1.1.  if not, ialth is set to 3 to prevent
c testing for that many steps.
c-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      mused = meth
      do 460 j = 1,l
        do 460 i = 1,n
 460      yh(i,j) = yh(i,j) + el(j)*acor(i)
      icount = icount - 1
      if (icount .ge. 0) go to 488
      if (meth .eq. 2) go to 480
c-----------------------------------------------------------------------
c we are currently using an adams method.  consider switching to bdf.
c if the current order is greater than 5, assume the problem is
c not stiff, and skip this section.
c if the lipschitz constant and error estimate are not polluted
c by roundoff, go to 470 and perform the usual test.
c otherwise, switch to the bdf methods if the last step was
c restricted to insure stability (irflag = 1), and stay with adams
c method if not.  when switching to bdf with polluted error estimates,
c in the absence of other information, double the step size.
c
c when the estimates are ok, we make the usual test by computing
c the step size we could have (ideally) used on this step,
c with the current (adams) method, and also that for the bdf.
c if nq .gt. mxords, we consider changing to order mxords on switching.
c compare the two step sizes to decide whether to switch.
c the step size advantage must be at least ratio = 5 to switch.
c-----------------------------------------------------------------------
      if (nq .gt. 5) go to 488
      if (dsm .gt. 100.0d0*pnorm*uround .and. pdest .ne. 0.0d0)
     1   go to 470
      if (irflag .eq. 0) go to 488
      rh2 = 2.0d0
      nqm2 = min0(nq,mxords)
      go to 478
 470  continue
      exsm = 1.0d0/dfloat(l)
      rh1 = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
      rh1it = 2.0d0*rh1
      pdh = pdlast*dabs(h)
      if (pdh*rh1 .gt. 0.00001d0) rh1it = sm1(nq)/pdh
      rh1 = dmin1(rh1,rh1it)
      if (nq .le. mxords) go to 474
         nqm2 = mxords
         lm2 = mxords + 1
         exm2 = 1.0d0/dfloat(lm2)
         lm2p1 = lm2 + 1
         dm2 = vmnorm (n, yh(1,lm2p1), ewt)/cm2(mxords)
         rh2 = 1.0d0/(1.2d0*dm2**exm2 + 0.0000012d0)
         go to 476
 474  dm2 = dsm*(cm1(nq)/cm2(nq))
      rh2 = 1.0d0/(1.2d0*dm2**exsm + 0.0000012d0)
      nqm2 = nq
 476  continue
      if (rh2 .lt. ratio*rh1) go to 488
c the switch test passed.  reset relevant quantities for bdf. ----------
 478  rh = rh2
      icount = 20
      meth = 2
      miter = jtyp
      pdlast = 0.0d0
      nq = nqm2
      l = nq + 1
      go to 170
c-----------------------------------------------------------------------
c we are currently using a bdf method.  consider switching to adams.
c compute the step size we could have (ideally) used on this step,
c with the current (bdf) method, and also that for the adams.
c if nq .gt. mxordn, we consider changing to order mxordn on switching.
c compare the two step sizes to decide whether to switch.
c the step size advantage must be at least 5/ratio = 1 to switch.
c if the step size for adams would be so small as to cause
c roundoff pollution, we stay with bdf.
c-----------------------------------------------------------------------
 480  continue
      exsm = 1.0d0/dfloat(l)
      if (mxordn .ge. nq) go to 484
         nqm1 = mxordn
         lm1 = mxordn + 1
         exm1 = 1.0d0/dfloat(lm1)
         lm1p1 = lm1 + 1
         dm1 = vmnorm (n, yh(1,lm1p1), ewt)/cm1(mxordn)
         rh1 = 1.0d0/(1.2d0*dm1**exm1 + 0.0000012d0)
         go to 486
 484  dm1 = dsm*(cm2(nq)/cm1(nq))
      rh1 = 1.0d0/(1.2d0*dm1**exsm + 0.0000012d0)
      nqm1 = nq
      exm1 = exsm
 486  rh1it = 2.0d0*rh1
      pdh = pdnorm*dabs(h)
      if (pdh*rh1 .gt. 0.00001d0) rh1it = sm1(nqm1)/pdh
      rh1 = dmin1(rh1,rh1it)
      rh2 = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
      if (rh1*ratio .lt. 5.0d0*rh2) go to 488
      alpha = dmax1(0.001d0,rh1)
      dm1 = (alpha**exm1)*dm1
      if (dm1 .le. 1000.0d0*uround*pnorm) go to 488
c the switch test passed.  reset relevant quantities for adams. --------
      rh = rh1
      icount = 20
      meth = 1
      miter = 0
      pdlast = 0.0d0
      nq = nqm1
      l = nq + 1
      go to 170
c
c no method switch is being made.  do the usual step/order selection. --
 488  continue
      ialth = ialth - 1
      if (ialth .eq. 0) go to 520
      if (ialth .gt. 1) go to 700
      if (l .eq. lmax) go to 700
      do 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 700
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c restore tn and the yh array to their previous values, and prepare
c to try the step again.  compute the optimum step size for this or
c one lower order.  after 2 or more failures, h is forced to decrease
c by a factor of 0.2 or less.
c-----------------------------------------------------------------------
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    continue
      rmax = 2.0d0
      if (dabs(h) .le. hmin*1.00001d0) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0d0
      go to 540
c-----------------------------------------------------------------------
c regardless of the success or failure of the step, factors
c rhdn, rhsm, and rhup are computed, by which h could be multiplied
c at order nq - 1, order nq, or order nq + 1, respectively.
c in the case of failure, rhup = 0.0 to avoid an order increase.
c the largest of these is determined and the new order chosen
c accordingly.  if the order is to be increased, we compute one
c additional scaled derivative.
c-----------------------------------------------------------------------
 520  rhup = 0.0d0
      if (l .eq. lmax) go to 540
      do 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vmnorm (n, savf, ewt)/tesco(3,nq)
      exup = 1.0d0/dfloat(l+1)
      rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)
 540  exsm = 1.0d0/dfloat(l)
      rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
      rhdn = 0.0d0
      if (nq .eq. 1) go to 550
      ddn = vmnorm (n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0d0/dfloat(nq)
      rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
c if meth = 1, limit rh according to the stability region also. --------
 550  if (meth .eq. 2) go to 560
      pdh = dmax1(dabs(h)*pdlast,0.000001d0)
      if (l .lt. lmax) rhup = dmin1(rhup,sm1(l)/pdh)
      rhsm = dmin1(rhsm,sm1(nq)/pdh)
      if (nq .gt. 1) rhdn = dmin1(rhdn,sm1(nq-1)/pdh)
      pdest = 0.0d0
 560  if (rhsm .ge. rhup) go to 570
      if (rhup .gt. rhdn) go to 590
      go to 580
 570  if (rhsm .lt. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      if (kflag .lt. 0 .and. rh .gt. 1.0d0) rh = 1.0d0
      go to 620
 590  newq = l
      rh = rhup
      if (rh .lt. 1.1d0) go to 610
      r = el(l)/dfloat(l)
      do 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 700
c if meth = 1 and h is restricted by stability, bypass 10 percent test.
 620  if (meth .eq. 2) go to 622
      if (rh*pdh*1.00001d0 .ge. sm1(newq)) go to 625
 622  if (kflag .eq. 0 .and. rh .lt. 1.1d0) go to 610
 625  if (kflag .le. -2) rh = dmin1(rh,0.2d0)
c-----------------------------------------------------------------------
c if there is a change of order, reset nq, l, and the coefficients.
c in any case h is reset according to rh and the yh array is rescaled.
c then exit from 690 if the step was ok, or redo the step otherwise.
c-----------------------------------------------------------------------
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
c-----------------------------------------------------------------------
c control reaches this section if 3 or more failures have occured.
c if 10 failures have occurred, exit with kflag = -1.
c it is assumed that the derivatives that have accumulated in the
c yh array have errors of the wrong order.  hence the first
c derivative is recomputed, and the order is set to 1.  then
c h is reduced by a factor of 10, and the step is retried,
c until it succeeds or h reaches hmin.
c-----------------------------------------------------------------------
 640  if (kflag .eq. -10) go to 660
      rh = 0.1d0
      rh = dmax1(hmin/dabs(h),rh)
      h = h*rh
      do 645 i = 1,n
 645    y(i) = yh(i,1)
      call srcma (rsav, isav, 1)
      call f (neq, tn, y, savf)
      call srcma (rsav, isav, 2)
      nfe = nfe + 1
      do 650 i = 1,n
 650    yh(i,2) = h*savf(i)
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
c-----------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
 660  kflag = -1
      go to 720
 670  kflag = -2
      go to 720
 680  kflag = -3
      go to 720
 690  rmax = 10.0d0
 700  r = 1.0d0/tesco(2,nqu)
      do 710 i = 1,n
 710    acor(i) = acor(i)*r
 720  hold = h
      jstart = 1
      return
c----------------------- end of subroutine stoda -----------------------
      end
