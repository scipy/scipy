      subroutine prepji (neq, y, yh, nyh, ewt, rtem, savr, s, wm, iwm,
     1   res, jac, adda)
clll. optimize
      external res, jac, adda
      integer neq, nyh, iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nre, nje, nqu
      integer i, i1, i2, ier, ii, ires, j, j1, jj, lenp,
     1   mba, mband, meb1, meband, ml, ml3, mu
      double precision y, yh, ewt, rtem, savr, s, wm
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision con, fac, hl0, r, srur, yi, yj, yjj
      dimension neq(1), y(1), yh(nyh,*), ewt(1), rtem(1),
     1   s(1), savr(1), wm(*), iwm(*)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nre, nje, nqu
c-----------------------------------------------------------------------
c prepji is called by stodi to compute and process the matrix
c p = a - h*el(1)*j , where j is an approximation to the jacobian dr/dy,
c where r = g(t,y) - a(t,y)*s. here j is computed by the user-supplied
c routine jac if miter = 1 or 4, or by finite differencing if miter =
c 2 or 5. j is stored in wm, rescaled, and adda is called to generate
c p. p is then subjected to lu decomposition in preparation
c for later solution of linear systems with p as coefficient
c matrix.  this is done by dgefa if miter = 1 or 2, and by
c dgbfa if miter = 4 or 5.
c
c in addition to variables described previously, communication
c with prepji uses the following..
c y     = array containing predicted values on entry.
c rtem  = work array of length n (acor in stodi).
c savr  = array used for output only.  on output it contains the
c         residual evaluated at current values of t and y.
c s     = array containing predicted values of dy/dt (savf in stodi).
c wm    = real work space for matrices.  on output it contains the
c         lu decomposition of p.
c         storage of matrix elements starts at wm(3).
c         wm also contains the following matrix-related data..
c         wm(1) = sqrt(uround), used in numerical jacobian increments.
c iwm   = integer work space containing pivot information, starting at
c         iwm(21).  iwm also contains the band parameters
c         ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c el0   = el(1) (input).
c ierpj = output error flag.
c         = 0 if no trouble occurred,
c         = 1 if the p matrix was found to be singular,
c         = ires (= 2 or 3) if res returned ires = 2 or 3.
c jcur  = output flag = 1 to indicate that the jacobian matrix
c         (or approximation) is now current.
c this routine also uses the common variables el0, h, tn, uround,
c miter, n, nre, and nje.
c-----------------------------------------------------------------------
      nje = nje + 1
      hl0 = h*el0
      ierpj = 0
      jcur = 1
      go to (100, 200, 300, 400, 500), miter
c if miter = 1, call res, then jac, and multiply by scalar. ------------
 100  ires = 1
      call res (neq, tn, y, s, savr, ires)
      nre = nre + 1
      if (ires .gt. 1) go to 600
      lenp = n*n
      do 110 i = 1,lenp
 110    wm(i+2) = 0.0d0
      call jac ( neq, tn, y, s, 0, 0, wm(3), n )
      con = -hl0
      do 120 i = 1,lenp
 120    wm(i+2) = wm(i+2)*con
      go to 240
c if miter = 2, make n + 1 calls to res to approximate j. --------------
 200  continue
      ires = -1
      call res (neq, tn, y, s, savr, ires)
      nre = nre + 1
      if (ires .gt. 1) go to 600
      srur = wm(1)
      j1 = 2
      do 230 j = 1,n
        yj = y(j)
        r = dmax1(srur*dabs(yj),0.01d0/ewt(j))
        y(j) = y(j) + r
        fac = -hl0/r
        call res ( neq, tn, y, s, rtem, ires )
        nre = nre + 1
        if (ires .gt. 1) go to 600
        do 220 i = 1,n
 220      wm(i+j1) = (rtem(i) - savr(i))*fac
        y(j) = yj
        j1 = j1 + n
 230    continue
      ires = 1
      call res (neq, tn, y, s, savr, ires)
      nre = nre + 1
      if (ires .gt. 1) go to 600
c add matrix a. --------------------------------------------------------
 240  continue
      call adda(neq, tn, y, 0, 0, wm(3), n)
c do lu decomposition on p. --------------------------------------------
      call dgefa (wm(3), n, n, iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c dummy section for miter = 3
 300  return
c if miter = 4, call res, then jac, and multiply by scalar. ------------
 400  ires = 1
      call res (neq, tn, y, s, savr, ires)
      nre = nre + 1
      if (ires .gt. 1) go to 600
      ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
      do 410 i = 1,lenp
 410    wm(i+2) = 0.0d0
      call jac ( neq, tn, y, s, ml, mu, wm(ml3), meband)
      con = -hl0
      do 420 i = 1,lenp
 420    wm(i+2) = wm(i+2)*con
      go to 570
c if miter = 5, make ml + mu + 2 calls to res to approximate j. --------
 500  continue
      ires = -1
      call res (neq, tn, y, s, savr, ires)
      nre = nre + 1
      if (ires .gt. 1) go to 600
      ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      mba = min0(mband,n)
      meband = mband + ml
      meb1 = meband - 1
      srur = wm(1)
      do 560 j = 1,mba
        do 530 i = j,n,mband
          yi = y(i)
          r = dmax1(srur*dabs(yi),0.01d0/ewt(i))
 530      y(i) = y(i) + r
        call res ( neq, tn, y, s, rtem, ires)
        nre = nre + 1
        if (ires .gt. 1) go to 600
        do 550 jj = j,n,mband
          y(jj) = yh(jj,1)
          yjj = y(jj)
          r = dmax1(srur*dabs(yjj),0.01d0/ewt(jj))
          fac = -hl0/r
          i1 = max0(jj-mu,1)
          i2 = min0(jj+ml,n)
          ii = jj*meb1 - ml + 2
          do 540 i = i1,i2
 540        wm(ii+i) = (rtem(i) - savr(i))*fac
 550      continue
 560    continue
      ires = 1
      call res (neq, tn, y, s, savr, ires)
      nre = nre + 1
      if (ires .gt. 1) go to 600
c add matrix a. --------------------------------------------------------
 570  continue
      call adda(neq, tn, y, ml, mu, wm(ml3), meband)
c do lu decomposition of p. --------------------------------------------
      call dgbfa (wm(3), meband, n, ml, mu, iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c error return for ires = 2 or ires = 3 return from res. ---------------
 600  ierpj = ires
      return
c----------------------- end of subroutine prepji ----------------------
      end
