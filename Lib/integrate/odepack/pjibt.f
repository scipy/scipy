      subroutine pjibt (neq, y, yh, nyh, ewt, rtem, savr, s, wm, iwm,
     1   res, jac, adda)
clll. optimize
      external res, jac, adda
      integer neq, nyh, iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nre, nje, nqu
      integer i, ier, iia, iib, iic, ipa, ipb, ipc, ires, j, j1, j2,
     1   k, k1, lenp, lblox, lpb, lpc, mb, mbsq, mwid, nb
      double precision y, yh, ewt, rtem, savr, s, wm
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision con, fac, hl0, r, srur
      dimension neq(1), y(1), yh(nyh,1), ewt(1), rtem(1),
     1   s(1), savr(1), wm(*), iwm(*)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nre, nje, nqu
c-----------------------------------------------------------------------
c pjibt is called by stodi to compute and process the matrix
c p = a - h*el(1)*j , where j is an approximation to the jacobian dr/dy,
c and r = g(t,y) - a(t,y)*s.  here j is computed by the user-supplied
c routine jac if miter = 1, or by finite differencing if miter = 2.
c j is stored in wm, rescaled, and adda is called to generate p.
c p is then subjected to lu decomposition by decbt in preparation
c for later solution of linear systems with p as coefficient matrix.
c
c in addition to variables described previously, communication
c with pjibt uses the following..
c y     = array containing predicted values on entry.
c rtem  = work array of length n (acor in stodi).
c savr  = array used for output only.  on output it contains the
c         residual evaluated at current values of t and y.
c s     = array containing predicted values of dy/dt (savf in stodi).
c wm    = real work space for matrices.  on output it contains the
c         lu decomposition of p.
c         storage of matrix elements starts at wm(3).
c         wm also contains the following matrix-related data..
c         wm(1) = dsqrt(uround), used in numerical jacobian increments.
c iwm   = integer work space containing pivot information, starting at
c         iwm(21).  iwm also contains block structure parameters
c         mb = iwm(1) and nb = iwm(2).
c el0   = el(1) (input).
c ierpj = output error flag.
c         = 0 if no trouble occurred,
c         = 1 if the p matrix was found to be unfactorable,
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
      mb = iwm(1)
      nb = iwm(2)
      mbsq = mb*mb
      lblox = mbsq*nb
      lpb = 3 + lblox
      lpc = lpb + lblox
      lenp = 3*lblox
      go to (100, 200), miter
c if miter = 1, call res, then jac, and multiply by scalar. ------------
 100  ires = 1
      call res (neq, tn, y, s, savr, ires)
      nre = nre + 1
      if (ires .gt. 1) go to 600
      do 110 i = 1,lenp
 110    wm(i+2) = 0.0d0
      call jac (neq, tn, y, s, mb, nb, wm(3), wm(lpb), wm(lpc))
      con = -hl0
      do 120 i = 1,lenp
 120    wm(i+2) = wm(i+2)*con
      go to 260
c
c if miter = 2, make 3*mb + 1 calls to res to approximate j. -----------
 200  continue
      ires = -1
      call res (neq, tn, y, s, savr, ires)
      nre = nre + 1
      if (ires .gt. 1) go to 600
      mwid = 3*mb
      srur = wm(1)
      do 205 i = 1,lenp
 205    wm(2+i) = 0.0d0
      do 250 k = 1,3
        do 240 j = 1,mb
c         increment y(i) for group of column indices, and call res. ----
          j1 = j+(k-1)*mb
          do 210 i = j1,n,mwid
            r = dmax1(srur*dabs(y(i)),0.01d0/ewt(i))
            y(i) = y(i) + r
 210      continue
          call res (neq, tn, y, s, rtem, ires)
          nre = nre + 1
          if (ires .gt. 1) go to 600
          do 215 i = 1,n
 215        rtem(i) = rtem(i) - savr(i)
          k1 = k
          do 230 i = j1,n,mwid
c           get jacobian elements in column i (block-column k1). -------
            y(i) = yh(i,1)
            r = dmax1(srur*dabs(y(i)),0.01d0/ewt(i))
            fac = -hl0/r
c           compute and load elements pa(*,j,k1). ----------------------
            iia = i - j
            ipa = 2 + (j-1)*mb + (k1-1)*mbsq
            do 221 j2 = 1,mb
 221          wm(ipa+j2) = rtem(iia+j2)*fac
            if (k1 .le. 1) go to 223
c           compute and load elements pb(*,j,k1-1). --------------------
            iib = iia - mb
            ipb = ipa + lblox - mbsq
            do 222 j2 = 1,mb
 222          wm(ipb+j2) = rtem(iib+j2)*fac
 223        continue
            if (k1 .ge. nb) go to 225
c           compute and load elements pc(*,j,k1+1). --------------------
            iic = iia + mb
            ipc = ipa + 2*lblox + mbsq
            do 224 j2 = 1,mb
 224          wm(ipc+j2) = rtem(iic+j2)*fac
 225        continue
            if (k1 .ne. 3) go to 227
c           compute and load elements pc(*,j,1). -----------------------
            ipc = ipa - 2*mbsq + 2*lblox
            do 226 j2 = 1,mb
 226          wm(ipc+j2) = rtem(j2)*fac
 227        continue
            if (k1 .ne. nb-2) go to 229
c           compute and load elements pb(*,j,nb). ----------------------
            iib = n - mb
            ipb = ipa + 2*mbsq + lblox
            do 228 j2 = 1,mb
 228          wm(ipb+j2) = rtem(iib+j2)*fac
 229      k1 = k1 + 3
 230      continue
 240    continue
 250  continue
c res call for first corrector iteration. ------------------------------
      ires = 1
      call res (neq, tn, y, s, savr, ires)
      nre = nre + 1
      if (ires .gt. 1) go to 600
c add matrix a. --------------------------------------------------------
 260  continue
      call adda (neq, tn, y, mb, nb, wm(3), wm(lpb), wm(lpc))
c do lu decomposition on p. --------------------------------------------
      call decbt (mb, nb, wm(3), wm(lpb), wm(lpc), iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c error return for ires = 2 or ires = 3 return from res. ---------------
 600  ierpj = ires
      return
c----------------------- end of subroutine pjibt -----------------------
      end
