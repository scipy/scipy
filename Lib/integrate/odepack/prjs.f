      subroutine prjs (neq,y,yh,nyh,ewt,ftem,savf,wk,iwk,f,jac)
clll. optimize
      external f,jac
      integer neq, nyh, iwk
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, imul, j, jj, jok, jmax, jmin, k, kmax, kmin, ng
      double precision y, yh, ewt, ftem, savf, wk
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision con0, conmin, ccmxj, psmall, rbig, seth
      double precision con, di, fac, hl0, pij, r, r0, rcon, rcont,
     1   srur, vnorm
      dimension neq(1), y(1), yh(nyh,*), ewt(1), ftem(1), savf(1),
     1   wk(*), iwk(*)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ con0, conmin, ccmxj, psmall, rbig, seth,
     1   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     4   nslj, ngp, nlu, nnz, nsp, nzl, nzu
c-----------------------------------------------------------------------
c prjs is called to compute and process the matrix
c p = i - h*el(1)*j , where j is an approximation to the jacobian.
c j is computed by columns, either by the user-supplied routine jac
c if miter = 1, or by finite differencing if miter = 2.
c if miter = 3, a diagonal approximation to j is used.
c if miter = 1 or 2, and if the existing value of the jacobian
c (as contained in p) is considered acceptable, then a new value of
c p is reconstructed from the old value.  in any case, when miter
c is 1 or 2, the p matrix is subjected to lu decomposition in cdrv.
c p and its lu decomposition are stored (separately) in wk.
c
c in addition to variables described previously, communication
c with prjs uses the following..
c y     = array containing predicted values on entry.
c ftem  = work array of length n (acor in stode).
c savf  = array containing f evaluated at predicted y.
c wk    = real work space for matrices.  on output it contains the
c         inverse diagonal matrix if miter = 3, and p and its sparse
c         lu decomposition if miter is 1 or 2.
c         storage of matrix elements starts at wk(3).
c         wk also contains the following matrix-related data..
c         wk(1) = sqrt(uround), used in numerical jacobian increments.
c         wk(2) = h*el0, saved for later use if miter = 3.
c iwk   = integer work space for matrix-related data, assumed to
c         be equivalenced to wk.  in addition, wk(iprsp) and iwk(ipisp)
c         are assumed to have identical locations.
c el0   = el(1) (input).
c ierpj = output error flag (in common).
c       = 0 if no error.
c       = 1  if zero pivot found in cdrv.
c       = 2  if a singular matrix arose with miter = 3.
c       = -1 if insufficient storage for cdrv (should not occur here).
c       = -2 if other error found in cdrv (should not occur here).
c jcur  = output flag = 1 to indicate that the jacobian matrix
c         (or approximation) is now current.
c this routine also uses other variables in common.
c-----------------------------------------------------------------------
      hl0 = h*el0
      con = -hl0
      if (miter .eq. 3) go to 300
c see whether j should be reevaluated (jok = 0) or not (jok = 1). ------
      jok = 1
      if (nst .eq. 0 .or. nst .ge. nslj+msbj) jok = 0
      if (icf .eq. 1 .and. dabs(rc - 1.0d0) .lt. ccmxj) jok = 0
      if (icf .eq. 2) jok = 0
      if (jok .eq. 1) go to 250
c
c miter = 1 or 2, and the jacobian is to be reevaluated. ---------------
 20   jcur = 1
      nje = nje + 1
      nslj = nst
      iplost = 0
      conmin = dabs(con)
      go to (100, 200), miter
c
c if miter = 1, call jac, multiply by scalar, and add identity. --------
 100  continue
      kmin = iwk(ipian)
      do 130 j = 1, n
        kmax = iwk(ipian+j) - 1
        do 110 i = 1,n
 110      ftem(i) = 0.0d0
        call jac (neq, tn, y, j, iwk(ipian), iwk(ipjan), ftem)
        do 120 k = kmin, kmax
          i = iwk(ibjan+k)
          wk(iba+k) = ftem(i)*con
          if (i .eq. j) wk(iba+k) = wk(iba+k) + 1.0d0
 120      continue
        kmin = kmax + 1
 130    continue
      go to 290
c
c if miter = 2, make ngp calls to f to approximate j and p. ------------
 200  continue
      fac = vnorm(n, savf, ewt)
      r0 = 1000.0d0 * dabs(h) * uround * dfloat(n) * fac
      if (r0 .eq. 0.0d0) r0 = 1.0d0
      srur = wk(1)
      jmin = iwk(ipigp)
      do 240 ng = 1,ngp
        jmax = iwk(ipigp+ng) - 1
        do 210 j = jmin,jmax
          jj = iwk(ibjgp+j)
          r = dmax1(srur*dabs(y(jj)),r0/ewt(jj))
 210      y(jj) = y(jj) + r
        call f (neq, tn, y, ftem)
        do 230 j = jmin,jmax
          jj = iwk(ibjgp+j)
          y(jj) = yh(jj,1)
          r = dmax1(srur*dabs(y(jj)),r0/ewt(jj))
          fac = -hl0/r
          kmin =iwk(ibian+jj)
          kmax =iwk(ibian+jj+1) - 1
          do 220 k = kmin,kmax
            i = iwk(ibjan+k)
            wk(iba+k) = (ftem(i) - savf(i))*fac
            if (i .eq. jj) wk(iba+k) = wk(iba+k) + 1.0d0
 220        continue
 230      continue
        jmin = jmax + 1
 240    continue
      nfe = nfe + ngp
      go to 290
c
c if jok = 1, reconstruct new p from old p. ----------------------------
 250  jcur = 0
      rcon = con/con0
      rcont = dabs(con)/conmin
      if (rcont .gt. rbig .and. iplost .eq. 1) go to 20
      kmin = iwk(ipian)
      do 275 j = 1,n
        kmax = iwk(ipian+j) - 1
        do 270 k = kmin,kmax
          i = iwk(ibjan+k)
          pij = wk(iba+k)
          if (i .ne. j) go to 260
          pij = pij - 1.0d0
          if (dabs(pij) .ge. psmall) go to 260
            iplost = 1
            conmin = dmin1(dabs(con0),conmin)
 260      pij = pij*rcon
          if (i .eq. j) pij = pij + 1.0d0
          wk(iba+k) = pij
 270      continue
        kmin = kmax + 1
 275    continue
c
c do numerical factorization of p matrix. ------------------------------
 290  nlu = nlu + 1
      con0 = con
      ierpj = 0
      do 295 i = 1,n
 295    ftem(i) = 0.0d0
      call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),
     1   wk(ipa),ftem,ftem,nsp,iwk(ipisp),wk(iprsp),iesp,2,iys)
      if (iys .eq. 0) return
      imul = (iys - 1)/n
      ierpj = -2
      if (imul .eq. 8) ierpj = 1
      if (imul .eq. 10) ierpj = -1
      return
c
c if miter = 3, construct a diagonal approximation to j and p. ---------
 300  continue
      jcur = 1
      nje = nje + 1
      wk(2) = hl0
      ierpj = 0
      r = el0*0.1d0
      do 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      call f (neq, tn, y, wk(3))
      nfe = nfe + 1
      do 320 i = 1,n
        r0 = h*savf(i) - yh(i,2)
        di = 0.1d0*r0 - h*(wk(i+2) - savf(i))
        wk(i+2) = 1.0d0
        if (dabs(r0) .lt. uround/ewt(i)) go to 320
        if (dabs(di) .eq. 0.0d0) go to 330
        wk(i+2) = 0.1d0*r0/di
 320    continue
      return
 330  ierpj = 2
      return
c----------------------- end of subroutine prjs ------------------------
      end
