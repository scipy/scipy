      subroutine prep (neq, y, yh, savf, ewt, ftem, ia, ja,
     1                     wk, iwk, ipper, f, jac)
clll. optimize
      external f,jac
      integer neq, ia, ja, iwk, ipper
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, jfound, k,
     1   knew, kmax, kmin, ldif, lenigp, liwk, maxg, np1, nzsut
      double precision y, yh, savf, ewt, ftem, wk
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision con0, conmin, ccmxj, psmall, rbig, seth
      double precision dq, dyj, erwt, fac, yj
      dimension neq(1), y(1), yh(1), savf(1), ewt(1), ftem(1),
     1   ia(1), ja(1), wk(1), iwk(1)
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
c this routine performs preprocessing related to the sparse linear
c systems that must be solved if miter = 1 or 2.
c the operations that are performed here are..
c  * compute sparseness structure of jacobian according to moss,
c  * compute grouping of column indices (miter = 2),
c  * compute a new ordering of rows and columns of the matrix,
c  * reorder ja corresponding to the new ordering,
c  * perform a symbolic lu factorization of the matrix, and
c  * set pointers for segments of the iwk/wk array.
c in addition to variables described previously, prep uses the
c following for communication..
c yh     = the history array.  only the first column, containing the
c          current y vector, is used.  used only if moss .ne. 0.
c savf   = a work array of length neq, used only if moss .ne. 0.
c ewt    = array of length neq containing (inverted) error weights.
c          used only if moss = 2 or if istate = moss = 1.
c ftem   = a work array of length neq, identical to acor in the driver,
c          used only if moss = 2.
c wk     = a real work array of length lenwk, identical to wm in
c          the driver.
c iwk    = integer work array, assumed to occupy the same space as wk.
c lenwk  = the length of the work arrays wk and iwk.
c istatc = a copy of the driver input argument istate (= 1 on the
c          first call, = 3 on a continuation call).
c iys    = flag value from odrv or cdrv.
c ipper  = output error flag with the following values and meanings..
c          0  no error.
c         -1  insufficient storage for internal structure pointers.
c         -2  insufficient storage for jgroup.
c         -3  insufficient storage for odrv.
c         -4  other error flag from odrv (should never occur).
c         -5  insufficient storage for cdrv.
c         -6  other error flag from cdrv.
c-----------------------------------------------------------------------
      ibian = lrat*2
      ipian = ibian + 1
      np1 = n + 1
      ipjan = ipian + np1
      ibjan = ipjan - 1
      liwk = lenwk*lrat
      if (ipjan+n-1 .gt. liwk) go to 210
      if (moss .eq. 0) go to 30
c
      if (istatc .eq. 3) go to 20
c istate = 1 and moss .ne. 0.  perturb y for structure determination. --
      do 10 i = 1,n
        erwt = 1.0d0/ewt(i)
        fac = 1.0d0 + 1.0d0/(dfloat(i)+1.0d0)
        y(i) = y(i) + fac*dsign(erwt,y(i))
 10     continue
      go to (70, 100), moss
c
 20   continue
c istate = 3 and moss .ne. 0.  load y from yh(*,1). --------------------
      do 25 i = 1,n
 25     y(i) = yh(i)
      go to (70, 100), moss
c
c moss = 0.  process user-s ia,ja.  add diagonal entries if necessary. -
 30   knew = ipjan
      kmin = ia(1)
      iwk(ipian) = 1
      do 60 j = 1,n
        jfound = 0
        kmax = ia(j+1) - 1
        if (kmin .gt. kmax) go to 45
        do 40 k = kmin,kmax
          i = ja(k)
          if (i .eq. j) jfound = 1
          if (knew .gt. liwk) go to 210
          iwk(knew) = i
          knew = knew + 1
 40       continue
        if (jfound .eq. 1) go to 50
 45     if (knew .gt. liwk) go to 210
        iwk(knew) = j
        knew = knew + 1
 50     iwk(ipian+j) = knew + 1 - ipjan
        kmin = kmax + 1
 60     continue
      go to 140
c
c moss = 1.  compute structure from user-supplied jacobian routine jac.
 70   continue
c a dummy call to f allows user to create temporaries for use in jac. --
      call f (neq, tn, y, savf)
      k = ipjan
      iwk(ipian) = 1
      do 90 j = 1,n
        if (k .gt. liwk) go to 210
        iwk(k) = j
        k = k + 1
        do 75 i = 1,n
 75       savf(i) = 0.0d0
        call jac (neq, tn, y, j, iwk(ipian), iwk(ipjan), savf)
        do 80 i = 1,n
          if (dabs(savf(i)) .le. seth) go to 80
          if (i .eq. j) go to 80
          if (k .gt. liwk) go to 210
          iwk(k) = i
          k = k + 1
 80       continue
        iwk(ipian+j) = k + 1 - ipjan
 90     continue
      go to 140
c
c moss = 2.  compute structure from results of n + 1 calls to f. -------
 100  k = ipjan
      iwk(ipian) = 1
      call f (neq, tn, y, savf)
      do 120 j = 1,n
        if (k .gt. liwk) go to 210
        iwk(k) = j
        k = k + 1
        yj = y(j)
        erwt = 1.0d0/ewt(j)
        dyj = dsign(erwt,yj)
        y(j) = yj + dyj
        call f (neq, tn, y, ftem)
        y(j) = yj
        do 110 i = 1,n
          dq = (ftem(i) - savf(i))/dyj
          if (dabs(dq) .le. seth) go to 110
          if (i .eq. j) go to 110
          if (k .gt. liwk) go to 210
          iwk(k) = i
          k = k + 1
 110      continue
        iwk(ipian+j) = k + 1 - ipjan
 120    continue
c
 140  continue
      if (moss .eq. 0 .or. istatc .ne. 1) go to 150
c if istate = 1 and moss .ne. 0, restore y from yh. --------------------
      do 145 i = 1,n
 145    y(i) = yh(i)
 150  nnz = iwk(ipian+n) - 1
      lenigp = 0
      ipigp = ipjan + nnz
      if (miter .ne. 2) go to 160
c
c compute grouping of column indices (miter = 2). ----------------------
      maxg = np1
      ipjgp = ipjan + nnz
      ibjgp = ipjgp - 1
      ipigp = ipjgp + n
      iptt1 = ipigp + np1
      iptt2 = iptt1 + n
      lreq = iptt2 + n - 1
      if (lreq .gt. liwk) go to 220
      call jgroup (n, iwk(ipian), iwk(ipjan), maxg, ngp, iwk(ipigp),
     1   iwk(ipjgp), iwk(iptt1), iwk(iptt2), ier)
      if (ier .ne. 0) go to 220
      lenigp = ngp + 1
c
c compute new ordering of rows/columns of jacobian. --------------------
 160  ipr = ipigp + lenigp
      ipc = ipr
      ipic = ipc + n
      ipisp = ipic + n
      iprsp = (ipisp - 2)/lrat + 2
      iesp = lenwk + 1 - iprsp
      if (iesp .lt. 0) go to 230
      ibr = ipr - 1
      do 170 i = 1,n
 170    iwk(ibr+i) = i
      nsp = liwk + 1 - ipisp
      call odrv (n, iwk(ipian), iwk(ipjan), wk, iwk(ipr), iwk(ipic),
     1   nsp, iwk(ipisp), 1, iys)
      if (iys .eq. 11*n+1) go to 240
      if (iys .ne. 0) go to 230
c
c reorder jan and do symbolic lu factorization of matrix. --------------
      ipa = lenwk + 1 - nnz
      nsp = ipa - iprsp
      lreq = max0(12*n/lrat, 6*n/lrat+2*n+nnz) + 3
      lreq = lreq + iprsp - 1 + nnz
      if (lreq .gt. lenwk) go to 250
      iba = ipa - 1
      do 180 i = 1,nnz
 180    wk(iba+i) = 0.0d0
      ipisp = lrat*(iprsp - 1) + 1
      call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),
     1   wk(ipa),wk(ipa),wk(ipa),nsp,iwk(ipisp),wk(iprsp),iesp,5,iys)
      lreq = lenwk - iesp
      if (iys .eq. 10*n+1) go to 250
      if (iys .ne. 0) go to 260
      ipil = ipisp
      ipiu = ipil + 2*n + 1
      nzu = iwk(ipil+n) - iwk(ipil)
      nzl = iwk(ipiu+n) - iwk(ipiu)
      if (lrat .gt. 1) go to 190
      call adjlr (n, iwk(ipisp), ldif)
      lreq = lreq + ldif
 190  continue
      if (lrat .eq. 2 .and. nnz .eq. n) lreq = lreq + 1
      nsp = nsp + lreq - lenwk
      ipa = lreq + 1 - nnz
      iba = ipa - 1
      ipper = 0
      return
c
 210  ipper = -1
      lreq = 2 + (2*n + 1)/lrat
      lreq = max0(lenwk+1,lreq)
      return
c
 220  ipper = -2
      lreq = (lreq - 1)/lrat + 1
      return
c
 230  ipper = -3
      call cntnzu (n, iwk(ipian), iwk(ipjan), nzsut)
      lreq = lenwk - iesp + (3*n + 4*nzsut - 1)/lrat + 1
      return
c
 240  ipper = -4
      return
c
 250  ipper = -5
      return
c
 260  ipper = -6
      lreq = lenwk
      return
c----------------------- end of subroutine prep ------------------------
      end
