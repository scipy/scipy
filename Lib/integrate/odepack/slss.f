      subroutine slss (wk, iwk, x, tem)
clll. optimize
      integer iwk
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i
      double precision wk, x, tem
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision rlss
      double precision di, hl0, phl0, r
      dimension wk(*), iwk(*), x(1), tem(1)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ rlss(6),
     1   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     4   nslj, ngp, nlu, nnz, nsp, nzl, nzu
c-----------------------------------------------------------------------
c this routine manages the solution of the linear system arising from
c a chord iteration.  it is called if miter .ne. 0.
c if miter is 1 or 2, it calls cdrv to accomplish this.
c if miter = 3 it updates the coefficient h*el0 in the diagonal
c matrix, and then computes the solution.
c communication with slss uses the following variables..
c wk    = real work space containing the inverse diagonal matrix if
c         miter = 3 and the lu decomposition of the matrix otherwise.
c         storage of matrix elements starts at wk(3).
c         wk also contains the following matrix-related data..
c         wk(1) = sqrt(uround) (not used here),
c         wk(2) = hl0, the previous value of h*el0, used if miter = 3.
c iwk   = integer work space for matrix-related data, assumed to
c         be equivalenced to wk.  in addition, wk(iprsp) and iwk(ipisp)
c         are assumed to have identical locations.
c x     = the right-hand side vector on input, and the solution vector
c         on output, of length n.
c tem   = vector of work space of length n, not used in this version.
c iersl = output flag (in common).
c         iersl = 0  if no trouble occurred.
c         iersl = -1 if cdrv returned an error flag (miter = 1 or 2).
c                    this should never occur and is considered fatal.
c         iersl = 1  if a singular matrix arose with miter = 3.
c this routine also uses other variables in common.
c-----------------------------------------------------------------------
      iersl = 0
      go to (100, 100, 300), miter
 100  call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),
     1   wk(ipa),x,x,nsp,iwk(ipisp),wk(iprsp),iesp,4,iersl)
      if (iersl .ne. 0) iersl = -1
      return
c
 300  phl0 = wk(2)
      hl0 = h*el0
      wk(2) = hl0
      if (hl0 .eq. phl0) go to 330
      r = hl0/phl0
      do 320 i = 1,n
        di = 1.0d0 - r*(1.0d0 - 1.0d0/wk(i+2))
        if (dabs(di) .eq. 0.0d0) go to 390
 320    wk(i+2) = 1.0d0/di
 330  do 340 i = 1,n
 340    x(i) = wk(i+2)*x(i)
      return
 390  iersl = 1
      return
c
c----------------------- end of subroutine slss ------------------------
      end
