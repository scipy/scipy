      subroutine iprep (neq, y, rwork, ia, ja, ipflag, f, jac)
clll. optimize
      external f, jac
      integer neq, ia, ja, ipflag
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, imax, lewtn, lyhd, lyhn
      double precision y, rwork
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision rlss
      dimension neq(1), y(1), rwork(1), ia(1), ja(1)
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ rlss(6),
     1   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     4   nslj, ngp, nlu, nnz, nsp, nzl, nzu
c-----------------------------------------------------------------------
c this routine serves as an interface between the driver and
c subroutine prep.  it is called only if miter is 1 or 2.
c tasks performed here are..
c  * call prep,
c  * reset the required wm segment length lenwk,
c  * move yh back to its final location (following wm in rwork),
c  * reset pointers for yh, savf, ewt, and acor, and
c  * move ewt to its new position if istate = 1.
c ipflag is an output error indication flag.  ipflag = 0 if there was
c no trouble, and ipflag is the value of the prep error flag ipper
c if there was trouble in subroutine prep.
c-----------------------------------------------------------------------
      ipflag = 0
c call prep to do matrix preprocessing operations. ---------------------
      call prep (neq, y, rwork(lyh), rwork(lsavf), rwork(lewt),
     1   rwork(lacor), ia, ja, rwork(lwm), rwork(lwm), ipflag, f, jac)
      lenwk = max0(lreq,lwmin)
      if (ipflag .lt. 0) return
c if prep was successful, move yh to end of required space for wm. -----
      lyhn = lwm + lenwk
      if (lyhn .gt. lyh) return
      lyhd = lyh - lyhn
      if (lyhd .eq. 0) go to 20
      imax = lyhn - 1 + lenyhm
      do 10 i = lyhn,imax
 10     rwork(i) = rwork(i+lyhd)
      lyh = lyhn
c reset pointers for savf, ewt, and acor. ------------------------------
 20   lsavf = lyh + lenyh
      lewtn = lsavf + n
      lacor = lewtn + n
      if (istatc .eq. 3) go to 40
c if istate = 1, move ewt (left) to its new position. ------------------
      if (lewtn .gt. lewt) return
      do 30 i = 1,n
 30     rwork(i+lewtn-1) = rwork(i+lewt-1)
 40   lewt = lewtn
      return
c----------------------- end of subroutine iprep -----------------------
      end
