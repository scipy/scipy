      subroutine adjlr (n, isp, ldif)
      integer n, isp, ldif
      dimension isp(1)
c-----------------------------------------------------------------------
c this routine computes an adjustment, ldif, to the required
c integer storage space in iwk (sparse matrix work space).
c it is called only if the word length ratio is lrat = 1.
c this is to account for the possibility that the symbolic lu phase
c may require more storage than the numerical lu and solution phases.
c-----------------------------------------------------------------------
      integer ip, jlmax, jumax, lnfc, lsfc, nzlu
c
      ip = 2*n + 1
c get jlmax = ijl(n) and jumax = iju(n) (sizes of jl and ju). ----------
      jlmax = isp(ip)
      jumax = isp(ip+ip)
c nzlu = (size of l) + (size of u) = (il(n+1)-il(1)) + (iu(n+1)-iu(1)).
      nzlu = isp(n+1) - isp(1) + isp(ip+n+1) - isp(ip+1)
      lsfc = 12*n + 3 + 2*max0(jlmax,jumax)
      lnfc = 9*n + 2 + jlmax + jumax + nzlu
      ldif = max0(0, lsfc - lnfc)
      return
c----------------------- end of subroutine adjlr -----------------------
      end
