      subroutine cntnzu (n, ia, ja, nzsut)
      integer n, ia, ja, nzsut
      dimension ia(1), ja(1)
c-----------------------------------------------------------------------
c this routine counts the number of nonzero elements in the strict
c upper triangle of the matrix m + m(transpose), where the sparsity
c structure of m is given by pointer arrays ia and ja.
c this is needed to compute the storage requirements for the
c sparse matrix reordering operation in odrv.
c-----------------------------------------------------------------------
      integer ii, jj, j, jmin, jmax, k, kmin, kmax, num
c
      num = 0
      do 50 ii = 1,n
        jmin = ia(ii)
        jmax = ia(ii+1) - 1
        if (jmin .gt. jmax) go to 50
        do 40 j = jmin,jmax
          if (ja(j) - ii) 10, 40, 30
 10       jj =ja(j)
          kmin = ia(jj)
          kmax = ia(jj+1) - 1
          if (kmin .gt. kmax) go to 30
          do 20 k = kmin,kmax
            if (ja(k) .eq. ii) go to 40
 20         continue
 30       num = num + 1
 40       continue
 50     continue
      nzsut = num
      return
c----------------------- end of subroutine cntnzu ----------------------
      end
