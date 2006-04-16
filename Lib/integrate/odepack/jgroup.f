      subroutine jgroup (n,ia,ja,maxg,ngrp,igp,jgp,incl,jdone,ier)
clll. optimize
      integer n, ia, ja, maxg, ngrp, igp, jgp, incl, jdone, ier
      dimension ia(1), ja(1), igp(1), jgp(n), incl(n), jdone(n)
c-----------------------------------------------------------------------
c this subroutine constructs groupings of the column indices of
c the jacobian matrix, used in the numerical evaluation of the
c jacobian by finite differences.
c
c input..
c n      = the order of the matrix.
c ia,ja  = sparse structure descriptors of the matrix by rows.
c maxg   = length of available storate in the igp array.
c
c output..
c ngrp   = number of groups.
c jgp    = array of length n containing the column indices by groups.
c igp    = pointer array of length ngrp + 1 to the locations in jgp
c          of the beginning of each group.
c ier    = error indicator.  ier = 0 if no error occurred, or 1 if
c          maxg was insufficient.
c
c incl and jdone are working arrays of length n.
c-----------------------------------------------------------------------
      integer i, j, k, kmin, kmax, ncol, ng
c
      ier = 0
      do 10 j = 1,n
 10     jdone(j) = 0
      ncol = 1
      do 60 ng = 1,maxg
        igp(ng) = ncol
        do 20 i = 1,n
 20       incl(i) = 0
        do 50 j = 1,n
c reject column j if it is already in a group.--------------------------
          if (jdone(j) .eq. 1) go to 50
          kmin = ia(j)
          kmax = ia(j+1) - 1
          do 30 k = kmin,kmax
c reject column j if it overlaps any column already in this group.------
            i = ja(k)
            if (incl(i) .eq. 1) go to 50
 30         continue
c accept column j into group ng.----------------------------------------
          jgp(ncol) = j
          ncol = ncol + 1
          jdone(j) = 1
          do 40 k = kmin,kmax
            i = ja(k)
 40         incl(i) = 1
 50       continue
c stop if this group is empty (grouping is complete).-------------------
        if (ncol .eq. igp(ng)) go to 70
 60     continue
c error return if not all columns were chosen (maxg too small).---------
      if (ncol .le. n) go to 80
      ng = maxg
 70   ngrp = ng - 1
      return
 80   ier = 1
      return
c----------------------- end of subroutine jgroup ----------------------
      end
