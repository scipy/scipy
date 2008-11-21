c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: ssesrt
c
c\Description:
c  Sort the array X in the order specified by WHICH and optionally 
c  apply the permutation to the columns of the matrix A.
c
c\Usage:
c  call ssesrt
c     ( WHICH, APPLY, N, X, NA, A, LDA)
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> X is sorted into increasing order of magnitude.
c          'SM' -> X is sorted into decreasing order of magnitude.
c          'LA' -> X is sorted into increasing order of algebraic.
c          'SA' -> X is sorted into decreasing order of algebraic.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to A.
c          APPLY = .FALSE. -> do not apply the sorted order to A.
c
c  N       Integer.  (INPUT)
c          Dimension of the array X.
c
c  X      Real array of length N.  (INPUT/OUTPUT)
c          The array to be sorted.
c
c  NA      Integer.  (INPUT)
c          Number of rows of the matrix A.
c
c  A      Real array of length NA by N.  (INPUT/OUTPUT)
c         
c  LDA     Integer.  (INPUT)
c          Leading dimension of A.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Routines
c     sswap  Level 1 BLAS that swaps the contents of two vectors.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     12/15/93: Version ' 2.1'.
c               Adapted from the sort routine in LANSO and 
c               the ARPACK code ssortr
c
c\SCCS Information: @(#) 
c FILE: sesrt.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine ssesrt (which, apply, n, x, na, a, lda)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    lda, n, na
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Real
     &           x(0:n-1), a(lda, 0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Real
     &           temp
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   sswap
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c 
      if (which .eq. 'SA') then
c
c        X is sorted into decreasing order of algebraic.
c
   10    continue
         if (igap .eq. 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            if (x(j).lt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call sswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        X is sorted into decreasing order of magnitude.
c
   40    continue
         if (igap .eq. 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j.lt.0) go to 60
c
            if (abs(x(j)).lt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call sswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
c
      else if (which .eq. 'LA') then
c
c        X is sorted into increasing order of algebraic.
c
   70    continue
         if (igap .eq. 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c           
            if (x(j).gt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call sswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
c 
      else if (which .eq. 'LM') then
c
c        X is sorted into increasing order of magnitude.
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (abs(x(j)).gt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call sswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
      end if
c
 9000 continue
      return
c
c     %---------------%
c     | End of ssesrt |
c     %---------------%
c
      end
