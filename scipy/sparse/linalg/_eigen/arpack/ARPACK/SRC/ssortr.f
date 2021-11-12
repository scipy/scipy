c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: ssortr
c
c\Description:
c  Sort the array X1 in the order specified by WHICH and optionally 
c  applies the permutation to the array X2.
c
c\Usage:
c  call ssortr
c     ( WHICH, APPLY, N, X1, X2 )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> X1 is sorted into increasing order of magnitude.
c          'SM' -> X1 is sorted into decreasing order of magnitude.
c          'LA' -> X1 is sorted into increasing order of algebraic.
c          'SA' -> X1 is sorted into decreasing order of algebraic.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to X2.
c          APPLY = .FALSE. -> do not apply the sorted order to X2.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  X1      Real array of length N.  (INPUT/OUTPUT)
c          The array to be sorted.
c
c  X2      Real array of length N.  (INPUT/OUTPUT)
c          Only referenced if APPLY = .TRUE.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     12/16/93: Version ' 2.1'.
c               Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#) 
c FILE: sortr.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine ssortr (which, apply, n, x1, x2)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Real
     &           x1(0:n-1), x2(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Real
     &           temp
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c 
      if (which .eq. 'SA') then
c
c        X1 is sorted into decreasing order of algebraic.
c
   10    continue
         if (igap .eq. 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            if (x1(j).lt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
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
c        X1 is sorted into decreasing order of magnitude.
c
   40    continue
         if (igap .eq. 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j.lt.0) go to 60
c
            if (abs(x1(j)).lt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
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
c        X1 is sorted into increasing order of algebraic.
c
   70    continue
         if (igap .eq. 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c           
            if (x1(j).gt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
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
c        X1 is sorted into increasing order of magnitude.
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (abs(x1(j)).gt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
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
c     | End of ssortr |
c     %---------------%
c
      end
