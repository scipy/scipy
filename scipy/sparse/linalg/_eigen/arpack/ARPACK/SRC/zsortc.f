c\BeginDoc
c
c\Name: zsortc
c
c\Description:
c  Sorts the Complex*16 array in X into the order
c  specified by WHICH and optionally applies the permutation to the
c  Double precision  array Y.
c
c\Usage:
c  call zsortc
c     ( WHICH, APPLY, N, X, Y )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> sort X into increasing order of magnitude.
c          'SM' -> sort X into decreasing order of magnitude.
c          'LR' -> sort X with real(X) in increasing algebraic order
c          'SR' -> sort X with real(X) in decreasing algebraic order
c          'LI' -> sort X with imag(X) in increasing algebraic order
c          'SI' -> sort X with imag(X) in decreasing algebraic order
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to array Y.
c          APPLY = .FALSE. -> do not apply the sorted order to array Y.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  X       Complex*16 array of length N.  (INPUT/OUTPUT)
c          This is the array to be sorted.
c
c  Y       Complex*16 array of length N.  (INPUT/OUTPUT)
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Routines called:
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c     Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#)
c FILE: sortc.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine zsortc (which, apply, n, x, y)
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
      Complex*16
     &           x(0:n-1), y(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Complex*16
     &           temp
      Double precision
     &           temp1, temp2
c
c     %--------------------%
c     | External functions |
c     %--------------------%
c
      Double precision
     &           dlapy2
c
c     %--------------------%
c     | Intrinsic Functions |
c     %--------------------%
       Intrinsic
     &           dble, aimag
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c
      if (which .eq. 'LM') then
c
c        %--------------------------------------------%
c        | Sort X into increasing order of magnitude. |
c        %--------------------------------------------%
c
   10    continue
         if (igap .eq. 0) go to 9000
c
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            temp1 = dlapy2(dble(x(j)),aimag(x(j)))
            temp2 = dlapy2(dble(x(j+igap)),aimag(x(j+igap)))
c
            if (temp1.gt.temp2) then
                temp = x(j)
                x(j) = x(j+igap)
                x(j+igap) = temp
c
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                go to 30
            end if
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        %--------------------------------------------%
c        | Sort X into decreasing order of magnitude. |
c        %--------------------------------------------%
c
   40    continue
         if (igap .eq. 0) go to 9000
c
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j .lt. 0) go to 60
c
            temp1 = dlapy2(dble(x(j)),aimag(x(j)))
            temp2 = dlapy2(dble(x(j+igap)),aimag(x(j+igap)))
c
            if (temp1.lt.temp2) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
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
      else if (which .eq. 'LR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into increasing order of algebraic. |
c        %------------------------------------------------%
c
   70    continue
         if (igap .eq. 0) go to 9000
c
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c
            if (dble(x(j)).gt.dble(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
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
      else if (which .eq. 'SR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into decreasing order of algebraic. |
c        %------------------------------------------------%
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (dble(x(j)).lt.dble(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
c
      else if (which .eq. 'LI') then
c
c        %--------------------------------------------%
c        | Sort XIMAG into increasing algebraic order |
c        %--------------------------------------------%
c
  130    continue
         if (igap .eq. 0) go to 9000
         do 150 i = igap, n-1
            j = i-igap
  140       continue
c
            if (j.lt.0) go to 150
c
            if (aimag(x(j)).gt.aimag(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 150
            endif
            j = j-igap
            go to 140
  150    continue
         igap = igap / 2
         go to 130
c
      else if (which .eq. 'SI') then
c
c        %---------------------------------------------%
c        | Sort XIMAG into decreasing algebraic order  |
c        %---------------------------------------------%
c
  160    continue
         if (igap .eq. 0) go to 9000
         do 180 i = igap, n-1
            j = i-igap
  170       continue
c
            if (j.lt.0) go to 180
c
            if (aimag(x(j)).lt.aimag(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 180
            endif
            j = j-igap
            go to 170
  180    continue
         igap = igap / 2
         go to 160
      end if
c
 9000 continue
      return
c
c     %---------------%
c     | End of zsortc |
c     %---------------%
c
      end
