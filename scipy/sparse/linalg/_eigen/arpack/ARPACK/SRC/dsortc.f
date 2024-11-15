c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsortc
c
c\Description:
c  Sorts the complex array in XREAL and XIMAG into the order
c  specified by WHICH and optionally applies the permutation to the
c  real array Y. It is assumed that if an element of XIMAG is
c  nonzero, then its negative is also an element. In other words,
c  both members of a complex conjugate pair are to be sorted and the
c  pairs are kept adjacent to each other.
c
c\Usage:
c  call dsortc
c     ( WHICH, APPLY, N, XREAL, XIMAG, Y )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> sort XREAL,XIMAG into increasing order of magnitude.
c          'SM' -> sort XREAL,XIMAG into decreasing order of magnitude.
c          'LR' -> sort XREAL into increasing order of algebraic.
c          'SR' -> sort XREAL into decreasing order of algebraic.
c          'LI' -> sort XIMAG into increasing order of magnitude.
c          'SI' -> sort XIMAG into decreasing order of magnitude.
c          NOTE: If an element of XIMAG is non-zero, then its negative
c                is also an element.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to array Y.
c          APPLY = .FALSE. -> do not apply the sorted order to array Y.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  XREAL,  Double precision array of length N.  (INPUT/OUTPUT)
c  XIMAG   Real and imaginary part of the array to be sorted.
c
c  Y       Double precision array of length N.  (INPUT/OUTPUT)
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
c     xx/xx/92: Version ' 2.1'
c               Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#)
c FILE: sortc.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsortc (which, apply, n, xreal, ximag, y)
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
      Double precision
     &           xreal(0:n-1), ximag(0:n-1), y(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Double precision
     &           temp, temp1, temp2
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2
      external   dlapy2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c
      if (which .eq. 'LM') then
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into increasing order of magnitude. |
c        %------------------------------------------------------%
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
            temp1 = dlapy2(xreal(j),ximag(j))
            temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
c
            if (temp1.gt.temp2) then
                temp = xreal(j)
                xreal(j) = xreal(j+igap)
                xreal(j+igap) = temp
c
                temp = ximag(j)
                ximag(j) = ximag(j+igap)
                ximag(j+igap) = temp
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
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------------%
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
            temp1 = dlapy2(xreal(j),ximag(j))
            temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
c
            if (temp1.lt.temp2) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
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
            if (xreal(j).gt.xreal(j+igap)) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
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
            if (xreal(j).lt.xreal(j+igap)) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
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
c        %------------------------------------------------%
c        | Sort XIMAG into increasing order of magnitude. |
c        %------------------------------------------------%
c
  130    continue
         if (igap .eq. 0) go to 9000
         do 150 i = igap, n-1
            j = i-igap
  140       continue
c
            if (j.lt.0) go to 150
c
            if (abs(ximag(j)).gt.abs(ximag(j+igap))) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
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
c        %------------------------------------------------%
c        | Sort XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------%
c
  160    continue
         if (igap .eq. 0) go to 9000
         do 180 i = igap, n-1
            j = i-igap
  170       continue
c
            if (j.lt.0) go to 180
c
            if (abs(ximag(j)).lt.abs(ximag(j+igap))) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
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
c     | End of dsortc |
c     %---------------%
c
      end
