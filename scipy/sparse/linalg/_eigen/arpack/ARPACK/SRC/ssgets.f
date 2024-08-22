c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: ssgets
c
c\Description:
c  Given the eigenvalues of the symmetric tridiagonal matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: This is called even in the case of user specified shifts in
c  order to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call ssgets
c     ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> KEV eigenvalues of largest magnitude are retained.
c          'SM' -> KEV eigenvalues of smallest magnitude are retained.
c          'LA' -> KEV eigenvalues of largest value are retained.
c          'SA' -> KEV eigenvalues of smallest value are retained.
c          'BE' -> KEV eigenvalues, half from each end of the spectrum.
c                  If KEV is odd, compute one more from the high end.
c
c  KEV      Integer.  (INPUT)
c          KEV+NP is the size of the matrix H.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be computed.
c
c  RITZ    Real array of length KEV+NP.  (INPUT/OUTPUT)
c          On INPUT, RITZ contains the eigenvalues of H.
c          On OUTPUT, RITZ are sorted so that the unwanted eigenvalues
c          are in the first NP locations and the wanted part is in
c          the last KEV locations.  When exact shifts are selected, the
c          unwanted part corresponds to the shifts to be applied.
c
c  BOUNDS  Real array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  SHIFTS  Real array of length NP.  (INPUT/OUTPUT)
c          On INPUT:  contains the user specified shifts if ISHIFT = 0.
c          On OUTPUT: contains the shifts sorted into decreasing order
c          of magnitude with respect to the Ritz estimates contained in
c          BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     ssortr  ARPACK utility sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     arscnd  ARPACK utility routine for timing.
c     svout   ARPACK utility routine that prints vectors.
c     scopy   Level 1 BLAS that copies one vector to another.
c     sswap   Level 1 BLAS that swaps the contents of two vectors.
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
c     xx/xx/93: Version ' 2.1'
c
c\SCCS Information: @(#)
c FILE: sgets.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine ssgets ( ishift, which, kev, np, ritz, bounds, shifts )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      integer    ishift, kev, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Real
     &           bounds(kev+np), ritz(kev+np), shifts(np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &           one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    kevd2, msglvl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   sswap, scopy, ssortr, arscnd
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    max, min
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call arscnd (t0)
      msglvl = msgets
c
      if (which .eq. 'BE') then
c
c        %-----------------------------------------------------%
c        | Both ends of the spectrum are requested.            |
c        | Sort the eigenvalues into algebraically increasing  |
c        | order first then swap high end of the spectrum next |
c        | to low end in appropriate locations.                |
c        | NOTE: when np < floor(kev/2) be careful not to swap |
c        | overlapping locations.                              |
c        %-----------------------------------------------------%
c
         call ssortr ('LA', .true., kev+np, ritz, bounds)
         kevd2 = kev / 2
         if ( kev .gt. 1 ) then
            call sswap ( min(kevd2,np), ritz, 1,
     &                   ritz( max(kevd2,np)+1 ), 1)
            call sswap ( min(kevd2,np), bounds, 1,
     &                   bounds( max(kevd2,np)+1 ), 1)
         end if
c
      else
c
c        %----------------------------------------------------%
c        | LM, SM, LA, SA case.                               |
c        | Sort the eigenvalues of H into the desired order   |
c        | and apply the resulting order to BOUNDS.           |
c        | The eigenvalues are sorted so that the wanted part |
c        | are always in the last KEV locations.               |
c        %----------------------------------------------------%
c
         call ssortr (which, .true., kev+np, ritz, bounds)
      end if
c
      if (ishift .eq. 1 .and. np .gt. 0) then
c
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first.       |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when the shifts  |
c        | are applied in subroutine ssapps.                     |
c        %-------------------------------------------------------%
c
         call ssortr ('SM', .true., np, bounds, ritz)
         call scopy (np, ritz, 1, shifts, 1)
      end if
c
      call arscnd (t1)
      tsgets = tsgets + (t1 - t0)
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, [kev], ndigit, '_sgets: KEV is')
         call ivout (logfil, 1, [np], ndigit, '_sgets: NP is')
         call svout (logfil, kev+np, ritz, ndigit,
     &        '_sgets: Eigenvalues of current H matrix')
         call svout (logfil, kev+np, bounds, ndigit,
     &        '_sgets: Associated Ritz estimates')
      end if
c
      return
c
c     %---------------%
c     | End of ssgets |
c     %---------------%
c
      end
