c\BeginDoc
c
c\Name: cngets
c
c\Description:
c  Given the eigenvalues of the upper Hessenberg matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: call this even in the case of user specified shifts in order
c  to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call cngets
c      ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> want the KEV eigenvalues of largest magnitude.
c          'SM' -> want the KEV eigenvalues of smallest magnitude.
c          'LR' -> want the KEV eigenvalues of largest REAL part.
c          'SR' -> want the KEV eigenvalues of smallest REAL part.
c          'LI' -> want the KEV eigenvalues of largest imaginary part.
c          'SI' -> want the KEV eigenvalues of smallest imaginary part.
c
c  KEV     Integer.  (INPUT)
c          The number of desired eigenvalues.
c
c  NP      Integer.  (INPUT)
c          The number of shifts to compute.
c
c  RITZ    Complex array of length KEV+NP.  (INPUT/OUTPUT)
c          On INPUT, RITZ contains the the eigenvalues of H.
c          On OUTPUT, RITZ are sorted so that the unwanted
c          eigenvalues are in the first NP locations and the wanted
c          portion is in the last KEV locations.  When exact shifts are
c          selected, the unwanted part corresponds to the shifts to
c          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
c          are further sorted so that the ones with largest Ritz values
c          are first.
c
c  BOUNDS  Complex array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex
c
c\Routines called:
c     csortc  ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     arscnd  ARPACK utility routine for timing.
c     cvout   ARPACK utility routine that prints vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: ngets.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. This routine does not keep complex conjugate pairs of
c        eigenvalues together.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine cngets ( ishift, which, kev, np, ritz, bounds)
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
      Complex
     &           bounds(kev+np), ritz(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex
     &           one, zero
      parameter (one = (1.0E+0, 0.0E+0), zero = (0.0E+0, 0.0E+0))
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    msglvl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   cvout,  csortc, arscnd
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
      msglvl = mcgets
c
      call csortc (which, .true., kev+np, ritz, bounds)
c
      if ( ishift .eq. 1 ) then
c
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first        |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when the shifts  |
c        | are applied in subroutine cnapps.                     |
c        | Be careful and use 'SM' since we want to sort BOUNDS! |
c        %-------------------------------------------------------%
c
         call csortc ( 'SM', .true., np, bounds, ritz )
c
      end if
c
      call arscnd (t1)
      tcgets = tcgets + (t1 - t0)
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, [kev], ndigit, '_ngets: KEV is')
         call ivout (logfil, 1, [np], ndigit, '_ngets: NP is')
         call cvout (logfil, kev+np, ritz, ndigit,
     &        '_ngets: Eigenvalues of current H matrix ')
         call cvout (logfil, kev+np, bounds, ndigit,
     &      '_ngets: Ritz estimates of the current KEV+NP Ritz values')
      end if
c
      return
c
c     %---------------%
c     | End of cngets |
c     %---------------%
c
      end
