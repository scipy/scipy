c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dngets
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
c  call dngets
c     ( ISHIFT, WHICH, KEV, NP, RITZR, RITZI, BOUNDS, SHIFTR, SHIFTI )
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
c          'LR' -> want the KEV eigenvalues of largest real part.
c          'SR' -> want the KEV eigenvalues of smallest real part.
c          'LI' -> want the KEV eigenvalues of largest imaginary part.
c          'SI' -> want the KEV eigenvalues of smallest imaginary part.
c
c  KEV      Integer.  (INPUT/OUTPUT)
c           INPUT: KEV+NP is the size of the matrix H.
c           OUTPUT: Possibly increases KEV by one to keep complex conjugate
c           pairs together.
c
c  NP       Integer.  (INPUT/OUTPUT)
c           Number of implicit shifts to be computed.
c           OUTPUT: Possibly decreases NP by one to keep complex conjugate
c           pairs together.
c
c  RITZR,  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c  RITZI   On INPUT, RITZR and RITZI contain the real and imaginary 
c          parts of the eigenvalues of H.
c          On OUTPUT, RITZR and RITZI are sorted so that the unwanted
c          eigenvalues are in the first NP locations and the wanted
c          portion is in the last KEV locations.  When exact shifts are 
c          selected, the unwanted part corresponds to the shifts to 
c          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
c          are further sorted so that the ones with largest Ritz values
c          are first.
c
c  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  SHIFTR, SHIFTI  *** USE deprecated as of version 2.1. ***
c  
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
c     dsortc  ARPACK sorting routine.
c     dcopy   Level 1 BLAS that copies one vector to another .
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
c
c\SCCS Information: @(#) 
c FILE: ngets.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. xxxx
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dngets ( ishift, which, kev, np, ritzr, ritzi, bounds,
     &                    shiftr, shifti )
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
      Double precision
     &           bounds(kev+np), ritzr(kev+np), ritzi(kev+np), 
     &           shiftr(1), shifti(1)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0, zero = 0.0)
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
      external   dcopy, dsortc, second
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs
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
      call second (t0)
      msglvl = mngets
c 
c     %----------------------------------------------------%
c     | LM, SM, LR, SR, LI, SI case.                       |
c     | Sort the eigenvalues of H into the desired order   |
c     | and apply the resulting order to BOUNDS.           |
c     | The eigenvalues are sorted so that the wanted part |
c     | are always in the last KEV locations.              |
c     | We first do a pre-processing sort in order to keep |
c     | complex conjugate pairs together                   |
c     %----------------------------------------------------%
c
      if (which .eq. 'LM') then
         call dsortc ('LR', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SM') then
         call dsortc ('SR', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'LR') then
         call dsortc ('LM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SR') then
         call dsortc ('SM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'LI') then
         call dsortc ('LM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SI') then
         call dsortc ('SM', .true., kev+np, ritzr, ritzi, bounds)
      end if
c      
      call dsortc (which, .true., kev+np, ritzr, ritzi, bounds)
c     
c     %-------------------------------------------------------%
c     | Increase KEV by one if the ( ritzr(np),ritzi(np) )    |
c     | = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero |
c     | Accordingly decrease NP by one. In other words keep   |
c     | complex conjugate pairs together.                     |
c     %-------------------------------------------------------%
c     
      if (       ( ritzr(np+1) - ritzr(np) ) .eq. zero
     &     .and. ( ritzi(np+1) + ritzi(np) ) .eq. zero ) then
         np = np - 1
         kev = kev + 1
      end if
c
      if ( ishift .eq. 1 ) then
c     
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first        |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when they shifts |
c        | are applied in subroutine dnapps.                     |
c        | Be careful and use 'SR' since we want to sort BOUNDS! |
c        %-------------------------------------------------------%
c     
         call dsortc ( 'SR', .true., np, bounds, ritzr, ritzi )
      end if
c     
      call second (t1)
      tngets = tngets + (t1 - t0)
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, kev, ndigit, '_ngets: KEV is')
         call ivout (logfil, 1, np, ndigit, '_ngets: NP is')
         call dvout (logfil, kev+np, ritzr, ndigit,
     &        '_ngets: Eigenvalues of current H matrix -- real part')
         call dvout (logfil, kev+np, ritzi, ndigit,
     &        '_ngets: Eigenvalues of current H matrix -- imag part')
         call dvout (logfil, kev+np, bounds, ndigit, 
     &      '_ngets: Ritz estimates of the current KEV+NP Ritz values')
      end if
c     
      return
c     
c     %---------------%
c     | End of dngets |
c     %---------------%
c     
      end
