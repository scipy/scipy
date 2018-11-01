c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: ssconv
c
c\Description: 
c  Convergence testing for the symmetric Arnoldi eigenvalue routine.
c
c\Usage:
c  call ssconv
c     ( N, RITZ, BOUNDS, TOL, NCONV )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Number of Ritz values to check for convergence.
c
c  RITZ    Real array of length N.  (INPUT)
c          The Ritz values to be checked for convergence.
c
c  BOUNDS  Real array of length N.  (INPUT)
c          Ritz estimates associated with the Ritz values in RITZ.
c
c  TOL     Real scalar.  (INPUT)
c          Desired relative accuracy for a Ritz value to be considered
c          "converged".
c
c  NCONV   Integer scalar.  (OUTPUT)
c          Number of "converged" Ritz values.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Routines called:
c     arscnd  ARPACK utility routine for timing.
c     slamch  LAPACK routine that determines machine constants. 
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
c FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
c
c\Remarks
c     1. Starting with version 2.4, this routine no longer uses the
c        Parlett strategy using the gap conditions. 
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine ssconv (n, ritz, bounds, tol, nconv)
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
      integer    n, nconv
      Real
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Real
     &           ritz(n), bounds(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i
      Real
     &           temp, eps23
c
c     %-------------------%
c     | External routines |
c     %-------------------%
c
      Real
     &           slamch
      external   slamch

c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      call arscnd (t0)
c
      eps23 = slamch('Epsilon-Machine') 
      eps23 = eps23**(2.0E+0 / 3.0E+0)
c
      nconv  = 0
      do 10 i = 1, n
c
c        %-----------------------------------------------------%
c        | The i-th Ritz value is considered "converged"       |
c        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |
c        %-----------------------------------------------------%
c
         temp = max( eps23, abs(ritz(i)) )
         if ( bounds(i) .le. tol*temp ) then
            nconv = nconv + 1
         end if
c
   10 continue
c 
      call arscnd (t1)
      tsconv = tsconv + (t1 - t0)
c 
      return
c
c     %---------------%
c     | End of ssconv |
c     %---------------%
c
      end
