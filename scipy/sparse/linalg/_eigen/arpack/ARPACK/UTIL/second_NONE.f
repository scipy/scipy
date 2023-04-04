      SUBROUTINE ARSCND( T )
*
      REAL       T
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     July 26, 1991
*
*  Purpose
*  =======
*
*  SECOND returns the user time for a process in arscnds.
*  This version gets the time from the system function ETIME.
*
*     .. Local Scalars ..
      REAL               T1
*     ..
*     .. Local Arrays ..
      REAL               TARRAY( 2 )
*     ..
*     .. External Functions ..
*     ..
*     .. Executable Statements ..
*

c      T1 = ETIME( TARRAY )
c      T  = TARRAY( 1 )
      T = 0
      RETURN
*
*     End of ARSCND
*
      END
