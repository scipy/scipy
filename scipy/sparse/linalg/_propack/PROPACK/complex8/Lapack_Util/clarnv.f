      SUBROUTINE CLARNV( IDIST, ISEED, N, X )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            IDIST, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      COMPLEX            X( * )
*     ..
*
*  Purpose
*  =======
*
*  CLARNV returns a vector of n random complex numbers from a uniform or
*  normal distribution.
*
*  Arguments
*  =========
*
*  IDIST   (input) INTEGER
*          Specifies the distribution of the random numbers:
*          = 1:  real and imaginary parts each uniform (0,1)
*          = 2:  real and imaginary parts each uniform (-1,1)
*          = 3:  real and imaginary parts each normal (0,1)
*          = 4:  uniformly distributed on the disc abs(z) < 1
*          = 5:  uniformly distributed on the circle abs(z) = 1
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  N       (input) INTEGER
*          The number of random numbers to be generated.
*
*  X       (output) COMPLEX array, dimension (N)
*          The generated random numbers.
*
*  Further Details
*  ===============
*
*  This routine calls the auxiliary routine SLARUV to generate random
*  real numbers from a uniform (0,1) distribution, in batches of up to
*  128 using vectorisable code. The Box-Muller method is used to
*  transform numbers from a uniform to a normal distribution.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
      INTEGER            LV
      PARAMETER          ( LV = 128 )
      REAL               TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IL, IV
*     ..
*     .. Local Arrays ..
      REAL               U( LV )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, EXP, LOG, MIN, SQRT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARUV
*     ..
*     .. Executable Statements ..
*
      DO 60 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
*
*        Call SLARUV to generate 2*IL real numbers from a uniform (0,1)
*        distribution (2*IL <= LV)
*
         CALL SLARUV( ISEED, 2*IL, U )
*
         IF( IDIST.EQ.1 ) THEN
*
*           Copy generated numbers
*
            DO 10 I = 1, IL
               X( IV+I-1 ) = CMPLX( U( 2*I-1 ), U( 2*I ) )
   10       CONTINUE
         ELSE IF( IDIST.EQ.2 ) THEN
*
*           Convert generated numbers to uniform (-1,1) distribution
*
            DO 20 I = 1, IL
               X( IV+I-1 ) = CMPLX( TWO*U( 2*I-1 )-ONE,
     $                       TWO*U( 2*I )-ONE )
   20       CONTINUE
         ELSE IF( IDIST.EQ.3 ) THEN
*
*           Convert generated numbers to normal (0,1) distribution
*
            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )*
     $                       EXP( CMPLX( ZERO, TWOPI*U( 2*I ) ) )
   30       CONTINUE
         ELSE IF( IDIST.EQ.4 ) THEN
*
*           Convert generated numbers to complex numbers uniformly
*           distributed on the unit disk
*
            DO 40 I = 1, IL
               X( IV+I-1 ) = SQRT( U( 2*I-1 ) )*
     $                       EXP( CMPLX( ZERO, TWOPI*U( 2*I ) ) )
   40       CONTINUE
         ELSE IF( IDIST.EQ.5 ) THEN
*
*           Convert generated numbers to complex numbers uniformly
*           distributed on the unit circle
*
            DO 50 I = 1, IL
               X( IV+I-1 ) = EXP( CMPLX( ZERO, TWOPI*U( 2*I ) ) )
   50       CONTINUE
         END IF
   60 CONTINUE
      RETURN
*
*     End of CLARNV
*
      END
