*-----------------------------------------------------------------------
*  Routine:    SMOUT
*
*  Purpose:    Real matrix output routine.
*
*  Usage:      CALL SMOUT (LOUT, M, N, A, LDA, IDIGIT, IFMT)
*
*  Arguments
*     M      - Number of rows of A.  (Input)
*     N      - Number of columns of A.  (Input)
*     A      - Real M by N matrix to be printed.  (Input)
*     LDA    - Leading dimension of A exactly as specified in the
*              dimension statement of the calling program.  (Input)
*     IFMT   - Format to be used in printing matrix A.  (Input)
*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
*              If IDIGIT .LT. 0, printing is done with 72 columns.
*              If IDIGIT .GT. 0, printing is done with 132 columns.
*
*-----------------------------------------------------------------------
*
      SUBROUTINE SMOUT( LOUT, M, N, A, LDA, IDIGIT, IFMT )
*     ...
*     ... SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M, N, IDIGIT, LDA, LOUT
      REAL               A( LDA, * )
      CHARACTER          IFMT*( * )
*     ...
*     ... SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I, J, NDIGIT, K1, K2, LLL
      CHARACTER*1        ICOL( 3 )
      CHARACTER*80       LINE
*     ...
*     ... SPECIFICATIONS INTRINSICS
      INTRINSIC          MIN
*
      DATA               ICOL( 1 ), ICOL( 2 ), ICOL( 3 ) / 'C', 'o',
     $                   'l' /
*     ...
*     ... FIRST EXECUTABLE STATEMENT
*
      LLL = MIN( LEN( IFMT ), 80 )
      DO 10 I = 1, LLL
         LINE( I: I ) = '-'
   10 CONTINUE
*
      DO 20 I = LLL + 1, 80
         LINE( I: I ) = ' '
   20 CONTINUE
*
      WRITE( LOUT, 9999 )IFMT, LINE( 1: LLL )
 9999 FORMAT( / 1X, A / 1X, A )
*
      IF( M.LE.0 .OR. N.LE.0 .OR. LDA.LE.0 )
     $   RETURN
      NDIGIT = IDIGIT
      IF( IDIGIT.EQ.0 )
     $   NDIGIT = 4
*
*=======================================================================
*             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
*=======================================================================
*
      IF( IDIGIT.LT.0 ) THEN
         NDIGIT = -IDIGIT
         IF( NDIGIT.LE.4 ) THEN
            DO 40 K1 = 1, N, 5
               K2 = MIN0( N, K1+4 )
               WRITE( LOUT, 9998 )( ICOL, I, I = K1, K2 )
               DO 30 I = 1, M
                  WRITE( LOUT, 9994 )I, ( A( I, J ), J = K1, K2 )
   30          CONTINUE
   40       CONTINUE
*
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 60 K1 = 1, N, 4
               K2 = MIN0( N, K1+3 )
               WRITE( LOUT, 9997 )( ICOL, I, I = K1, K2 )
               DO 50 I = 1, M
                  WRITE( LOUT, 9993 )I, ( A( I, J ), J = K1, K2 )
   50          CONTINUE
   60       CONTINUE
*
         ELSE IF( NDIGIT.LE.10 ) THEN
            DO 80 K1 = 1, N, 3
               K2 = MIN0( N, K1+2 )
               WRITE( LOUT, 9996 )( ICOL, I, I = K1, K2 )
               DO 70 I = 1, M
                  WRITE( LOUT, 9992 )I, ( A( I, J ), J = K1, K2 )
   70          CONTINUE
   80       CONTINUE
*
         ELSE
            DO 100 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               WRITE( LOUT, 9995 )( ICOL, I, I = K1, K2 )
               DO 90 I = 1, M
                  WRITE( LOUT, 9991 )I, ( A( I, J ), J = K1, K2 )
   90          CONTINUE
  100       CONTINUE
         END IF
*
*=======================================================================
*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
*=======================================================================
*
      ELSE
         IF( NDIGIT.LE.4 ) THEN
            DO 120 K1 = 1, N, 10
               K2 = MIN0( N, K1+9 )
               WRITE( LOUT, 9998 )( ICOL, I, I = K1, K2 )
               DO 110 I = 1, M
                  WRITE( LOUT, 9994 )I, ( A( I, J ), J = K1, K2 )
  110          CONTINUE
  120       CONTINUE
*
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 140 K1 = 1, N, 8
               K2 = MIN0( N, K1+7 )
               WRITE( LOUT, 9997 )( ICOL, I, I = K1, K2 )
               DO 130 I = 1, M
                  WRITE( LOUT, 9993 )I, ( A( I, J ), J = K1, K2 )
  130          CONTINUE
  140       CONTINUE
*
         ELSE IF( NDIGIT.LE.10 ) THEN
            DO 160 K1 = 1, N, 6
               K2 = MIN0( N, K1+5 )
               WRITE( LOUT, 9996 )( ICOL, I, I = K1, K2 )
               DO 150 I = 1, M
                  WRITE( LOUT, 9992 )I, ( A( I, J ), J = K1, K2 )
  150          CONTINUE
  160       CONTINUE
*
         ELSE
            DO 180 K1 = 1, N, 5
               K2 = MIN0( N, K1+4 )
               WRITE( LOUT, 9995 )( ICOL, I, I = K1, K2 )
               DO 170 I = 1, M
                  WRITE( LOUT, 9991 )I, ( A( I, J ), J = K1, K2 )
  170          CONTINUE
  180       CONTINUE
         END IF
      END IF
      WRITE( LOUT, 9990 )
*
 9998 FORMAT( 10X, 10( 4X, 3A1, I4, 1X ) )
 9997 FORMAT( 10X, 8( 5X, 3A1, I4, 2X ) )
 9996 FORMAT( 10X, 6( 7X, 3A1, I4, 4X ) )
 9995 FORMAT( 10X, 5( 9X, 3A1, I4, 6X ) )
 9994 FORMAT( 1X, ' Row', I4, ':', 1X, 1P10E12.3 )
 9993 FORMAT( 1X, ' Row', I4, ':', 1X, 1P8E14.5 )
 9992 FORMAT( 1X, ' Row', I4, ':', 1X, 1P6E18.9 )
 9991 FORMAT( 1X, ' Row', I4, ':', 1X, 1P5E22.13 )
 9990 FORMAT( 1X, ' ' )
*
      RETURN
      END
