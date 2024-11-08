*
*  Routine:    CMOUT
*
*  Purpose:    Complex matrix output routine.
*
*  Usage:      CALL CMOUT (LOUT, M, N, A, LDA, IDIGIT, IFMT)
*
*  Arguments
*     M      - Number of rows of A.  (Input)
*     N      - Number of columns of A.  (Input)
*     A      - Complex M by N matrix to be printed.  (Input)
*     LDA    - Leading dimension of A exactly as specified in the
*              dimension statement of the calling program.  (Input)
*     IFMT   - Format to be used in printing matrix A.  (Input)
*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
*              If IDIGIT .LT. 0, printing is done with 72 columns.
*              If IDIGIT .GT. 0, printing is done with 132 columns.
*
*\SCCS Information: @(#)
* FILE: cmout.f   SID: 2.1   DATE OF SID: 11/16/95   RELEASE: 2
*
*-----------------------------------------------------------------------
*
      SUBROUTINE CMOUT( LOUT, M, N, A, LDA, IDIGIT, IFMT )
*     ...
*     ... SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M, N, IDIGIT, LDA, LOUT
      Complex
     &                   A( LDA, * )
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
            DO 40 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               WRITE( LOUT, 9998 )( ICOL, I, I = K1, K2 )
               DO 30 I = 1, M
                  IF (K1.NE.N) THEN
                     WRITE( LOUT, 9994 )I, ( A( I, J ), J = K1, K2 )
                  ELSE
                     WRITE( LOUT, 9984 )I, ( A( I, J ), J = K1, K2 )
                  END IF
   30          CONTINUE
   40       CONTINUE
*
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 60 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               WRITE( LOUT, 9997 )( ICOL, I, I = K1, K2 )
               DO 50 I = 1, M
                  IF (K1.NE.N) THEN
                     WRITE( LOUT, 9993 )I, ( A( I, J ), J = K1, K2 )
                  ELSE
                     WRITE( LOUT, 9983 )I, ( A( I, J ), J = K1, K2 )
                  END IF
   50          CONTINUE
   60       CONTINUE
*
         ELSE IF( NDIGIT.LE.8 ) THEN
            DO 80 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               WRITE( LOUT, 9996 )( ICOL, I, I = K1, K2 )
               DO 70 I = 1, M
                  IF (K1.NE.N) THEN
                     WRITE( LOUT, 9992 )I, ( A( I, J ), J = K1, K2 )
                  ELSE
                     WRITE( LOUT, 9982 )I, ( A( I, J ), J = K1, K2 )
                  END IF
   70          CONTINUE
   80       CONTINUE
*
         ELSE
            DO 100 K1 = 1, N
               WRITE( LOUT, 9995 ) ICOL, K1
               DO 90 I = 1, M
                  WRITE( LOUT, 9991 )I, A( I, K1 )
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
            DO 120 K1 = 1, N, 4
               K2 = MIN0( N, K1+3 )
               WRITE( LOUT, 9998 )( ICOL, I, I = K1, K2 )
               DO 110 I = 1, M
                  IF ((K1+3).LE.N) THEN
                     WRITE( LOUT, 9974 )I, ( A( I, J ), J = K1, K2 )
                  ELSE IF ((K1+3-N).EQ.1) THEN
                     WRITE( LOUT, 9964 )I, ( A( I, J ), J = k1, K2 )
                  ELSE IF ((K1+3-N).EQ.2) THEN
                     WRITE( LOUT, 9954 )I, ( A( I, J ), J = K1, K2 )
                  ELSE IF ((K1+3-N).EQ.3) THEN
                     WRITE( LOUT, 9944 )I, ( A( I, J ), J = K1, K2 )
                  END IF
  110          CONTINUE
  120       CONTINUE
*
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 140 K1 = 1, N, 3
               K2 = MIN0( N, K1+ 2)
               WRITE( LOUT, 9997 )( ICOL, I, I = K1, K2 )
               DO 130 I = 1, M
                  IF ((K1+2).LE.N) THEN
                     WRITE( LOUT, 9973 )I, ( A( I, J ), J = K1, K2 )
                  ELSE IF ((K1+2-N).EQ.1) THEN
                     WRITE( LOUT, 9963 )I, ( A( I, J ), J = K1, K2 )
                  ELSE IF ((K1+2-N).EQ.2) THEN
                     WRITE( LOUT, 9953 )I, ( A( I, J ), J = K1, K2 )
                  END IF
  130          CONTINUE
  140       CONTINUE
*
         ELSE IF( NDIGIT.LE.8 ) THEN
            DO 160 K1 = 1, N, 3
               K2 = MIN0( N, K1+2 )
                  WRITE( LOUT, 9996 )( ICOL, I, I = K1, K2 )
               DO 150 I = 1, M
                  IF ((K1+2).LE.N) THEN
                     WRITE( LOUT, 9972 )I, ( A( I, J ), J = K1, K2 )
                  ELSE IF ((K1+2-N).EQ.1) THEN
                     WRITE( LOUT, 9962 )I, ( A( I, J ), J = K1, K2 )
                  ELSE IF ((K1+2-N).EQ.2) THEN
                     WRITE( LOUT, 9952 )I, ( A( I, J ), J = K1, K2 )
                  END IF
  150          CONTINUE
  160       CONTINUE
*
         ELSE
            DO 180 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               WRITE( LOUT, 9995 )( ICOL, I, I = K1, K2 )
               DO 170 I = 1, M
                  IF ((K1+1).LE.N) THEN
                     WRITE( LOUT, 9971 )I, ( A( I, J ), J = K1, K2 )
                  ELSE
                     WRITE( LOUT, 9961 )I, ( A( I, J ), J = K1, K2 )
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      END IF
      WRITE( LOUT, 9990 )
*
 9998 FORMAT( 11X, 4( 9X, 3A1, I4, 9X ) )
 9997 FORMAT( 10X, 4( 11X, 3A1, I4, 11X ) )
 9996 FORMAT( 10X, 3( 13X, 3A1, I4, 13X ) )
 9995 FORMAT( 12X, 2( 18x, 3A1, I4, 18X ) )
*
*========================================================
*              FORMAT FOR 72 COLUMN
*========================================================
*
*            DISPLAY 4 SIGNIFICANT DIGITS
*
 9994 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,2('(',E10.3,',',E10.3,')  ') )
 9984 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,1('(',E10.3,',',E10.3,')  ') )
*
*            DISPLAY 6 SIGNIFICANT DIGITS
*
 9993 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,2('(',E12.5,',',E12.5,')  ') )
 9983 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,1('(',E12.5,',',E12.5,')  ') )
*
*            DISPLAY 8 SIGNIFICANT DIGITS
*
 9992 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,2('(',E14.7,',',E14.7,')  ') )
 9982 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,1('(',E14.7,',',E14.7,')  ') )
*
*            DISPLAY 13 SIGNIFICANT DIGITS
*
 9991 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,1('(',E20.13,',',E20.13,')') )
 9990 FORMAT( 1X, ' ' )
*
*
*========================================================
*              FORMAT FOR 132 COLUMN
*========================================================
*
*            DISPLAY 4 SIGNIFICANT DIGIT
*
 9974 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,4('(',E10.3,',',E10.3,')  ') )
 9964 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,3('(',E10.3,',',E10.3,')  ') )
 9954 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,2('(',E10.3,',',E10.3,')  ') )
 9944 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,1('(',E10.3,',',E10.3,')  ') )
*
*            DISPLAY 6 SIGNIFICANT DIGIT
*
 9973 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,3('(',E12.5,',',E12.5,')  ') )
 9963 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,2('(',E12.5,',',E12.5,')  ') )
 9953 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,1('(',E12.5,',',E12.5,')  ') )
*
*            DISPLAY 8 SIGNIFICANT DIGIT
*
 9972 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,3('(',E14.7,',',E14.7,')  ') )
 9962 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,2('(',E14.7,',',E14.7,')  ') )
 9952 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,1('(',E14.7,',',E14.7,')  ') )
*
*            DISPLAY 13 SIGNIFICANT DIGIT
*
 9971 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,2('(',E20.13,',',E20.13,
     &        ')  '))
 9961 FORMAT( 1X, ' Row', I4, ':', 1X, 1P,1('(',E20.13,',',E20.13,
     &        ')  '))

*
*
*
*
      RETURN
      END
