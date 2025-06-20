c-----------------------------------------------------------------------
c
c\SCCS Information: @(#)
c FILE: cvout.f   SID: 2.1   DATE OF SID: 11/16/95   RELEASE: 2
c
*-----------------------------------------------------------------------
*  Routine:    CVOUT
*
*  Purpose:    Complex vector output routine.
*
*  Usage:      CALL CVOUT (LOUT, N, CX, IDIGIT, IFMT)
*
*  Arguments
*     N      - Length of array CX.  (Input)
*     CX     - Complex array to be printed.  (Input)
*     IFMT   - Format to be used in printing array CX.  (Input)
*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
*              If IDIGIT .LT. 0, printing is done with 72 columns.
*              If IDIGIT .GT. 0, printing is done with 132 columns.
*
*-----------------------------------------------------------------------
*
      SUBROUTINE CVOUT( LOUT, N, CX, IDIGIT, IFMT )
*     ...
*     ... SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N, IDIGIT, LOUT
      Complex
     &                   CX( * )
      CHARACTER          IFMT*( * )
*     ...
*     ... SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I, NDIGIT, K1, K2, LLL
      CHARACTER*80       LINE
*     ...
*     ... FIRST EXECUTABLE STATEMENT
*
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
      IF( N.LE.0 )
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
            DO 30 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               IF (K1.NE.N) THEN
                  WRITE( LOUT, 9998 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE
                  WRITE( LOUT, 9997 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               END IF
   30       CONTINUE
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 40 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               IF (K1.NE.N) THEN
                  WRITE( LOUT, 9988 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE
                  WRITE( LOUT, 9987 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               END IF
   40       CONTINUE
         ELSE IF( NDIGIT.LE.8 ) THEN
            DO 50 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               IF (K1.NE.N) THEN
                  WRITE( LOUT, 9978 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE
                  WRITE( LOUT, 9977 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               END IF
   50       CONTINUE
         ELSE
            DO 60 K1 = 1, N
               WRITE( LOUT, 9968 )K1, K1, CX( I )
   60       CONTINUE
         END IF
*
*=======================================================================
*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
*=======================================================================
*
      ELSE
         IF( NDIGIT.LE.4 ) THEN
            DO 70 K1 = 1, N, 4
               K2 = MIN0( N, K1+3 )
               IF ((K1+3).LE.N) THEN
                  WRITE( LOUT, 9958 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE IF ((K1+3-N) .EQ. 1) THEN
                  WRITE( LOUT, 9957 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE IF ((K1+3-N) .EQ. 2) THEN
                  WRITE( LOUT, 9956 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE IF ((K1+3-N) .EQ. 1) THEN
                  WRITE( LOUT, 9955 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               END IF
   70       CONTINUE
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 80 K1 = 1, N, 3
               K2 = MIN0( N, K1+2 )
               IF ((K1+2).LE.N) THEN
                  WRITE( LOUT, 9948 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE IF ((K1+2-N) .EQ. 1) THEN
                  WRITE( LOUT, 9947 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE IF ((K1+2-N) .EQ. 2) THEN
                  WRITE( LOUT, 9946 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               END IF
   80       CONTINUE
         ELSE IF( NDIGIT.LE.8 ) THEN
            DO 90 K1 = 1, N, 3
               K2 = MIN0( N, K1+2 )
               IF ((K1+2).LE.N) THEN
                  WRITE( LOUT, 9938 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE IF ((K1+2-N) .EQ. 1) THEN
                  WRITE( LOUT, 9937 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE IF ((K1+2-N) .EQ. 2) THEN
                  WRITE( LOUT, 9936 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               END IF
   90       CONTINUE
         ELSE
            DO 100 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               IF ((K1+2).LE.N) THEN
                  WRITE( LOUT, 9928 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               ELSE IF ((K1+2-N) .EQ. 1) THEN
                  WRITE( LOUT, 9927 )K1, K2, ( CX( I ),
     $                   I = K1, K2 )
               END IF
  100       CONTINUE
         END IF
      END IF
      WRITE( LOUT, 9994 )
      RETURN
*
*=======================================================================
*                   FORMAT FOR 72 COLUMNS
*=======================================================================
*
*                 DISPLAY 4 SIGNIFICANT DIGITS
*
 9998 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,2('(',E10.3,',',E10.3,')  ') )
 9997 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,1('(',E10.3,',',E10.3,')  ') )
*
*                 DISPLAY 6 SIGNIFICANT DIGITS
*
 9988 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,2('(',E12.5,',',E12.5,')  ') )
 9987 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,1('(',E12.5,',',E12.5,')  ') )
*
*                 DISPLAY 8 SIGNIFICANT DIGITS
*
 9978 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,2('(',E14.7,',',E14.7,')  ') )
 9977 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,1('(',E14.7,',',E14.7,')  ') )
*
*                 DISPLAY 13 SIGNIFICANT DIGITS
*
 9968 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,1('(',E20.13,',',E20.13,')  ') )
*
*=========================================================================
*                   FORMAT FOR 132 COLUMNS
*=========================================================================
*
*                 DISPLAY 4 SIGNIFICANT DIGITS
*
 9958 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,4('(',E10.3,',',E10.3,')  ') )
 9957 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,3('(',E10.3,',',E10.3,')  ') )
 9956 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,2('(',E10.3,',',E10.3,')  ') )
 9955 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,1('(',E10.3,',',E10.3,')  ') )
*
*                 DISPLAY 6 SIGNIFICANT DIGITS
*
 9948 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,3('(',E12.5,',',E12.5,')  ') )
 9947 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,2('(',E12.5,',',E12.5,')  ') )
 9946 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,1('(',E12.5,',',E12.5,')  ') )
*
*                 DISPLAY 8 SIGNIFICANT DIGITS
*
 9938 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,3('(',E14.7,',',E14.7,')  ') )
 9937 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,2('(',E14.7,',',E14.7,')  ') )
 9936 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,1('(',E14.7,',',E14.7,')  ') )
*
*                 DISPLAY 13 SIGNIFICANT DIGITS
*
 9928 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,2('(',E20.13,',',E20.13,')  ') )
 9927 FORMAT( 1X, I4, ' - ', I4, ':', 1X,
     $        1P,1('(',E20.13,',',E20.13,')  ') )
*
*
*
 9994 FORMAT( 1X, ' ' )
      END
