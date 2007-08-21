      SUBROUTINE cumgam(x,a,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE CUMGAM(X,A,CUM,CCUM)
C           Double precision cUMulative incomplete GAMma distribution
C
C
C                              Function
C
C
C     Computes   the  cumulative        of    the     incomplete   gamma
C     distribution, i.e., the integral from 0 to X of
C          (1/GAM(A))*EXP(-T)*T**(A-1) DT
C     where GAM(A) is the complete gamma function of A, i.e.,
C          GAM(A) = integral from 0 to infinity of
C                    EXP(-T)*T**(A-1) DT
C
C
C                              Arguments
C
C
C     X --> The upper limit of integration of the incomplete gamma.
C                                                X is DOUBLE PRECISION
C
C     A --> The shape parameter of the incomplete gamma.
C                                                A is DOUBLE PRECISION
C
C     CUM <-- Cumulative incomplete gamma distribution.
C                                        CUM is DOUBLE PRECISION
C
C     CCUM <-- Compliment of Cumulative incomplete gamma distribution.
C                                                CCUM is DOUBLE PRECISIO
C
C
C                              Method
C
C
C     Calls the routine GRATIO.
C
C**********************************************************************
C
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,x,cum,ccum
C     ..
C     .. External Routines ..
      EXTERNAL gratio
C     ..
C     .. Executable Statements ..
      IF (.NOT. (x.LE.0.0D0)) GO TO 10
      cum = 0.0D0
      ccum = 1.0D0
      RETURN

   10 CALL gratio(a,x,cum,ccum,0)

C     Call gratio routine

      RETURN

      END
