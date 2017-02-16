      SUBROUTINE cumchi(x,df,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE FUNCTION CUMCHI(X,DF,CUM,CCUM)
C             CUMulative of the CHi-square distribution
C
C
C                              Function
C
C
C     Calculates the cumulative chi-square distribution.
C
C
C                              Arguments
C
C
C     X       --> Upper limit of integration of the
C                 chi-square distribution.
C                                                 X is DOUBLE PRECISION
C
C     DF      --> Degrees of freedom of the
C                 chi-square distribution.
C                                                 DF is DOUBLE PRECISION
C
C     CUM <-- Cumulative chi-square distribution.
C                                                 CUM is DOUBLE PRECISIO
C
C     CCUM <-- Compliment of Cumulative chi-square distribution.
C                                                 CCUM is DOUBLE PRECISI
C
C
C                              Method
C
C
C     Calls incomplete gamma function (CUMGAM)
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION df,x,cum,ccum
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,xx
C     ..
C     .. External Subroutines ..
      EXTERNAL cumgam
C     ..
C     .. Executable Statements ..
      a = df*0.5D0
      xx = x*0.5D0
      CALL cumgam(xx,a,cum,ccum)
      RETURN

      END
