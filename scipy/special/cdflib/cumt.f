      SUBROUTINE cumt(t,df,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE CUMT(T,DF,CUM,CCUM)
C                    CUMulative T-distribution
C
C
C                              Function
C
C
C     Computes the integral from -infinity to T of the t-density.
C
C
C                              Arguments
C
C
C     T --> Upper limit of integration of the t-density.
C                                                  T is DOUBLE PRECISION
C
C     DF --> Degrees of freedom of the t-distribution.
C                                                  DF is DOUBLE PRECISIO
C
C     CUM <-- Cumulative t-distribution.
C                                                  CCUM is DOUBLE PRECIS
C
C     CCUM <-- Compliment of Cumulative t-distribution.
C                                                  CCUM is DOUBLE PRECIS
C
C
C                              Method
C
C
C     Formula 26.5.27   of     Abramowitz  and   Stegun,    Handbook  of
C     Mathematical Functions  is   used   to  reduce the  t-distribution
C     to an incomplete beta.
C
C**********************************************************************

C     .. Scalar Arguments ..
      DOUBLE PRECISION df,t,cum,ccum
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION xx,a,oma,tt,yy,dfptt
C     ..
C     .. External Subroutines ..
      EXTERNAL cumbet
C     ..
C     .. Executable Statements ..
      tt = t*t
      dfptt = df + tt
      xx = df/dfptt
      yy = tt/dfptt
      CALL cumbet(xx,yy,0.5D0*df,0.5D0,a,oma)
      IF (.NOT. (t.LE.0.0D0)) GO TO 10
      cum = 0.5D0*a
      ccum = oma + cum
      GO TO 20

   10 ccum = 0.5D0*a
      cum = oma + ccum
   20 RETURN

      END
