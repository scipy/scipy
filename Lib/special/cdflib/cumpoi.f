      SUBROUTINE cumpoi(s,xlam,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE CUMPOI(S,XLAM,CUM,CCUM)
C                    CUMulative POIsson distribution
C
C
C                              Function
C
C
C     Returns the  probability  of  S   or  fewer events in  a   Poisson
C     distribution with mean XLAM.
C
C
C                              Arguments
C
C
C     S --> Upper limit of cumulation of the Poisson.
C                                                  S is DOUBLE PRECISION
C
C     XLAM --> Mean of the Poisson distribution.
C                                                  XLAM is DOUBLE PRECIS
C
C     CUM <-- Cumulative poisson distribution.
C                                        CUM is DOUBLE PRECISION
C
C     CCUM <-- Compliment of Cumulative poisson distribution.
C                                                  CCUM is DOUBLE PRECIS
C
C
C                              Method
C
C
C     Uses formula  26.4.21   of   Abramowitz and  Stegun,  Handbook  of
C     Mathematical   Functions  to reduce   the   cumulative Poisson  to
C     the cumulative chi-square distribution.
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION s,xlam,cum,ccum
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION chi,df
C     ..
C     .. External Subroutines ..
      EXTERNAL cumchi
C     ..
C     .. Executable Statements ..
      df = 2.0D0* (s+1.0D0)
      chi = 2.0D0*xlam
      CALL cumchi(chi,df,ccum,cum)
      RETURN

      END
