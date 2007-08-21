      SUBROUTINE cumbin(s,xn,pr,ompr,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE CUMBIN(S,XN,PBIN,OMPR,CUM,CCUM)
C                    CUmulative BINomial distribution
C
C
C                              Function
C
C
C     Returns the probability   of 0  to  S  successes in  XN   binomial
C     trials, each of which has a probability of success, PBIN.
C
C
C                              Arguments
C
C
C     S --> The upper limit of cumulation of the binomial distribution.
C                                                  S is DOUBLE PRECISION
C
C     XN --> The number of binomial trials.
C                                                  XN is DOUBLE PRECISIO
C
C     PBIN --> The probability of success in each binomial trial.
C                                                  PBIN is DOUBLE PRECIS
C
C     OMPR --> 1 - PBIN
C                                                  OMPR is DOUBLE PRECIS
C
C     CUM <-- Cumulative binomial distribution.
C                                                  CUM is DOUBLE PRECISI
C
C     CCUM <-- Compliment of Cumulative binomial distribution.
C                                                  CCUM is DOUBLE PRECIS

C
C
C                              Method
C
C
C     Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of
C     Mathematical   Functions (1966) is   used  to reduce the  binomial
C     distribution  to  the  cumulative    beta distribution.
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION pr,ompr,s,xn,cum,ccum
C     ..
C     .. External Subroutines ..
      EXTERNAL cumbet
C     ..
C     .. Executable Statements ..
      IF (.NOT. (s.LT.xn)) GO TO 10
      CALL cumbet(pr,ompr,s+1.0D0,xn-s,ccum,cum)
      GO TO 20

   10 cum = 1.0D0
      ccum = 0.0D0
   20 RETURN

      END
