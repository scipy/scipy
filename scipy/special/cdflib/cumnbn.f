      SUBROUTINE cumnbn(s,xn,pr,ompr,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE CUMNNBN(S,XN,PR,OMPR,CUM,CCUM)
C                    CUmulative Negative BINomial distribution
C
C
C                              Function
C
C
C     Returns the probability that it there will be S or fewer failures
C     before there are XN successes, with each binomial trial having
C     a probability of success PR.
C
C     Prob(# failures = S | XN successes, PR)  =
C                        ( XN + S - 1 )
C                        (            ) * PR^XN * (1-PR)^S
C                        (      S     )
C
C
C                              Arguments
C
C
C     S --> The number of failures
C                                                  S is DOUBLE PRECISION
C
C     XN --> The number of successes
C                                                  XN is DOUBLE PRECISIO
C
C     PR --> The probability of success in each binomial trial.
C                                                  PR is DOUBLE PRECISIO
C
C     OMPR --> 1 - PR
C                                                  OMPR is DOUBLE PRECIS
C
C     CUM <-- Cumulative negative binomial distribution.
C                                                  CUM is DOUBLE PRECISI
C
C     CCUM <-- Compliment of Cumulative negative binomial distribution.
C                                                  CCUM is DOUBLE PRECIS
C
C
C                              Method
C
C
C     Formula  26.5.26    of   Abramowitz  and    Stegun,  Handbook   of
C     Mathematical   Functions (1966) is   used  to reduce the  negative
C     binomial distribution to the cumulative beta distribution.
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION pr,ompr,s,xn,cum,ccum
C     ..
C     .. External Subroutines ..
      EXTERNAL cumbet
C     ..
C     .. Executable Statements ..
      CALL cumbet(pr,ompr,xn,s+1.D0,cum,ccum)
      RETURN

      END
