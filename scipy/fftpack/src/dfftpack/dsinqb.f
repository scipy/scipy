      SUBROUTINE DSINQB (N,X,WSAVE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       X(1)       ,WSAVE(1)
      IF (N .GT. 1) GO TO 101
      X(1) = 4.0D0*X(1)
      RETURN
  101 NS2 = N/2
      DO 102 K=2,N,2
         X(K) = -X(K)
  102 CONTINUE
      CALL DCOSQB (N,X,WSAVE)
      DO 103 K=1,NS2
         KC = N-K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
  103 CONTINUE
      RETURN
      END
