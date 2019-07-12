      SUBROUTINE COSQI (N,WSAVE)
      IMPLICIT NONE
      INTEGER :: N, K
      REAL :: WSAVE(*), DT, FK, PIH
      DATA PIH /1.57079632679491/
      DT = PIH/FLOAT(N)
      FK = 0.
      DO 101 K=1,N
         FK = FK+1.
         WSAVE(K) = COS(FK*DT)
  101 CONTINUE
      CALL RFFTI (N,WSAVE(N+1))
      RETURN
      END
