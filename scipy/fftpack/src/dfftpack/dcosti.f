      SUBROUTINE DCOSTI (N,WSAVE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       WSAVE(1)
      DATA PI /3.14159265358979323846D0/
      IF (N .LE. 3) RETURN
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      DT = PI/FLOAT(NM1)
      FK = 0.0D0
      DO 101 K=2,NS2
         KC = NP1-K
         FK = FK+1.0D0
         WSAVE(K) = 2.0D0*SIN(FK*DT)
         WSAVE(KC) = 2.0D0*COS(FK*DT)
  101 CONTINUE
      CALL DFFTI (NM1,WSAVE(N+1))
      RETURN
      END
