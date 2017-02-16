      SUBROUTINE SINTI (N,WSAVE)
      DIMENSION       WSAVE(*)
      DATA PI /3.14159265358979/
      IF (N .LE. 1) RETURN
      NS2 = N/2
      NP1 = N+1
      DT = PI/FLOAT(NP1)
      DO 101 K=1,NS2
         WSAVE(K) = 2.*SIN(K*DT)
  101 CONTINUE
      CALL RFFTI (NP1,WSAVE(NS2+1))
      RETURN
      END
