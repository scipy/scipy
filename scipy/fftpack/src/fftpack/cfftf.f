      SUBROUTINE CFFTF (N,C,WSAVE)
      IMPLICIT NONE
      INTEGER :: N, IW1, IW2
      REAL :: C(*), WSAVE(*)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
