      SUBROUTINE RFFTB (N,R,WSAVE)
      DIMENSION       R(*)       ,WSAVE(*)
      IF (N .EQ. 1) RETURN
      CALL RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
