      DOUBLE PRECISION FUNCTION rexp(x)
C-----------------------------------------------------------------------
C            EVALUATION OF THE FUNCTION EXP(X) - 1
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION p1,p2,q1,q2,q3,q4,w
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp
C     ..
C     .. Data statements ..
      DATA p1/.914041914819518D-09/,p2/.238082361044469D-01/,
     +     q1/-.499999999085958D+00/,q2/.107141568980644D+00/,
     +     q3/-.119041179760821D-01/,q4/.595130811860248D-03/
C     ..
C     .. Executable Statements ..
C-----------------------
      IF (abs(x).GT.0.15D0) GO TO 10
      rexp = x* (((p2*x+p1)*x+1.0D0)/ ((((q4*x+q3)*x+q2)*x+q1)*x+1.0D0))
      RETURN
C
   10 w = exp(x)
      IF (x.GT.0.0D0) GO TO 20
      rexp = (w-0.5D0) - 0.5D0
      RETURN

   20 rexp = w* (0.5D0+ (0.5D0-1.0D0/w))
      RETURN

      END
