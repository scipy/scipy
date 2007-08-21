      DOUBLE PRECISION FUNCTION rlog(x)
C     -------------------
C     COMPUTATION OF  X - 1 - LN(X)
C     -------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,b,p0,p1,p2,q1,q2,r,t,u,w,w1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,dlog
C     ..
C     .. Data statements ..
C     -------------------
      DATA a/.566749439387324D-01/
      DATA b/.456512608815524D-01/
      DATA p0/.333333333333333D+00/,p1/-.224696413112536D+00/,
     +     p2/.620886815375787D-02/
      DATA q1/-.127408923933623D+01/,q2/.354508718369557D+00/
C     ..
C     .. Executable Statements ..
C     -------------------
      IF (x.LT.0.61D0 .OR. x.GT.1.57D0) GO TO 40
      IF (x.LT.0.82D0) GO TO 10
      IF (x.GT.1.18D0) GO TO 20
C
C              ARGUMENT REDUCTION
C
      u = (x-0.5D0) - 0.5D0
      w1 = 0.0D0
      GO TO 30
C
   10 u = dble(x) - 0.7D0
      u = u/0.7D0
      w1 = a - u*0.3D0
      GO TO 30
C
   20 u = 0.75D0*dble(x) - 1.D0
      w1 = b + u/3.0D0
C
C               SERIES EXPANSION
C
   30 r = u/ (u+2.0D0)
      t = r*r
      w = ((p2*t+p1)*t+p0)/ ((q2*t+q1)*t+1.0D0)
      rlog = 2.0D0*t* (1.0D0/ (1.0D0-r)-r*w) + w1
      RETURN
C
C
   40 r = (x-0.5D0) - 0.5D0
      rlog = r - dlog(x)
      RETURN

      END
