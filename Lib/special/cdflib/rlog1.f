      DOUBLE PRECISION FUNCTION rlog1(x)
C-----------------------------------------------------------------------
C             EVALUATION OF THE FUNCTION X - LN(1 + X)
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,b,h,p0,p1,p2,q1,q2,r,t,w,w1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,dlog
C     ..
C     .. Data statements ..
C------------------------
      DATA a/.566749439387324D-01/
      DATA b/.456512608815524D-01/
      DATA p0/.333333333333333D+00/,p1/-.224696413112536D+00/,
     +     p2/.620886815375787D-02/
      DATA q1/-.127408923933623D+01/,q2/.354508718369557D+00/
C     ..
C     .. Executable Statements ..
C------------------------
      IF (x.LT.-0.39D0 .OR. x.GT.0.57D0) GO TO 40
      IF (x.LT.-0.18D0) GO TO 10
      IF (x.GT.0.18D0) GO TO 20
C
C              ARGUMENT REDUCTION
C
      h = x
      w1 = 0.0D0
      GO TO 30
C
   10 h = dble(x) + 0.3D0
      h = h/0.7D0
      w1 = a - h*0.3D0
      GO TO 30
C
   20 h = 0.75D0*dble(x) - 0.25D0
      w1 = b + h/3.0D0
C
C               SERIES EXPANSION
C
   30 r = h/ (h+2.0D0)
      t = r*r
      w = ((p2*t+p1)*t+p0)/ ((q2*t+q1)*t+1.0D0)
      rlog1 = 2.0D0*t* (1.0D0/ (1.0D0-r)-r*w) + w1
      RETURN
C
C
   40 w = (x+0.5D0) + 0.5D0
      rlog1 = x - dlog(w)
      RETURN

      END
