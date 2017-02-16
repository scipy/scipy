      DOUBLE PRECISION FUNCTION alnrel(a)
C-----------------------------------------------------------------------
C            EVALUATION OF THE FUNCTION LN(1 + A)
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION p1,p2,p3,q1,q2,q3,t,t2,w,x
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog
C     ..
C     .. Data statements ..
      DATA p1/-.129418923021993D+01/,p2/.405303492862024D+00/,
     +     p3/-.178874546012214D-01/
      DATA q1/-.162752256355323D+01/,q2/.747811014037616D+00/,
     +     q3/-.845104217945565D-01/
C     ..
C     .. Executable Statements ..
C--------------------------
      IF (abs(a).GT.0.375D0) GO TO 10
      t = a/ (a+2.0D0)
      t2 = t*t
      w = (((p3*t2+p2)*t2+p1)*t2+1.0D0)/ (((q3*t2+q2)*t2+q1)*t2+1.0D0)
      alnrel = 2.0D0*t*w
      RETURN
C
   10 x = 1.D0 + dble(a)
      alnrel = dlog(x)
      RETURN

      END
