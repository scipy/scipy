      SUBROUTINE bratio(a,b,x,y,w,w1,ierr)
C-----------------------------------------------------------------------
C
C            EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B)
C
C                     --------------------
C
C     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X .LE. 1
C     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES
C
C                      W  = IX(A,B)
C                      W1 = 1 - IX(A,B)
C
C     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
C     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
C     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED,
C     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO
C     ONE OF THE FOLLOWING VALUES ...
C
C        IERR = 1  IF A OR B IS NEGATIVE
C        IERR = 2  IF A = B = 0
C        IERR = 3  IF X .LT. 0 OR X .GT. 1
C        IERR = 4  IF Y .LT. 0 OR Y .GT. 1
C        IERR = 5  IF X + Y .NE. 1
C        IERR = 6  IF X = A = 0
C        IERR = 7  IF Y = B = 0
C
C--------------------
C     WRITTEN BY ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN, VIRGINIA
C     REVISED ... NOV 1991
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,w,w1,x,y
      INTEGER ierr
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a0,b0,eps,lambda,t,x0,y0,z
      INTEGER ierr1,ind,n
C     ..
C     .. External Functions ..
      DOUBLE PRECISION apser,basym,bfrac,bpser,bup,fpser,spmpar
      EXTERNAL apser,basym,bfrac,bpser,bup,fpser,spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL bgrat
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dmax1,dmin1
C     ..
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
C            FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
C
      eps = spmpar(1)
C
C-----------------------------------------------------------------------
      w = 0.0D0
      w1 = 0.0D0
      IF (a.LT.0.0D0 .OR. b.LT.0.0D0) GO TO 270
      IF (a.EQ.0.0D0 .AND. b.EQ.0.0D0) GO TO 280
      IF (x.LT.0.0D0 .OR. x.GT.1.0D0) GO TO 290
      IF (y.LT.0.0D0 .OR. y.GT.1.0D0) GO TO 300
      z = ((x+y)-0.5D0) - 0.5D0
      IF (abs(z).GT.3.0D0*eps) GO TO 310
C
      ierr = 0
      IF (x.EQ.0.0D0) GO TO 210
      IF (y.EQ.0.0D0) GO TO 230
      IF (a.EQ.0.0D0) GO TO 240
      IF (b.EQ.0.0D0) GO TO 220
C
      eps = dmax1(eps,1.D-15)
      IF (dmax1(a,b).LT.1.D-3*eps) GO TO 260
C
      ind = 0
      a0 = a
      b0 = b
      x0 = x
      y0 = y
      IF (dmin1(a0,b0).GT.1.0D0) GO TO 40
C
C             PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
C
      IF (x.LE.0.5D0) GO TO 10
      ind = 1
      a0 = b
      b0 = a
      x0 = y
      y0 = x
C
   10 IF (b0.LT.dmin1(eps,eps*a0)) GO TO 90
      IF (a0.LT.dmin1(eps,eps*b0) .AND. b0*x0.LE.1.0D0) GO TO 100
      IF (dmax1(a0,b0).GT.1.0D0) GO TO 20
      IF (a0.GE.dmin1(0.2D0,b0)) GO TO 110
      IF (x0**a0.LE.0.9D0) GO TO 110
      IF (x0.GE.0.3D0) GO TO 120
      n = 20
      GO TO 140
C
   20 IF (b0.LE.1.0D0) GO TO 110
      IF (x0.GE.0.3D0) GO TO 120
      IF (x0.GE.0.1D0) GO TO 30
      IF ((x0*b0)**a0.LE.0.7D0) GO TO 110
   30 IF (b0.GT.15.0D0) GO TO 150
      n = 20
      GO TO 140
C
C             PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
C
   40 IF (a.GT.b) GO TO 50
      lambda = a - (a+b)*x
      GO TO 60

   50 lambda = (a+b)*y - b
   60 IF (lambda.GE.0.0D0) GO TO 70
      ind = 1
      a0 = b
      b0 = a
      x0 = y
      y0 = x
      lambda = abs(lambda)
C
   70 IF (b0.LT.40.0D0 .AND. b0*x0.LE.0.7D0) GO TO 110
      IF (b0.LT.40.0D0) GO TO 160
      IF (a0.GT.b0) GO TO 80
      IF (a0.LE.100.0D0) GO TO 130
      IF (lambda.GT.0.03D0*a0) GO TO 130
      GO TO 200

   80 IF (b0.LE.100.0D0) GO TO 130
      IF (lambda.GT.0.03D0*b0) GO TO 130
      GO TO 200
C
C            EVALUATION OF THE APPROPRIATE ALGORITHM
C
   90 w = fpser(a0,b0,x0,eps)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
  100 w1 = apser(a0,b0,x0,eps)
      w = 0.5D0 + (0.5D0-w1)
      GO TO 250
C
  110 w = bpser(a0,b0,x0,eps)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
  120 w1 = bpser(b0,a0,y0,eps)
      w = 0.5D0 + (0.5D0-w1)
      GO TO 250
C
  130 w = bfrac(a0,b0,x0,y0,lambda,15.0D0*eps)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
  140 w1 = bup(b0,a0,y0,x0,n,eps)
      b0 = b0 + n
  150 CALL bgrat(b0,a0,y0,x0,w1,15.0D0*eps,ierr1)
      w = 0.5D0 + (0.5D0-w1)
      GO TO 250
C
  160 n = b0
      b0 = b0 - n
      IF (b0.NE.0.0D0) GO TO 170
      n = n - 1
      b0 = 1.0D0
  170 w = bup(b0,a0,y0,x0,n,eps)
      IF (x0.GT.0.7D0) GO TO 180
      w = w + bpser(a0,b0,x0,eps)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
  180 IF (a0.GT.15.0D0) GO TO 190
      n = 20
      w = w + bup(a0,b0,x0,y0,n,eps)
      a0 = a0 + n
  190 CALL bgrat(a0,b0,x0,y0,w,15.0D0*eps,ierr1)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
  200 w = basym(a0,b0,lambda,100.0D0*eps)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
C               TERMINATION OF THE PROCEDURE
C
  210 IF (a.EQ.0.0D0) GO TO 320
  220 w = 0.0D0
      w1 = 1.0D0
      RETURN
C
  230 IF (b.EQ.0.0D0) GO TO 330
  240 w = 1.0D0
      w1 = 0.0D0
      RETURN
C
  250 IF (ind.EQ.0) RETURN
      t = w
      w = w1
      w1 = t
      RETURN
C
C           PROCEDURE FOR A AND B .LT. 1.E-3*EPS
C
  260 w = b/ (a+b)
      w1 = a/ (a+b)
      RETURN
C
C                       ERROR RETURN
C
  270 ierr = 1
      RETURN

  280 ierr = 2
      RETURN

  290 ierr = 3
      RETURN

  300 ierr = 4
      RETURN

  310 ierr = 5
      RETURN

  320 ierr = 6
      RETURN

  330 ierr = 7
      RETURN

      END
