      DOUBLE PRECISION FUNCTION betaln(a0,b0)
C-----------------------------------------------------------------------
C     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
C-----------------------------------------------------------------------
C     E = 0.5*LN(2*PI)
C--------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a0,b0
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,b,c,e,h,u,v,w,z
      INTEGER i,n
C     ..
C     .. External Functions ..
      DOUBLE PRECISION algdiv,alnrel,bcorr,gamln,gsumln
      EXTERNAL algdiv,alnrel,bcorr,gamln,gsumln
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dlog,dmax1,dmin1
C     ..
C     .. Data statements ..
      DATA e/.918938533204673D0/
C     ..
C     .. Executable Statements ..
C--------------------------
      a = dmin1(a0,b0)
      b = dmax1(a0,b0)
      IF (a.GE.8.0D0) GO TO 100
      IF (a.GE.1.0D0) GO TO 20
C-----------------------------------------------------------------------
C                   PROCEDURE WHEN A .LT. 1
C-----------------------------------------------------------------------
      IF (b.GE.8.0D0) GO TO 10
      betaln = gamln(a) + (gamln(b)-gamln(a+b))
      RETURN

   10 betaln = gamln(a) + algdiv(a,b)
      RETURN
C-----------------------------------------------------------------------
C                PROCEDURE WHEN 1 .LE. A .LT. 8
C-----------------------------------------------------------------------
   20 IF (a.GT.2.0D0) GO TO 40
      IF (b.GT.2.0D0) GO TO 30
      betaln = gamln(a) + gamln(b) - gsumln(a,b)
      RETURN

   30 w = 0.0D0
      IF (b.LT.8.0D0) GO TO 60
      betaln = gamln(a) + algdiv(a,b)
      RETURN
C
C                REDUCTION OF A WHEN B .LE. 1000
C
   40 IF (b.GT.1000.0D0) GO TO 80
      n = a - 1.0D0
      w = 1.0D0
      DO 50 i = 1,n
          a = a - 1.0D0
          h = a/b
          w = w* (h/ (1.0D0+h))
   50 CONTINUE
      w = dlog(w)
      IF (b.LT.8.0D0) GO TO 60
      betaln = w + gamln(a) + algdiv(a,b)
      RETURN
C
C                 REDUCTION OF B WHEN B .LT. 8
C
   60 n = b - 1.0D0
      z = 1.0D0
      DO 70 i = 1,n
          b = b - 1.0D0
          z = z* (b/ (a+b))
   70 CONTINUE
      betaln = w + dlog(z) + (gamln(a)+ (gamln(b)-gsumln(a,b)))
      RETURN
C
C                REDUCTION OF A WHEN B .GT. 1000
C
   80 n = a - 1.0D0
      w = 1.0D0
      DO 90 i = 1,n
          a = a - 1.0D0
          w = w* (a/ (1.0D0+a/b))
   90 CONTINUE
      betaln = (dlog(w)-n*dlog(b)) + (gamln(a)+algdiv(a,b))
      RETURN
C-----------------------------------------------------------------------
C                   PROCEDURE WHEN A .GE. 8
C-----------------------------------------------------------------------
  100 w = bcorr(a,b)
      h = a/b
      c = h/ (1.0D0+h)
      u = - (a-0.5D0)*dlog(c)
      v = b*alnrel(h)
      IF (u.LE.v) GO TO 110
      betaln = (((-0.5D0*dlog(b)+e)+w)-v) - u
      RETURN

  110 betaln = (((-0.5D0*dlog(b)+e)+w)-u) - v
      RETURN

      END
