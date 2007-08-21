      DOUBLE PRECISION FUNCTION psi(xx)
C---------------------------------------------------------------------
C
C                 EVALUATION OF THE DIGAMMA FUNCTION
C
C                           -----------
C
C     PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT
C     BE COMPUTED.
C
C     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV
C     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY
C     CODY, STRECOK AND THACHER.
C
C---------------------------------------------------------------------
C     PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK
C     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY
C     A.H. MORRIS (NSWC).
C---------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION xx
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION aug,den,dx0,piov4,sgn,upper,w,x,xmax1,xmx0,
     +                 xsmall,z
      INTEGER i,m,n,nq
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION p1(7),p2(4),q1(6),q2(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION spmpar
      INTEGER ipmpar
      EXTERNAL spmpar,ipmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,cos,dble,dlog,dmin1,int,sin
C     ..
C     .. Data statements ..
C---------------------------------------------------------------------
C
C     PIOV4 = PI/4
C     DX0 = ZERO OF PSI TO EXTENDED PRECISION
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
C     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
C     PSI(X) / (X - X0),  0.5 .LE. X .LE. 3.0
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
C     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
C     PSI(X) - LN(X) + 1 / (2*X),  X .GT. 3.0
C
C---------------------------------------------------------------------
      DATA piov4/.785398163397448D0/
      DATA dx0/1.461632144968362341262659542325721325D0/
      DATA p1(1)/.895385022981970D-02/,p1(2)/.477762828042627D+01/,
     +     p1(3)/.142441585084029D+03/,p1(4)/.118645200713425D+04/,
     +     p1(5)/.363351846806499D+04/,p1(6)/.413810161269013D+04/,
     +     p1(7)/.130560269827897D+04/
      DATA q1(1)/.448452573429826D+02/,q1(2)/.520752771467162D+03/,
     +     q1(3)/.221000799247830D+04/,q1(4)/.364127349079381D+04/,
     +     q1(5)/.190831076596300D+04/,q1(6)/.691091682714533D-05/
      DATA p2(1)/-.212940445131011D+01/,p2(2)/-.701677227766759D+01/,
     +     p2(3)/-.448616543918019D+01/,p2(4)/-.648157123766197D+00/
      DATA q2(1)/.322703493791143D+02/,q2(2)/.892920700481861D+02/,
     +     q2(3)/.546117738103215D+02/,q2(4)/.777788548522962D+01/
C     ..
C     .. Executable Statements ..
C---------------------------------------------------------------------
C
C     MACHINE DEPENDENT CONSTANTS ...
C
C        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
C                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
C                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
C                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
C                 PSI MAY BE REPRESENTED AS ALOG(X).
C
C        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
C                 MAY BE REPRESENTED BY 1/X.
C
C---------------------------------------------------------------------
      xmax1 = ipmpar(3)
      xmax1 = dmin1(xmax1,1.0D0/spmpar(1))
      xsmall = 1.D-9
C---------------------------------------------------------------------
      x = xx
      aug = 0.0D0
      IF (x.GE.0.5D0) GO TO 50
C---------------------------------------------------------------------
C     X .LT. 0.5,  USE REFLECTION FORMULA
C     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
C---------------------------------------------------------------------
      IF (abs(x).GT.xsmall) GO TO 10
      IF (x.EQ.0.0D0) GO TO 100
C---------------------------------------------------------------------
C     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
C     FOR  PI*COTAN(PI*X)
C---------------------------------------------------------------------
      aug = -1.0D0/x
      GO TO 40
C---------------------------------------------------------------------
C     REDUCTION OF ARGUMENT FOR COTAN
C---------------------------------------------------------------------
   10 w = -x
      sgn = piov4
      IF (w.GT.0.0D0) GO TO 20
      w = -w
      sgn = -sgn
C---------------------------------------------------------------------
C     MAKE AN ERROR EXIT IF X .LE. -XMAX1
C---------------------------------------------------------------------
   20 IF (w.GE.xmax1) GO TO 100
      nq = int(w)
      w = w - dble(nq)
      nq = int(w*4.0D0)
      w = 4.0D0* (w-dble(nq)*.25D0)
C---------------------------------------------------------------------
C     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
C     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
C     QUADRANT AND DETERMINE SIGN
C---------------------------------------------------------------------
      n = nq/2
      IF ((n+n).NE.nq) w = 1.0D0 - w
      z = piov4*w
      m = n/2
      IF ((m+m).NE.n) sgn = -sgn
C---------------------------------------------------------------------
C     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
C---------------------------------------------------------------------
      n = (nq+1)/2
      m = n/2
      m = m + m
      IF (m.NE.n) GO TO 30
C---------------------------------------------------------------------
C     CHECK FOR SINGULARITY
C---------------------------------------------------------------------
      IF (z.EQ.0.0D0) GO TO 100
C---------------------------------------------------------------------
C     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
C     SIN/COS AS A SUBSTITUTE FOR TAN
C---------------------------------------------------------------------
      aug = sgn* ((cos(z)/sin(z))*4.0D0)
      GO TO 40

   30 aug = sgn* ((sin(z)/cos(z))*4.0D0)
   40 x = 1.0D0 - x
   50 IF (x.GT.3.0D0) GO TO 70
C---------------------------------------------------------------------
C     0.5 .LE. X .LE. 3.0
C---------------------------------------------------------------------
      den = x
      upper = p1(1)*x
C
      DO 60 i = 1,5
          den = (den+q1(i))*x
          upper = (upper+p1(i+1))*x
   60 CONTINUE
C
      den = (upper+p1(7))/ (den+q1(6))
      xmx0 = dble(x) - dx0
      psi = den*xmx0 + aug
      RETURN
C---------------------------------------------------------------------
C     IF X .GE. XMAX1, PSI = LN(X)
C---------------------------------------------------------------------
   70 IF (x.GE.xmax1) GO TO 90
C---------------------------------------------------------------------
C     3.0 .LT. X .LT. XMAX1
C---------------------------------------------------------------------
      w = 1.0D0/ (x*x)
      den = w
      upper = p2(1)*w
C
      DO 80 i = 1,3
          den = (den+q2(i))*w
          upper = (upper+p2(i+1))*w
   80 CONTINUE
C
      aug = upper/ (den+q2(4)) - 0.5D0/x + aug
   90 psi = aug + dlog(x)
      RETURN
C---------------------------------------------------------------------
C     ERROR RETURN
C---------------------------------------------------------------------
  100 psi = 0.0D0
      RETURN

      END
