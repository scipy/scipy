      DOUBLE PRECISION FUNCTION alngam(x)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION ALNGAM(X)
C                 double precision LN of the GAMma function
C
C
C                              Function
C
C
C     Returns the natural logarithm of GAMMA(X).
C
C
C                              Arguments
C
C
C     X --> value at which scaled log gamma is to be returned
C                    X is DOUBLE PRECISION
C
C
C                              Method
C
C
C     If X .le. 6.0, then use recursion to get X below 3
C     then apply rational approximation number 5236 of
C     Hart et al, Computer Approximations, John Wiley and
C     Sons, NY, 1968.
C
C     If X .gt. 6.0, then use recursion to get X to at least 12 and
C     then use formula 5423 of the same source.
C
C**********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION hln2pi
      PARAMETER (hln2pi=0.91893853320467274178D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION offset,prod,xx
      INTEGER i,n
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION coef(5),scoefd(4),scoefn(9)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION devlpl
      EXTERNAL devlpl
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC log,dble,int
C     ..
C     .. Data statements ..
      DATA scoefn(1)/0.62003838007127258804D2/,
     +     scoefn(2)/0.36036772530024836321D2/,
     +     scoefn(3)/0.20782472531792126786D2/,
     +     scoefn(4)/0.6338067999387272343D1/,
     +     scoefn(5)/0.215994312846059073D1/,
     +     scoefn(6)/0.3980671310203570498D0/,
     +     scoefn(7)/0.1093115956710439502D0/,
     +     scoefn(8)/0.92381945590275995D-2/,
     +     scoefn(9)/0.29737866448101651D-2/
      DATA scoefd(1)/0.62003838007126989331D2/,
     +     scoefd(2)/0.9822521104713994894D1/,
     +     scoefd(3)/-0.8906016659497461257D1/,
     +     scoefd(4)/0.1000000000000000000D1/
      DATA coef(1)/0.83333333333333023564D-1/,
     +     coef(2)/-0.27777777768818808D-2/,
     +     coef(3)/0.79365006754279D-3/,coef(4)/-0.594997310889D-3/,
     +     coef(5)/0.8065880899D-3/
C     ..
C     .. Executable Statements ..
      IF (.NOT. (x.LE.6.0D0)) GO TO 70
      prod = 1.0D0
      xx = x
      IF (.NOT. (x.GT.3.0D0)) GO TO 30
   10 IF (.NOT. (xx.GT.3.0D0)) GO TO 20
      xx = xx - 1.0D0
      prod = prod*xx
      GO TO 10

   20 CONTINUE
   30 IF (.NOT. (x.LT.2.0D0)) GO TO 60
   40 IF (.NOT. (xx.LT.2.0D0)) GO TO 50
      prod = prod/xx
      xx = xx + 1.0D0
      GO TO 40

   50 CONTINUE
   60 alngam = devlpl(scoefn,9,xx-2.0D0)/devlpl(scoefd,4,xx-2.0D0)
C
C
C     COMPUTE RATIONAL APPROXIMATION TO GAMMA(X)
C
C
      alngam = alngam*prod
      alngam = log(alngam)
      GO TO 110

   70 offset = hln2pi
C
C
C     IF NECESSARY MAKE X AT LEAST 12 AND CARRY CORRECTION IN OFFSET
C
C
      n = int(12.0D0-x)
      IF (.NOT. (n.GT.0)) GO TO 90
      prod = 1.0D0
      DO 80,i = 1,n
          prod = prod* (x+dble(i-1))
   80 CONTINUE
      offset = offset - log(prod)
      xx = x + dble(n)
      GO TO 100

   90 xx = x
C
C
C     COMPUTE POWER SERIES
C
C
  100 alngam = devlpl(coef,5,1.0D0/ (xx**2))/xx
      alngam = alngam + offset + (xx-0.5D0)*log(xx) - xx
  110 RETURN

      END
