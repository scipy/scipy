      SUBROUTINE cumfnc(f,dfn,dfd,pnonc,cum,ccum,status)
C**********************************************************************
C
C               F -NON- -C-ENTRAL F DISTRIBUTION
C
C
C
C                              Function
C
C
C     COMPUTES NONCENTRAL F DISTRIBUTION WITH DFN AND DFD
C     DEGREES OF FREEDOM AND NONCENTRALITY PARAMETER PNONC
C
C
C                              Arguments
C
C
C     X --> UPPER LIMIT OF INTEGRATION OF NONCENTRAL F IN EQUATION
C
C     DFN --> DEGREES OF FREEDOM OF NUMERATOR
C
C     DFD -->  DEGREES OF FREEDOM OF DENOMINATOR
C
C     PNONC --> NONCENTRALITY PARAMETER.
C
C     CUM <-- CUMULATIVE NONCENTRAL F DISTRIBUTION
C
C     CCUM <-- COMPLIMENT OF CUMMULATIVE
C
C
C                              Method
C
C
C     USES FORMULA 26.6.20 OF REFERENCE FOR INFINITE SERIES.
C     SERIES IS CALCULATED BACKWARD AND FORWARD FROM J = LAMBDA/2
C     (THIS IS THE TERM WITH THE LARGEST POISSON WEIGHT) UNTIL
C     THE CONVERGENCE CRITERION IS MET.
C
C     FOR SPEED, THE INCOMPLETE BETA FUNCTIONS ARE EVALUATED
C     BY FORMULA 26.5.16.
C
C
C               REFERENCE
C
C
C     HANDBOOD OF MATHEMATICAL FUNCTIONS
C     EDITED BY MILTON ABRAMOWITZ AND IRENE A. STEGUN
C     NATIONAL BUREAU OF STANDARDS APPLIED MATEMATICS SERIES - 55
C     MARCH 1965
C     P 947, EQUATIONS 26.6.17, 26.6.18
C
C
C                              Note
C
C
C     THE SUM CONTINUES UNTIL A SUCCEEDING TERM IS LESS THAN EPS
C     TIMES THE SUM (OR THE SUM IS LESS THAN 1.0E-20).  EPS IS
C     SET TO 1.0E-4 IN A DATA STATEMENT WHICH CAN BE CHANGED.
C
C**********************************************************************

C     .. Scalar Arguments ..
      DOUBLE PRECISION dfd,dfn,pnonc,f,cum,ccum
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION dsum,dummy,prod,xx,yy
      DOUBLE PRECISION adn,aup,b,betdn,betup,centwt,dnterm,eps,sum,
     +                 upterm,xmult,xnonc,x,abstol
      INTEGER i,icent,ierr,status
C     ..
C     .. External Functions ..
      DOUBLE PRECISION alngam, betaln
      EXTERNAL alngam, betaln
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC log,dble,exp
C     ..
C     .. Statement Functions ..
      LOGICAL qsmall
C     ..
C     .. External Subroutines ..
      EXTERNAL bratio,cumf
C     ..
C     .. Parameters ..
      DOUBLE PRECISION half
      PARAMETER (half=0.5D0)
      DOUBLE PRECISION done
      PARAMETER (done=1.0D0)
C     ..
C     .. Data statements ..
      DATA eps/1.0D-4/
      DATA abstol/1.0D-300/
C     ..
C     .. Statement Function definitions ..
      qsmall(x) = .NOT. (sum .GE. abstol .AND. x .GE. eps*sum)
C     ..
C     .. Executable Statements ..
C
      status = 0
      IF (.NOT. (f.LE.0.0D0)) GO TO 10
      cum = 0.0D0
      ccum = 1.0D0
      RETURN

   10 IF (.NOT. (pnonc.LT.1.0D-10)) GO TO 20
C
C     Handle case in which the non-centrality parameter is
C     (essentially) zero.

      CALL cumf(f,dfn,dfd,cum,ccum)
      RETURN

   20 xnonc = pnonc/2.0D0

C     Calculate the central term of the poisson weighting factor.

      icent = xnonc
      IF (.NOT.(DABS(xnonc-icent).LT.1)) THEN
         status = 1
         RETURN
      ENDIF
      IF (icent.EQ.0) icent = 1

C     Compute central weight term

      centwt = exp(-xnonc+icent*log(xnonc)-alngam(dble(icent+1)))

C     Compute central incomplete beta term
C     Assure that minimum of arg to beta and 1 - arg is computed
C          accurately.

      prod = dfn*f
      dsum = dfd + prod
      yy = dfd/dsum
      IF (yy.GT.half) THEN
          xx = prod/dsum
          yy = done - xx

      ELSE
          xx = done - yy
      END IF

      CALL bratio(dfn*half+dble(icent),dfd*half,xx,yy,betdn,dummy,ierr)
      adn = dfn/2.0D0 + dble(icent)
      aup = adn
      b = dfd/2.0D0
      betup = betdn
      sum = centwt*betdn

C     Now sum terms backward from icent until convergence or all done

      xmult = centwt
      i = icent
      IF (adn.LT.2D0) THEN
          dnterm = exp(alngam(adn+b)-alngam(adn+1.0D0)-alngam(b)+
     +         adn*log(xx)+b*log(yy))
      ELSE
C         Same expression, but avoid problems for large adn
          dnterm = exp(-betaln(adn,b)-log(adn)+
     +         adn*log(xx)+b*log(yy))
      END IF
   30 IF (qsmall(xmult*betdn) .OR. i.LE.0) GO TO 40
      xmult = xmult* (i/xnonc)
      i = i - 1
      adn = adn - 1
      dnterm = (adn+1)/ ((adn+b)*xx)*dnterm
      betdn = betdn + dnterm
      sum = sum + xmult*betdn
      GO TO 30

   40 i = icent + 1

C     Now sum forwards until convergence

      xmult = centwt
      IF ((aup-1+b).EQ.0) THEN
          upterm = exp(-alngam(aup)-alngam(b)+ (aup-1)*log(xx)+
     +             b*log(yy))

      ELSE
          IF (aup.LT.2D0) THEN
              upterm = exp(alngam(aup-1+b)-alngam(aup)-alngam(b)+
     +             (aup-1)*log(xx)+b*log(yy))
          ELSE
C             Same expression, but avoid problems for large aup
              upterm = exp(-betaln(aup-1,b)-log(aup-1)+
     +             (aup-1)*log(xx)+b*log(yy))
          END IF
      END IF

      GO TO 60

   50 IF (qsmall(xmult*betup)) GO TO 70
   60 xmult = xmult* (xnonc/i)
      i = i + 1
      aup = aup + 1
      upterm = (aup+b-2.0D0)*xx/ (aup-1)*upterm
      betup = betup - upterm
      sum = sum + xmult*betup
      GO TO 50

   70 cum = sum

      ccum = 0.5D0 + (0.5D0-cum)
      RETURN

      END
