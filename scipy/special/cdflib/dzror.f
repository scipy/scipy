      SUBROUTINE dzror(status,x,fx,xlo,xhi,qleft,qhi)
C**********************************************************************
C
C     SUBROUTINE DZROR(STATUS, X, FX, XLO, XHI, QLEFT, QHI)
C     Double precision ZeRo of a function -- Reverse Communication
C
C
C                              Function
C
C
C     Performs the zero finding.  STZROR must have been called before
C     this routine in order to set its parameters.
C
C
C                              Arguments
C
C
C     STATUS <--> At the beginning of a zero finding problem, STATUS
C                 should be set to 0 and ZROR invoked.  (The value
C                 of other parameters will be ignored on this call.)
C
C                 When ZROR needs the function evaluated, it will set
C                 STATUS to 1 and return.  The value of the function
C                 should be set in FX and ZROR again called without
C                 changing any of its other parameters.
C
C                 When ZROR has finished without error, it will return
C                 with STATUS 0.  In that case (XLO,XHI) bound the answe
C
C                 If ZROR finds an error (which implies that F(XLO)-Y an
C                 F(XHI)-Y have the same sign, it returns STATUS -1.  In
C                 this case, XLO and XHI are undefined.
C                         INTEGER STATUS
C
C     X <-- The value of X at which F(X) is to be evaluated.
C                         DOUBLE PRECISION X
C
C     FX --> The value of F(X) calculated when ZROR returns with
C            STATUS = 1.
C                         DOUBLE PRECISION FX
C
C     XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
C             inverval in X containing the solution below.
C                         DOUBLE PRECISION XLO
C
C     XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
C             inverval in X containing the solution above.
C                         DOUBLE PRECISION XHI
C
C     QLEFT <-- .TRUE. if the stepping search terminated unsuccessfully
C                at XLO.  If it is .FALSE. the search terminated
C                unsuccessfully at XHI.
C                    QLEFT is LOGICAL
C
C     QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
C              search and .FALSE. if F(X) .LT. Y at the
C              termination of the search.
C                    QHI is LOGICAL
C
C**********************************************************************
      IMPLICIT NONE
C     .. Scalar Arguments ..
      DOUBLE PRECISION fx,x,xhi,xlo,zabstl,zreltl,zxhi,zxlo
      INTEGER status
      LOGICAL qhi,qleft
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,abstol,b,c,d,fa,fb,fc,fd,fda,fdb,m,mb,p,q,
     +                 reltol,tol,w,xxhi,xxlo,zx
      INTEGER ext,i99999
      LOGICAL first,qrzero
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,sign
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION ftol
C     ..
C     .. Statement Function definitions ..
      ftol(zx) = 0.5D0*max(abstol,reltol*abs(zx))
C     ..
C     .. Executable Statements ..

      IF (status.GT.0) GO TO 280
      xlo = xxlo
      xhi = xxhi
      b = xlo
      x = xlo
C     GET-FUNCTION-VALUE
      ASSIGN 10 TO i99999
      GO TO 270

   10 fb = fx
      xlo = xhi
      a = xlo
      x = xlo
C     GET-FUNCTION-VALUE
      ASSIGN 20 TO i99999
      GO TO 270
C
C     Check that F(ZXLO) < 0 < F(ZXHI)  or
C                F(ZXLO) > 0 > F(ZXHI)
C
   20 IF (.NOT. (fb.LT.0.0D0)) GO TO 40
      IF (.NOT. (fx.LT.0.0D0)) GO TO 30
      status = -1
      qleft = fx .LT. fb
      qhi = .FALSE.
      RETURN

   30 CONTINUE
   40 IF (.NOT. (fb.GT.0.0D0)) GO TO 60
      IF (.NOT. (fx.GT.0.0D0)) GO TO 50
      status = -1
      qleft = fx .GT. fb
      qhi = .TRUE.
      RETURN

   50 CONTINUE
   60 fa = fx
C
      first = .TRUE.
   70 c = a
      fc = fa
      ext = 0
   80 IF (.NOT. (abs(fc).LT.abs(fb))) GO TO 100
      IF (.NOT. (c.NE.a)) GO TO 90
      d = a
      fd = fa
   90 a = b
      fa = fb
      xlo = c
      b = xlo
      fb = fc
      c = a
      fc = fa
  100 tol = ftol(xlo)
      m = (c+b)*.5D0
      mb = m - b
      IF (.NOT. (abs(mb).GT.tol)) GO TO 240
      IF (.NOT. (ext.GT.3)) GO TO 110
      w = mb
      GO TO 190

  110 tol = sign(tol,mb)
      p = (b-a)*fb
      IF (.NOT. (first)) GO TO 120
      q = fa - fb
      first = .FALSE.
      GO TO 130

  120 fdb = (fd-fb)/ (d-b)
      fda = (fd-fa)/ (d-a)
      p = fda*p
      q = fdb*fa - fda*fb
  130 IF (.NOT. (p.LT.0.0D0)) GO TO 140
      p = -p
      q = -q
  140 IF (ext.EQ.3) p = p*2.0D0
      IF (.NOT. ((p*1.0D0).EQ.0.0D0.OR.p.LE. (q*tol))) GO TO 150
      w = tol
      GO TO 180

  150 IF (.NOT. (p.LT. (mb*q))) GO TO 160
      w = p/q
      GO TO 170

  160 w = mb
  170 CONTINUE
  180 CONTINUE
  190 d = a
      fd = fa
      a = b
      fa = fb
      b = b + w
      xlo = b
      x = xlo
C     GET-FUNCTION-VALUE
      ASSIGN 200 TO i99999
      GO TO 270

  200 fb = fx
      IF (.NOT. ((fc*fb).GE.0.0D0)) GO TO 210
      GO TO 70

  210 IF (.NOT. (w.EQ.mb)) GO TO 220
      ext = 0
      GO TO 230

  220 ext = ext + 1
  230 GO TO 80

  240 xhi = c
      qrzero = (fc.GE.0.0D0 .AND. fb.LE.0.0D0) .OR.
     +         (fc.LT.0.0D0 .AND. fb.GE.0.0D0)
      IF (.NOT. (qrzero)) GO TO 250
      status = 0
      GO TO 260

  250 status = -1
  260 RETURN

      ENTRY dstzr(zxlo,zxhi,zabstl,zreltl)
C**********************************************************************
C
C     SUBROUTINE DSTZR( XLO, XHI, ABSTOL, RELTOL )
C     Double precision SeT ZeRo finder - Reverse communication version
C
C
C                              Function
C
C
C
C     Sets quantities needed by ZROR.  The function of ZROR
C     and the quantities set is given here.
C
C     Concise Description - Given a function F
C     find XLO such that F(XLO) = 0.
C
C          More Precise Description -
C
C     Input condition. F is a double precision function of a single
C     double precision argument and XLO and XHI are such that
C          F(XLO)*F(XHI)  .LE.  0.0
C
C     If the input condition is met, QRZERO returns .TRUE.
C     and output values of XLO and XHI satisfy the following
C          F(XLO)*F(XHI)  .LE. 0.
C          ABS(F(XLO)  .LE. ABS(F(XHI)
C          ABS(XLO-XHI)  .LE. TOL(X)
C     where
C          TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
C
C     If this algorithm does not find XLO and XHI satisfying
C     these conditions then QRZERO returns .FALSE.  This
C     implies that the input condition was not met.
C
C
C                              Arguments
C
C
C     XLO --> The left endpoint of the interval to be
C           searched for a solution.
C                    XLO is DOUBLE PRECISION
C
C     XHI --> The right endpoint of the interval to be
C           for a solution.
C                    XHI is DOUBLE PRECISION
C
C     ABSTOL, RELTOL --> Two numbers that determine the accuracy
C                      of the solution.  See function for a
C                      precise definition.
C                    ABSTOL is DOUBLE PRECISION
C                    RELTOL is DOUBLE PRECISION
C
C
C                              Method
C
C
C     Algorithm R of the paper 'Two Efficient Algorithms with
C     Guaranteed Convergence for Finding a Zero of a Function'
C     by J. C. P. Bus and T. J. Dekker in ACM Transactions on
C     Mathematical Software, Volume 1, no. 4 page 330
C     (Dec. '75) is employed to find the zero of F(X)-Y.
C
C**********************************************************************

C     Reset all saved variables to known values
      a = 0d0
      abstol = 0d0
      b = 0d0
      c = 0d0
      d = 0d0
      fa = 0d0
      fb = 0d0
      fc = 0d0
      fd = 0d0
      fda = 0d0
      fdb = 0d0
      m = 0d0
      mb = 0d0
      p = 0d0
      q = 0d0
      reltol = 0d0
      tol = 0d0
      w = 0d0
      xxhi = 0d0
      xxlo = 0d0
      zx = 0d0
      ext = 0
      i99999 = 0
      first = .FALSE.
      qrzero = .FALSE.

C     Set initial values
      xxlo = zxlo
      xxhi = zxhi
      abstol = zabstl
      reltol = zreltl
      RETURN

      STOP '*** EXECUTION FLOWING INTO FLECS PROCEDURES ***'
C     TO GET-FUNCTION-VALUE
  270 status = 1
      RETURN

  280 CONTINUE
      GO TO i99999

      END
