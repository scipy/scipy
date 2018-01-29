      SUBROUTINE dinvr(status,x,fx,qleft,qhi)
C**********************************************************************
C
C     SUBROUTINE DINVR(STATUS, X, FX, QLEFT, QHI)
C          Double precision
C          bounds the zero of the function and invokes zror
C                    Reverse Communication
C
C
C                              Function
C
C
C     Bounds the    function  and  invokes  ZROR   to perform the   zero
C     finding.  STINVR  must  have   been  called  before this   routine
C     in order to set its parameters.
C
C
C                              Arguments
C
C
C     STATUS <--> At the beginning of a zero finding problem, STATUS
C                 should be set to 0 and INVR invoked.  (The value
C                 of parameters other than X will be ignored on this cal
C
C                 When INVR needs the function evaluated, it will set
C                 STATUS to 1 and return.  The value of the function
C                 should be set in FX and INVR again called without
C                 changing any of its other parameters.
C
C                 When INVR has finished without error, it will return
C                 with STATUS 0.  In that case X is approximately a root
C                 of F(X).
C
C                 If INVR cannot bound the function, it returns status
C                 -1 and sets QLEFT and QHI.
C                         INTEGER STATUS
C
C     X <-- The value of X at which F(X) is to be evaluated.
C                         DOUBLE PRECISION X
C
C     FX --> The value of F(X) calculated when INVR returns with
C            STATUS = 1.
C                         DOUBLE PRECISION FX
C
C     QLEFT <-- Defined only if QMFINV returns .FALSE.  In that
C          case it is .TRUE. If the stepping search terminated
C          unsuccessfully at SMALL.  If it is .FALSE. the search
C          terminated unsuccessfully at BIG.
C                    QLEFT is LOGICAL
C
C     QHI <-- Defined only if QMFINV returns .FALSE.  In that
C          case it is .TRUE. if F(X) .GT. Y at the termination
C          of the search and .FALSE. if F(X) .LT. Y at the
C          termination of the search.
C                    QHI is LOGICAL

C
C**********************************************************************
      IMPLICIT NONE
C     .. Scalar Arguments ..
      DOUBLE PRECISION fx,x,zabsst,zabsto,zbig,zrelst,zrelto,zsmall,
     +                 zstpmu
      INTEGER status
      LOGICAL qhi,qleft
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION absstp,abstol,big,fbig,fsmall,relstp,reltol,
     +                 small,step,stpmul,xhi,xlb,xlo,xsave,xub,yy,zx,zy,
     +                 zz
      INTEGER i99999
      LOGICAL qbdd,qcond,qdum1,qdum2,qincr,qlim,qok,qup
C     ..
C     .. External Subroutines ..
      EXTERNAL dstzr,dzror
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,min
C     ..
C     .. Statement Functions ..
      LOGICAL qxmon
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Statement Function definitions ..
      qxmon(zx,zy,zz) = zx .LE. zy .AND. zy .LE. zz
C     ..
C     .. Executable Statements ..

      IF (status.GT.0) GO TO 310

      qcond = .NOT. qxmon(small,x,big)
      IF (qcond) STOP ' SMALL, X, BIG not monotone in INVR'
      xsave = x
C
C     See that SMALL and BIG bound the zero and set QINCR
C
      x = small
C     GET-FUNCTION-VALUE
      ASSIGN 10 TO i99999
      GO TO 300

   10 fsmall = fx
      x = big
C     GET-FUNCTION-VALUE
      ASSIGN 20 TO i99999
      GO TO 300

   20 fbig = fx
      qincr = fbig .GT. fsmall
      IF (.NOT. (qincr)) GO TO 50
      IF (fsmall.LE.0.0D0) GO TO 30
      status = -1
      qleft = .TRUE.
      qhi = .TRUE.
      RETURN

   30 IF (fbig.GE.0.0D0) GO TO 40
      status = -1
      qleft = .FALSE.
      qhi = .FALSE.
      RETURN

   40 GO TO 80

   50 IF (fsmall.GE.0.0D0) GO TO 60
      status = -1
      qleft = .TRUE.
      qhi = .FALSE.
      RETURN

   60 IF (fbig.LE.0.0D0) GO TO 70
      status = -1
      qleft = .FALSE.
      qhi = .TRUE.
      RETURN

   70 CONTINUE
   80 x = xsave
      step = max(absstp,relstp*abs(x))
C      YY = F(X) - Y
C     GET-FUNCTION-VALUE
      ASSIGN 90 TO i99999
      GO TO 300

   90 yy = fx
      IF (.NOT. (yy.EQ.0.0D0)) GO TO 100
      status = 0
      qok = .TRUE.
      RETURN

  100 qup = (qincr .AND. (yy.LT.0.0D0)) .OR.
     +      (.NOT.qincr .AND. (yy.GT.0.0D0))
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     HANDLE CASE IN WHICH WE MUST STEP HIGHER
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (.NOT. (qup)) GO TO 170
      xlb = xsave
      xub = min(xlb+step,big)
      GO TO 120

  110 IF (qcond) GO TO 150
C      YY = F(XUB) - Y
  120 x = xub
C     GET-FUNCTION-VALUE
      ASSIGN 130 TO i99999
      GO TO 300

  130 yy = fx
      qbdd = (qincr .AND. (yy.GE.0.0D0)) .OR.
     +       (.NOT.qincr .AND. (yy.LE.0.0D0))
      qlim = xub .GE. big
      qcond = qbdd .OR. qlim
      IF (qcond) GO TO 140
      step = stpmul*step
      xlb = xub
      xub = min(xlb+step,big)
  140 GO TO 110

  150 IF (.NOT. (qlim.AND..NOT.qbdd)) GO TO 160
      status = -1
      qleft = .FALSE.
      qhi = .NOT. qincr
      x = big
      RETURN

  160 GO TO 240
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     HANDLE CASE IN WHICH WE MUST STEP LOWER
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  170 xub = xsave
      xlb = max(xub-step,small)
      GO TO 190

  180 IF (qcond) GO TO 220
C      YY = F(XLB) - Y
  190 x = xlb
C     GET-FUNCTION-VALUE
      ASSIGN 200 TO i99999
      GO TO 300

  200 yy = fx
      qbdd = (qincr .AND. (yy.LE.0.0D0)) .OR.
     +       (.NOT.qincr .AND. (yy.GE.0.0D0))
      qlim = xlb .LE. small
      qcond = qbdd .OR. qlim
      IF (qcond) GO TO 210
      step = stpmul*step
      xub = xlb
      xlb = max(xub-step,small)
  210 GO TO 180

  220 IF (.NOT. (qlim.AND..NOT.qbdd)) GO TO 230
      status = -1
      qleft = .TRUE.
      qhi = qincr
      x = small
      RETURN

  230 CONTINUE
  240 CALL dstzr(xlb,xub,abstol,reltol)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      status = 0
      GO TO 260

  250 IF (.NOT. (status.EQ.1)) GO TO 290
  260 CALL dzror(status,x,fx,xlo,xhi,qdum1,qdum2)
      IF (.NOT. (status.EQ.1)) GO TO 280
C     GET-FUNCTION-VALUE
      ASSIGN 270 TO i99999
      GO TO 300

  270 CONTINUE
  280 GO TO 250

  290 x = xlo
      status = 0
      RETURN

      ENTRY dstinv(zsmall,zbig,zabsst,zrelst,zstpmu,zabsto,zrelto)
C**********************************************************************
C
C      SUBROUTINE DSTINV( SMALL, BIG, ABSSTP, RELSTP, STPMUL,
C     +                   ABSTOL, RELTOL )
C      Double Precision - SeT INverse finder - Reverse Communication
C
C
C                              Function
C
C
C     Concise Description - Given a monotone function F finds X
C     such that F(X) = Y.  Uses Reverse communication -- see invr.
C     This routine sets quantities needed by INVR.
C
C          More Precise Description of INVR -
C
C     F must be a monotone function, the results of QMFINV are
C     otherwise undefined.  QINCR must be .TRUE. if F is non-
C     decreasing and .FALSE. if F is non-increasing.
C
C     QMFINV will return .TRUE. if and only if F(SMALL) and
C     F(BIG) bracket Y, i. e.,
C          QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or
C          QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL)
C
C     if QMFINV returns .TRUE., then the X returned satisfies
C     the following condition.  let
C               TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
C     then if QINCR is .TRUE.,
C          F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X))
C     and if QINCR is .FALSE.
C          F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
C
C
C                              Arguments
C
C
C     SMALL --> The left endpoint of the interval to be
C          searched for a solution.
C                    SMALL is DOUBLE PRECISION
C
C     BIG --> The right endpoint of the interval to be
C          searched for a solution.
C                    BIG is DOUBLE PRECISION
C
C     ABSSTP, RELSTP --> The initial step size in the search
C          is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
C                    ABSSTP is DOUBLE PRECISION
C                    RELSTP is DOUBLE PRECISION
C
C     STPMUL --> When a step doesn't bound the zero, the step
C                size is multiplied by STPMUL and another step
C                taken.  A popular value is 2.0
C                    DOUBLE PRECISION STPMUL
C
C     ABSTOL, RELTOL --> Two numbers that determine the accuracy
C          of the solution.  See function for a precise definition.
C                    ABSTOL is DOUBLE PRECISION
C                    RELTOL is DOUBLE PRECISION
C
C
C                              Method
C
C
C     Compares F(X) with Y for the input value of X then uses QINCR
C     to determine whether to step left or right to bound the
C     desired x.  the initial step size is
C          MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
C     Iteratively steps right or left until it bounds X.
C     At each step which doesn't bound X, the step size is doubled.
C     The routine is careful never to step beyond SMALL or BIG.  If
C     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
C     after setting QLEFT and QHI.
C
C     If X is successfully bounded then Algorithm R of the paper
C     'Two Efficient Algorithms with Guaranteed Convergence for
C     Finding a Zero of a Function' by J. C. P. Bus and
C     T. J. Dekker in ACM Transactions on Mathematical
C     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
C     to find the zero of the function F(X)-Y. This is routine
C     QRZERO.
C
C**********************************************************************

C     Reset all saved variables to known values
      absstp = 0d0
      abstol = 0d0
      big = 0d0
      fbig = 0d0
      fsmall = 0d0
      relstp = 0d0
      reltol = 0d0
      small = 0d0
      step = 0d0
      stpmul = 0d0
      xhi = 0d0
      xlb = 0d0
      xlo = 0d0
      xsave = 0d0
      xub = 0d0
      yy = 0d0
      zx = 0d0
      zy = 0d0
      zz = 0d0
      i99999 = 0
      qbdd = .FALSE.
      qcond = .FALSE.
      qdum1 = .FALSE.
      qdum2 = .FALSE.
      qincr = .FALSE.
      qlim = .FALSE.
      qok = .FALSE.
      qup = .FALSE.

C     Set initial values
      small = zsmall
      big = zbig
      absstp = zabsst
      relstp = zrelst
      stpmul = zstpmu
      abstol = zabsto
      reltol = zrelto
      RETURN

      STOP '*** EXECUTION FLOWING INTO FLECS PROCEDURES ***'
C     TO GET-FUNCTION-VALUE
  300 status = 1
      RETURN

  310 CONTINUE
      GO TO i99999

      END
