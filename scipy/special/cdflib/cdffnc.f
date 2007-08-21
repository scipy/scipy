      SUBROUTINE cdffnc(which,p,q,f,dfn,dfd,phonc,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFFNC( WHICH, P, Q, F, DFN, DFD, PNONC, STATUS, BOUND )
C               Cumulative Distribution Function
C               Non-central F distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the Non-central F
C     distribution given values for the others.
C
C
C                              Arguments
C
C
C     WHICH --> Integer indicating which of the next five argument
C               values is to be calculated from the others.
C               Legal range: 1..5
C               iwhich = 1 : Calculate P and Q from F,DFN,DFD and PNONC
C               iwhich = 2 : Calculate F from P,Q,DFN,DFD and PNONC
C               iwhich = 3 : Calculate DFN from P,Q,F,DFD and PNONC
C               iwhich = 4 : Calculate DFD from P,Q,F,DFN and PNONC
C               iwhich = 5 : Calculate PNONC from P,Q,F,DFN and DFD
C                    INTEGER WHICH
C
C       P <--> The integral from 0 to F of the non-central f-density.
C              Input range: [0,1-1E-16).
C                    DOUBLE PRECISION P
C
C       Q <--> 1-P.
C            Q is not used by this subroutine and is only included
C            for similarity with other cdf* routines.
C                    DOUBLE PRECISION Q
C
C       F <--> Upper limit of integration of the non-central f-density.
C              Input range: [0, +infinity).
C              Search range: [0,1E100]
C                    DOUBLE PRECISION F
C
C     DFN < --> Degrees of freedom of the numerator sum of squares.
C               Input range: (0, +infinity).
C               Search range: [ 1E-100, 1E100]
C                    DOUBLE PRECISION DFN
C
C     DFD < --> Degrees of freedom of the denominator sum of squares.
C               Must be in range: (0, +infinity).
C               Input range: (0, +infinity).
C               Search range: [ 1E-100, 1E100]
C                    DOUBLE PRECISION DFD
C
C     PNONC <-> The non-centrality parameter
C               Input range: [0,infinity)
C               Search range: [0,1E4]
C                    DOUBLE PRECISION PHONC
C
C     STATUS <-- 0 if calculation completed correctly
C               -I if input parameter number I is out of range
C                1 if answer appears to be lower than lowest
C                  search bound
C                2 if answer appears to be higher than greatest
C                  search bound
C                3 if P + Q .ne. 1
C                    INTEGER STATUS
C
C     BOUND <-- Undefined if STATUS is 0
C
C               Bound exceeded by parameter number I if STATUS
C               is negative.
C
C               Lower search bound if STATUS is 1.
C
C               Upper search bound if STATUS is 2.
C
C
C                              Method
C
C
C     Formula  26.6.20   of   Abramowitz   and   Stegun,  Handbook  of
C     Mathematical  Functions (1966) is used to compute the cumulative
C     distribution function.
C
C     Computation of other parameters involve a seach for a value that
C     produces  the desired  value  of P.   The search relies  on  the
C     monotinicity of P with the other parameter.
C
C                            WARNING
C
C     The computation time  required for this  routine is proportional
C     to the noncentrality  parameter  (PNONC).  Very large  values of
C     this parameter can consume immense  computer resources.  This is
C     why the search range is bounded by 10,000.
C
C                              WARNING
C
C     The  value  of the  cumulative  noncentral F distribution is not
C     necessarily monotone in either degrees  of freedom.  There  thus
C     may be two values that provide a given  CDF value.  This routine
C     assumes monotonicity  and will find  an arbitrary one of the two
C     values.
C
C**********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION tent4
      PARAMETER (tent4=1.0D4)
      DOUBLE PRECISION tol
      PARAMETER (tol=1.0D-8)
      DOUBLE PRECISION atol
      PARAMETER (atol=1.0D-50)
      DOUBLE PRECISION zero,one,inf
      PARAMETER (zero=1.0D-100,one=1.0D0-1.0D-16,inf=1.0D100)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,dfd,dfn,f,p,phonc,q
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ccum,cum,fx
      LOGICAL qhi,qleft
C     ..
C     .. External Subroutines ..
      EXTERNAL cumfnc,dinvr,dstinv
C     ..
      IF (.NOT. ((which.LT.1).OR. (which.GT.5))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 5.0D0
   20 status = -1
      RETURN

   30 IF (which.EQ.1) GO TO 70
      IF (.NOT. ((p.LT.0.0D0).OR. (p.GT.one))) GO TO 60
      IF (.NOT. (p.LT.0.0D0)) GO TO 40
      bound = 0.0D0
      GO TO 50

   40 bound = one
   50 status = -2
      RETURN

   60 CONTINUE
   70 IF (which.EQ.2) GO TO 90
      IF (.NOT. (f.LT.0.0D0)) GO TO 80
      bound = 0.0D0
      status = -4
      RETURN

   80 CONTINUE
   90 IF (which.EQ.3) GO TO 110
      IF (.NOT. (dfn.LE.0.0D0)) GO TO 100
      bound = 0.0D0
      status = -5
      RETURN

  100 CONTINUE
  110 IF (which.EQ.4) GO TO 130
      IF (.NOT. (dfd.LE.0.0D0)) GO TO 120
      bound = 0.0D0
      status = -6
      RETURN

  120 CONTINUE
  130 IF (which.EQ.5) GO TO 150
      IF (.NOT. (phonc.LT.0.0D0)) GO TO 140
      bound = 0.0D0
      status = -7
      RETURN

  140 CONTINUE
  150 IF ((1).EQ. (which)) THEN
          CALL cumfnc(f,dfn,dfd,phonc,p,q)
          status = 0

      ELSE IF ((2).EQ. (which)) THEN
          f = 5.0D0
          CALL dstinv(0.0D0,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,f,fx,qleft,qhi)
  160     IF (.NOT. (status.EQ.1)) GO TO 170
          CALL cumfnc(f,dfn,dfd,phonc,cum,ccum)
          fx = cum - p
          CALL dinvr(status,f,fx,qleft,qhi)
          GO TO 160

  170     IF (.NOT. (status.EQ.-1)) GO TO 200
          IF (.NOT. (qleft)) GO TO 180
          status = 1
          bound = 0.0D0
          GO TO 190

  180     status = 2
          bound = inf
  190     CONTINUE
  200     CONTINUE

      ELSE IF ((3).EQ. (which)) THEN
          dfn = 5.0D0
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,dfn,fx,qleft,qhi)
  210     IF (.NOT. (status.EQ.1)) GO TO 220
          CALL cumfnc(f,dfn,dfd,phonc,cum,ccum)
          fx = cum - p
          CALL dinvr(status,dfn,fx,qleft,qhi)
          GO TO 210

  220     IF (.NOT. (status.EQ.-1)) GO TO 250
          IF (.NOT. (qleft)) GO TO 230
          status = 1
          bound = zero
          GO TO 240

  230     status = 2
          bound = inf
  240     CONTINUE
  250     CONTINUE

      ELSE IF ((4).EQ. (which)) THEN
          dfd = 5.0D0
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,dfd,fx,qleft,qhi)
  260     IF (.NOT. (status.EQ.1)) GO TO 270
          CALL cumfnc(f,dfn,dfd,phonc,cum,ccum)
          fx = cum - p
          CALL dinvr(status,dfd,fx,qleft,qhi)
          GO TO 260

  270     IF (.NOT. (status.EQ.-1)) GO TO 300
          IF (.NOT. (qleft)) GO TO 280
          status = 1
          bound = zero
          GO TO 290

  280     status = 2
          bound = inf
  290     CONTINUE
  300     CONTINUE

      ELSE IF ((5).EQ. (which)) THEN
          phonc = 5.0D0
          CALL dstinv(0.0D0,tent4,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,phonc,fx,qleft,qhi)
  310     IF (.NOT. (status.EQ.1)) GO TO 320
          CALL cumfnc(f,dfn,dfd,phonc,cum,ccum)
          fx = cum - p
          CALL dinvr(status,phonc,fx,qleft,qhi)
          GO TO 310

  320     IF (.NOT. (status.EQ.-1)) GO TO 350
          IF (.NOT. (qleft)) GO TO 330
          status = 1
          bound = 0.0D0
          GO TO 340

  330     status = 2
          bound = tent4
  340     CONTINUE
  350 END IF

      RETURN

      END
