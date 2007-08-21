C------------------------------------------------------------------------------
      SUBROUTINE TRSTLP (N,M,A,B,RHO,DX,IFULL,IACT,Z,ZDOTA,VMULTC,
     1  SDIRN,DXNEW,VMULTD,IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION TEMP
      DIMENSION A(N,*),B(*),DX(*),IACT(*),Z(N,*),ZDOTA(*),
     1  VMULTC(*),SDIRN(*),DXNEW(*),VMULTD(*)
C
C     This subroutine calculates an N-component vector DX by applying the
C     following two stages. In the first stage, DX is set to the shortest
C     vector that minimizes the greatest violation of the constraints
C       A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K=2,3,...,M,
C     subject to the Euclidean length of DX being at most RHO. If its length is
C     strictly less than RHO, then we use the resultant freedom in DX to
C     minimize the objective function
C              -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N)
C     subject to no increase in any greatest constraint violation. This
C     notation allows the gradient of the objective function to be regarded as
C     the gradient of a constraint. Therefore the two stages are distinguished
C     by MCON .EQ. M and MCON .GT. M respectively. It is possible that a
C     degeneracy may prevent DX from attaining the target length RHO. Then the
C     value IFULL=0 would be set, but usually IFULL=1 on return.
C
C     In general NACT is the number of constraints in the active set and
C     IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT
C     contains a permutation of the remaining constraint indices. Further, Z is
C     an orthogonal matrix whose first NACT columns can be regarded as the
C     result of Gram-Schmidt applied to the active constraint gradients. For
C     J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th
C     column of Z with the gradient of the J-th active constraint. DX is the
C     current vector of variables and here the residuals of the active
C     constraints should be zero. Further, the active constraints have
C     nonnegative Lagrange multipliers that are held at the beginning of
C     VMULTC. The remainder of this vector holds the residuals of the inactive
C     constraints at DX, the ordering of the components of VMULTC being in
C     agreement with the permutation of the indices of the constraints that is
C     in IACT. All these residuals are nonnegative, which is achieved by the
C     shift RESMAX that makes the least residual zero.
C
C     Initialize Z and some other variables. The value of RESMAX will be
C     appropriate to DX=0, while ICON will be the index of a most violated
C     constraint if RESMAX is positive. Usually during the first stage the
C     vector SDIRN gives a search direction that reduces all the active
C     constraint violations by one simultaneously.
C

      IF (IPRINT .EQ. 3) THEN
         print *, ' '
         print *, 'BEFORE trstlp:'
         PRINT *, '  **DX = ', (DX(I),I=1,N)
         PRINT *, '  **IACT = ', (IACT(I),I=1,M+1)
         PRINT *, 'M,N,RHO,IFULL =', M, N, RHO, IFULL
         PRINT *, '  **A = ', ((A(I,K),I=1,N),K=1,M+1)
         PRINT *, '  **B = ', (B(I),I=1,M)
         PRINT *, '  **Z = ', ((Z(I,K),I=1,N),K=1,N)
         PRINT *, '  **ZDOTA = ', (ZDOTA(I),I=1,N)
         PRINT *, '  **VMULTC = ', (VMULTC(I),I=1,M+1)
         PRINT *, '  **SDIRN = ', (SDIRN(I),I=1,N)
         PRINT *, '  **DXNEW = ', (DXNEW(I),I=1,N)
         PRINT *, '  **VMULTD = ', (VMULTD(I),I=1,M+1)
         PRINT *, ' '
      END IF

      ICON=0
      NACTX=0
      RESOLD=0

      IFULL=1
      MCON=M
      NACT=0
      RESMAX=0.0d0
      DO 20 I=1,N
      DO 10 J=1,N
   10 Z(I,J)=0.0d0
      Z(I,I)=1.0d0
   20 DX(I)=0.0d0
      IF (M .GE. 1) THEN
          DO 30 K=1,M
          IF (B(K) .GT. RESMAX) THEN
              RESMAX=B(K)
              ICON=K
          END IF
   30     CONTINUE
          DO 40 K=1,M
          IACT(K)=K
   40     VMULTC(K)=RESMAX-B(K)
      END IF
      IF (IPRINT .EQ. 3) THEN
         PRINT *, '  1. VMULTC = ', (VMULTC(I),I=1,M+1)
      END IF
      IF (RESMAX .EQ. 0.0d0) GOTO 480
      DO 50 I=1,N
   50 SDIRN(I)=0.0d0
C
C     End the current stage of the calculation if 3 consecutive iterations
C     have either failed to reduce the best calculated value of the objective
C     function or to increase the number of active constraints since the best
C     value was calculated. This strategy prevents cycling, but there is a
C     remote possibility that it will cause premature termination.
C
   60 OPTOLD=0.0d0
      ICOUNT=0
   70 IF (MCON .EQ. M) THEN
          OPTNEW=RESMAX
      ELSE
          OPTNEW=0.0d0
          DO 80 I=1,N
   80     OPTNEW=OPTNEW-DX(I)*A(I,MCON)
      END IF
      IF (IPRINT .EQ. 3) THEN
         PRINT *, ' ICOUNT, OPTNEW, OPTOLD = ', ICOUNT, OPTNEW, OPTOLD
      END IF
      IF (ICOUNT .EQ. 0 .OR. OPTNEW .LT. OPTOLD) THEN
          OPTOLD=OPTNEW
          NACTX=NACT
          ICOUNT=3
      ELSE IF (NACT .GT. NACTX) THEN
          NACTX=NACT
          ICOUNT=3
      ELSE
          ICOUNT=ICOUNT-1
          IF (ICOUNT .EQ. 0) GOTO 490
      END IF
C
C     If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
C     the active set. Apply Givens rotations so that the last N-NACT-1 columns
C     of Z are orthogonal to the gradient of the new constraint, a scalar
C     product being set to zero if its nonzero value could be due to computer
C     rounding errors. The array DXNEW is used for working space.
C
      IF (ICON .LE. NACT) GOTO 260
      KK=IACT(ICON)
      DO 90 I=1,N
   90 DXNEW(I)=A(I,KK)
      TOT=0.0D0
      K=N
  100 IF (K .GT. NACT) THEN
          SP=0.0d0
          SPABS=0.0d0
          DO 110 I=1,N
          TEMP=Z(I,K)*DXNEW(I)
          SP=SP+TEMP
  110     SPABS=SPABS+DABS(TEMP)
          ACCA=SPABS+0.1d0*DABS(SP)
          ACCB=SPABS+0.2d0*DABS(SP)
          IF ((SPABS .GE. ACCA) .OR. (ACCA .GE. ACCB)) SP=0.0D0
          IF (TOT .EQ. 0.0D0) THEN
              TOT=SP
          ELSE
              KP=K+1
              TEMP=DSQRT(SP*SP+TOT*TOT)
              ALPHA=SP/TEMP
              BETA=TOT/TEMP
              TOT=TEMP
              DO 120 I=1,N
              TEMP=ALPHA*Z(I,K)+BETA*Z(I,KP)
              Z(I,KP)=ALPHA*Z(I,KP)-BETA*Z(I,K)
  120         Z(I,K)=TEMP
          END IF
          K=K-1
          GOTO 100
      END IF
C
C     Add the new constraint if this can be done without a deletion from the
C     active set.
C
      IF (IPRINT .EQ. 3) THEN
         PRINT *, '*TOT, NACT, ICON = ', TOT, NACT, ICON
      END IF
      IF (TOT .NE. 0.0d0) THEN
          NACT=NACT+1
          ZDOTA(NACT)=TOT
          VMULTC(ICON)=VMULTC(NACT)
          VMULTC(NACT)=0.0d0
          GOTO 210
      END IF
C
C     The next instruction is reached if a deletion has to be made from the
C     active set in order to make room for the new active constraint, because
C     the new constraint gradient is a linear combination of the gradients of
C     the old active constraints. Set the elements of VMULTD to the multipliers
C     of the linear combination. Further, set IOUT to the index of the
C     constraint to be deleted, but branch if no suitable index can be found.
C
      RATIO=-1.0d0
      K=NACT
  130 ZDOTV=0.0d0
      ZDVABS=0.0d0
      DO 140 I=1,N
      TEMP=Z(I,K)*DXNEW(I)
      ZDOTV=ZDOTV+TEMP
  140 ZDVABS=ZDVABS+DABS(TEMP)
      ACCA=ZDVABS+0.1d0*DABS(ZDOTV)
      ACCB=ZDVABS+0.2d0*DABS(ZDOTV)
      IF (ZDVABS .LT. ACCA .AND. ACCA .LT. ACCB) THEN
          TEMP=ZDOTV/ZDOTA(K)
          IF (TEMP .GT. 0.0d0 .AND. IACT(K) .LE. M) THEN
              TEMPA=VMULTC(K)/TEMP
              IF (RATIO .LT. 0.0d0 .OR. TEMPA .LT. RATIO) THEN
                  RATIO=TEMPA
                  IOUT=K
              END IF
           END IF
          IF (K .GE. 2) THEN
              KW=IACT(K)
              DO 150 I=1,N
  150         DXNEW(I)=DXNEW(I)-TEMP*A(I,KW)
          END IF
          VMULTD(K)=TEMP
      ELSE
          VMULTD(K)=0.0d0
      END IF
      K=K-1
      IF (K .GT. 0) GOTO 130
      IF (IPRINT .EQ. 3) THEN
         PRINT *, '  1. VMULTD = ', (VMULTD(I),I=1,M+1)
      END IF
      IF (RATIO .LT. 0.0d0) GOTO 490
C
C     Revise the Lagrange multipliers and reorder the active constraints so
C     that the one to be replaced is at the end of the list. Also calculate the
C     new value of ZDOTA(NACT) and branch if it is not acceptable.
C
      DO 160 K=1,NACT
  160 VMULTC(K)=DMAX1(0.0d0,VMULTC(K)-RATIO*VMULTD(K))
      IF (IPRINT .EQ. 3) THEN
         PRINT *, '  2. VMULTC = ', (VMULTC(I),I=1,M+1)
      END IF
      IF (ICON .LT. NACT) THEN
          ISAVE=IACT(ICON)
          VSAVE=VMULTC(ICON)
          K=ICON
  170     KP=K+1
          KW=IACT(KP)
          SP=0.0d0
          DO 180 I=1,N
  180     SP=SP+Z(I,K)*A(I,KW)
          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
          ALPHA=ZDOTA(KP)/TEMP
          BETA=SP/TEMP
          ZDOTA(KP)=ALPHA*ZDOTA(K)
          ZDOTA(K)=TEMP
          DO 190 I=1,N
          TEMP=ALPHA*Z(I,KP)+BETA*Z(I,K)
          Z(I,KP)=ALPHA*Z(I,K)-BETA*Z(I,KP)
  190     Z(I,K)=TEMP
          IACT(K)=KW
          VMULTC(K)=VMULTC(KP)
          K=KP
          IF (K .LT. NACT) GOTO 170
          IACT(K)=ISAVE
          VMULTC(K)=VSAVE
      END IF
      TEMP=0.0d0
      DO 200 I=1,N
  200 TEMP=TEMP+Z(I,NACT)*A(I,KK)
      IF (TEMP .EQ. 0.0d0) GOTO 490
      ZDOTA(NACT)=TEMP
      VMULTC(ICON)=0.0d0
      VMULTC(NACT)=RATIO
C
C     Update IACT and ensure that the objective function continues to be
C     treated as the last active constraint when MCON>M.
C
  210 IACT(ICON)=IACT(NACT)
      IACT(NACT)=KK
      IF (MCON .GT. M .AND. KK .NE. MCON) THEN
          K=NACT-1
          SP=0.0d0
          DO 220 I=1,N
  220     SP=SP+Z(I,K)*A(I,KK)
          TEMP=SQRT(SP*SP+ZDOTA(NACT)**2)
          ALPHA=ZDOTA(NACT)/TEMP
          BETA=SP/TEMP
          ZDOTA(NACT)=ALPHA*ZDOTA(K)
          ZDOTA(K)=TEMP
          DO 230 I=1,N
          TEMP=ALPHA*Z(I,NACT)+BETA*Z(I,K)
          Z(I,NACT)=ALPHA*Z(I,K)-BETA*Z(I,NACT)
  230     Z(I,K)=TEMP
          IACT(NACT)=IACT(K)
          IACT(K)=KK
          TEMP=VMULTC(K)
          VMULTC(K)=VMULTC(NACT)
          VMULTC(NACT)=TEMP
      END IF
C
C     If stage one is in progress, then set SDIRN to the direction of the next
C     change to the current vector of variables.
C
      IF (MCON .GT. M) GOTO 320
      KK=IACT(NACT)
      TEMP=0.0d0
      DO 240 I=1,N
  240 TEMP=TEMP+SDIRN(I)*A(I,KK)
      TEMP=TEMP-1.0d0
      TEMP=TEMP/ZDOTA(NACT)
      DO 250 I=1,N
  250 SDIRN(I)=SDIRN(I)-TEMP*Z(I,NACT)
      GOTO 340
C
C     Delete the constraint that has the index IACT(ICON) from the active set.
C
  260 IF (ICON .LT. NACT) THEN
          ISAVE=IACT(ICON)
          VSAVE=VMULTC(ICON)
          K=ICON
  270     KP=K+1
          KK=IACT(KP)
          SP=0.0d0
          DO 280 I=1,N
  280     SP=SP+Z(I,K)*A(I,KK)
          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
          ALPHA=ZDOTA(KP)/TEMP
          BETA=SP/TEMP
          ZDOTA(KP)=ALPHA*ZDOTA(K)
          ZDOTA(K)=TEMP
          DO 290 I=1,N
          TEMP=ALPHA*Z(I,KP)+BETA*Z(I,K)
          Z(I,KP)=ALPHA*Z(I,K)-BETA*Z(I,KP)
  290     Z(I,K)=TEMP
          IACT(K)=KK
          VMULTC(K)=VMULTC(KP)
          K=KP
          IF (K .LT. NACT) GOTO 270
          IACT(K)=ISAVE
          VMULTC(K)=VSAVE
      END IF
      NACT=NACT-1
C
C     If stage one is in progress, then set SDIRN to the direction of the next
C     change to the current vector of variables.
C
      IF (MCON .GT. M) GOTO 320
      TEMP=0.0d0
      DO 300 I=1,N
  300 TEMP=TEMP+SDIRN(I)*Z(I,NACT+1)
      DO 310 I=1,N
  310 SDIRN(I)=SDIRN(I)-TEMP*Z(I,NACT+1)
      GO TO 340
C
C     Pick the next search direction of stage two.
C
  320 TEMP=1.0d0/ZDOTA(NACT)
      DO 330 I=1,N
  330 SDIRN(I)=TEMP*Z(I,NACT)
C
C     Calculate the step to the boundary of the trust region or take the step
C     that reduces RESMAX to zero. The two statements below that include the
C     factor 1.0E-6 prevent some harmless underflows that occurred in a test
C     calculation. Further, we skip the step if it could be zero within a
C     reasonable tolerance for computer rounding errors.
C
  340 DD=RHO*RHO
      SD=0.0d0
      SS=0.0d0
      DO 350 I=1,N
      IF (ABS(DX(I)) .GE. 1.0E-6*RHO) DD=DD-DX(I)**2
      SD=SD+DX(I)*SDIRN(I)
  350 SS=SS+SDIRN(I)**2
      IF (DD .LE. 0.0d0) GOTO 490
      TEMP=SQRT(SS*DD)
      IF (ABS(SD) .GE. 1.0E-6*TEMP) TEMP=SQRT(SS*DD+SD*SD)
      STPFUL=DD/(TEMP+SD)
      STEP=STPFUL
      IF (MCON .EQ. M) THEN
          ACCA=STEP+0.1d0*RESMAX
          ACCB=STEP+0.2d0*RESMAX
          IF (STEP .GE. ACCA .OR. ACCA .GE. ACCB) GOTO 480
          STEP=DMIN1(STEP,RESMAX)
      END IF
C
C     Set DXNEW to the new variables if STEP is the steplength, and reduce
C     RESMAX to the corresponding maximum residual if stage one is being done.
C     Because DXNEW will be changed during the calculation of some Lagrange
C     multipliers, it will be restored to the following value later.
      call s360_380(DXNEW,DX,STEP,SDIRN,N,M,MCON,RESMAX,
     1 NACT,IACT,B,A,RESOLD)

C
C     Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
C     device is included to force VMULTD(K)=0.0 if deviations from this value
C     can be attributed to computer rounding errors. First calculate the new
C     Lagrange multipliers.
C
      K=NACT
  390 ZDOTW=0.0d0
      ZDWABS=0.0d0
      DO 400 I=1,N
      TEMP=Z(I,K)*DXNEW(I)
      ZDOTW=ZDOTW+TEMP
  400 ZDWABS=ZDWABS+ABS(TEMP)
      ACCA=ZDWABS+0.1d0*ABS(ZDOTW)
      ACCB=ZDWABS+0.2d0*ABS(ZDOTW)
      IF (ZDWABS .GE. ACCA .OR. ACCA .GE. ACCB) ZDOTW=0.0d0
      VMULTD(K)=ZDOTW/ZDOTA(K)
      IF (K .GE. 2) THEN
          KK=IACT(K)
          DO 410 I=1,N
  410     DXNEW(I)=DXNEW(I)-VMULTD(K)*A(I,KK)
          K=K-1
          GOTO 390
      END IF
      IF (MCON .GT. M) VMULTD(NACT)=DMAX1(0.0d0,VMULTD(NACT))
      IF (IPRINT .EQ. 3) THEN     
         PRINT *, '  2. VMULTD = ', (VMULTD(I),I=1,M+1)
      END IF
C
C     Complete VMULTC by finding the new constraint residuals.
C
      DO 420 I=1,N
  420 DXNEW(I)=DX(I)+STEP*SDIRN(I)
      IF (MCON .GT. NACT) THEN
          KL=NACT+1
          DO 440 K=KL,MCON
          KK=IACT(K)
          SUM=RESMAX-B(KK)
          SUMABS=RESMAX+DABS(B(KK))
          DO 430 I=1,N
          TEMP=A(I,KK)*DXNEW(I)
          SUM=SUM+TEMP
  430     SUMABS=SUMABS+DABS(TEMP)
          ACCA=SUMABS+0.1*DABS(SUM)
          ACCB=SUMABS+0.2*DABS(SUM)
          IF (SUMABS .GE. ACCA .OR. ACCA .GE. ACCB) SUM=0.0
  440     VMULTD(K)=SUM
      END IF
      IF (IPRINT .EQ. 3) THEN
         PRINT *, '  3. VMULTD = ', (VMULTD(I),I=1,M+1)
      END IF
C
C     Calculate the fraction of the step from DX to DXNEW that will be taken.
C
      RATIO=1.0d0
      ICON=0
C     
      EPS = 2.2E-16
      DO 450 K=1,MCON
      IF (VMULTD(K) .GT. -EPS .AND. VMULTD(K) .LT. EPS) VMULTD(K)=0.0D0
      IF (VMULTD(K) .LT. 0.0D0) THEN
          TEMP=VMULTC(K)/(VMULTC(K)-VMULTD(K))
          IF (TEMP .LT. RATIO) THEN
              RATIO=TEMP
              ICON=K
          END IF
      END IF
  450 CONTINUE
C
C     Update DX, VMULTC and RESMAX.
C
      TEMP=1.0d0-RATIO
      DO 460 I=1,N
  460 DX(I)=TEMP*DX(I)+RATIO*DXNEW(I)
      DO 470 K=1,MCON
  470 VMULTC(K)=DMAX1(0.0d0,TEMP*VMULTC(K)+RATIO*VMULTD(K))
      IF (IPRINT .EQ. 3) THEN
         PRINT *, '  3. VMULTC = ', (VMULTC(I),I=1,M+1)
      END IF
      IF (MCON .EQ. M) RESMAX=RESOLD+RATIO*(RESMAX-RESOLD)
      IF (IPRINT .EQ. 3) THEN
         PRINT *, ' RESMAX, MCON, M, ICON = ', 
     1        RESMAX, MCON, M, ICON
      END IF
C
C     If the full step is not acceptable then begin another iteration.
C     Otherwise switch to stage two or end the calculation.
C
      IF (ICON .GT. 0) GOTO 70
      IF (STEP .EQ. STPFUL) GOTO 500
  480 MCON=M+1
      ICON=MCON
      IACT(MCON)=MCON
      VMULTC(MCON)=0.0d0
      GOTO 60
C
C     We employ any freedom that may be available to reduce the objective
C     function before returning a DX whose length is less than RHO.
C
  490 IF (MCON .EQ. M) GOTO 480
      IFULL=0
  500 CONTINUE
      IF (IPRINT .EQ. 3) THEN
         print *, ' '
         print *, 'AFTER trstlp:'
         PRINT *, '  **DX = ', (DX(I),I=1,N)
         PRINT *, '  **IACT = ', (IACT(I),I=1,M+1)
         PRINT *, 'M,N,RHO,IFULL =', M, N, RHO, IFULL
         PRINT *, '  **A = ', ((A(I,K),I=1,N),K=1,M+1)
         PRINT *, '  **B = ', (B(I),I=1,M)
         PRINT *, '  **Z = ', ((Z(I,K),I=1,N),K=1,N)
         PRINT *, '  **ZDOTA = ', (ZDOTA(I),I=1,N)
         PRINT *, '  **VMULTC = ', (VMULTC(I),I=1,M+1)
         PRINT *, '  **SDIRN = ', (SDIRN(I),I=1,N)
         PRINT *, '  **DXNEW = ', (DXNEW(I),I=1,N)
         PRINT *, '  **VMULTD = ', (VMULTD(I),I=1,M+1)
         PRINT *, ' '
      END IF
C  500 RETURN
      END

      subroutine s360_380(DXNEW,DX,STEP,SDIRN,N,M,MCON,RESMAX,
     1 NACT,IACT,B,A,RESOLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,*),B(*),DX(*),IACT(*), SDIRN(*),DXNEW(*)
      DO 360 I=1,N
  360 DXNEW(I)=DX(I)+STEP*SDIRN(I)
      IF (MCON .EQ. M) THEN
          RESOLD=RESMAX
          RESMAX=0.0d0
          DO 380 K=1,NACT
          KK=IACT(K)
          TEMP=B(KK)
          DO 370 I=1,N
  370     TEMP=TEMP-A(I,KK)*DXNEW(I)
          RESMAX=DMAX1(RESMAX,TEMP)
  380     CONTINUE
      END IF
      end
