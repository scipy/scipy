      SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, INFO )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      LOGICAL            WANTT, WANTZ
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
*     ..
*     .. Array Arguments ..
      COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CLAHQR is an auxiliary routine called by CHSEQR to update the
*  eigenvalues and Schur decomposition already computed by CHSEQR, by
*  dealing with the Hessenberg submatrix in rows and columns ILO to IHI.
*
*  Arguments
*  =========
*
*  WANTT   (input) LOGICAL
*          = .TRUE. : the full Schur form T is required;
*          = .FALSE.: only eigenvalues are required.
*
*  WANTZ   (input) LOGICAL
*          = .TRUE. : the matrix of Schur vectors Z is required;
*          = .FALSE.: Schur vectors are not required.
*
*  N       (input) INTEGER
*          The order of the matrix H.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that H is already upper triangular in rows and
*          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
*          CLAHQR works primarily with the Hessenberg submatrix in rows
*          and columns ILO to IHI, but applies transformations to all of
*          H if WANTT is .TRUE..
*          1 <= ILO <= max(1,IHI); IHI <= N.
*
*  H       (input/output) COMPLEX array, dimension (LDH,N)
*          On entry, the upper Hessenberg matrix H.
*          On exit, if WANTT is .TRUE., H is upper triangular in rows
*          and columns ILO:IHI, with any 2-by-2 diagonal blocks in
*          standard form. If WANTT is .FALSE., the contents of H are
*          unspecified on exit.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H. LDH >= max(1,N).
*
*  W       (output) COMPLEX array, dimension (N)
*          The computed eigenvalues ILO to IHI are stored in the
*          corresponding elements of W. If WANTT is .TRUE., the
*          eigenvalues are stored in the same order as on the diagonal
*          of the Schur form returned in H, with W(i) = H(i,i).
*
*  ILOZ    (input) INTEGER
*  IHIZ    (input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE..
*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*
*  Z       (input/output) COMPLEX array, dimension (LDZ,N)
*          If WANTZ is .TRUE., on entry Z must contain the current
*          matrix Z of transformations accumulated by CHSEQR, and on
*          exit Z has been updated; transformations are applied only to
*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*          If WANTZ is .FALSE., Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          > 0: if INFO = i, CLAHQR failed to compute all the
*               eigenvalues ILO to IHI in a total of 30*(IHI-ILO+1)
*               iterations; elements i+1:ihi of W contain those
*               eigenvalues which have been successfully computed.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ),
     $                   ONE = ( 1.0E+0, 0.0E+0 ) )
      REAL               RZERO, HALF
      PARAMETER          ( RZERO = 0.0E+0, HALF = 0.5E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I1, I2, ITN, ITS, J, K, L, M, NH, NZ
      REAL               H10, H21, RTEMP, S, SMLNUM, T2, TST1, ULP
      COMPLEX            CDUM, H11, H11S, H22, SUM, T, T1, TEMP, U, V2,
     $                   X, Y
*     ..
*     .. Local Arrays ..
      REAL               RWORK( 1 )
      COMPLEX            V( 2 )
*     ..
*     .. External Functions ..
      REAL               CLANHS, SLAMCH
      COMPLEX            WCLADIV
      EXTERNAL           CLANHS, SLAMCH, WCLADIV
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CLARFG, CSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CONJG, MAX, MIN, REAL, SQRT
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( ILO.EQ.IHI ) THEN
         W( ILO ) = H( ILO, ILO )
         RETURN
      END IF
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
*
*     Set machine-dependent constants for the stopping criterion.
*     If norm(H) <= sqrt(OVFL), overflow should not occur.
*
      ULP = SLAMCH( 'Precision' )
      SMLNUM = SLAMCH( 'Safe minimum' ) / ULP
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
*
*     ITN is the total number of QR iterations allowed.
*
      ITN = 30*NH
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of 1. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
*     H(L,L-1) is negligible so that the matrix splits.
*
      I = IHI
   10 CONTINUE
      IF( I.LT.ILO )
     $   GO TO 130
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      L = ILO
      DO 110 ITS = 0, ITN
*
*        Look for a single small subdiagonal element.
*
         DO 20 K = I, L + 1, -1
            TST1 = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
            IF( TST1.EQ.RZERO )
     $         TST1 = CLANHS( '1', I-L+1, H( L, L ), LDH, RWORK )
            IF( ABS( REAL( H( K, K-1 ) ) ).LE.MAX( ULP*TST1, SMLNUM ) )
     $         GO TO 30
   20    CONTINUE
   30    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
*
*           H(L,L-1) is negligible
*
            H( L, L-1 ) = ZERO
         END IF
*
*        Exit from loop if a submatrix of order 1 has split off.
*
         IF( L.GE.I )
     $      GO TO 120
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
*
         IF( ITS.EQ.10 .OR. ITS.EQ.20 ) THEN
*
*           Exceptional shift.
*
            T = ABS( REAL( H( I, I-1 ) ) ) +
     $          ABS( REAL( H( I-1, I-2 ) ) )
         ELSE
*
*           Wilkinson's shift.
*
            T = H( I, I )
            U = H( I-1, I )*REAL( H( I, I-1 ) )
            IF( U.NE.ZERO ) THEN
               X = HALF*( H( I-1, I-1 )-T )
               Y = SQRT( X*X+U )
               IF( REAL( X )*REAL( Y )+AIMAG( X )*AIMAG( Y ).LT.RZERO )
     $            Y = -Y
               T = T - WCLADIV( U, ( X+Y ) )
            END IF
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 40 M = I - 1, L + 1, -1
*
*           Determine the effect of starting the single-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.
*
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H11S = H11 - T
            H21 = H( M+1, M )
            S = CABS1( H11S ) + ABS( H21 )
            H11S = H11S / S
            H21 = H21 / S
            V( 1 ) = H11S
            V( 2 ) = H21
            H10 = H( M, M-1 )
            TST1 = CABS1( H11S )*( CABS1( H11 )+CABS1( H22 ) )
            IF( ABS( H10*H21 ).LE.ULP*TST1 )
     $         GO TO 50
   40    CONTINUE
         H11 = H( L, L )
         H22 = H( L+1, L+1 )
         H11S = H11 - T
         H21 = H( L+1, L )
         S = CABS1( H11S ) + ABS( H21 )
         H11S = H11S / S
         H21 = H21 / S
         V( 1 ) = H11S
         V( 2 ) = H21
   50    CONTINUE
*
*        Single-shift QR step
*
         DO 100 K = M, I - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix.
*
*           V(2) is always real before the call to CLARFG, and hence
*           after the call T2 ( = T1*V(2) ) is also real.
*
            IF( K.GT.M )
     $         CALL CCOPY( 2, H( K, K-1 ), 1, V, 1 )
            CALL CLARFG( 2, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
            END IF
            V2 = V( 2 )
            T2 = REAL( T1*V2 )
*
*           Apply G from the left to transform the rows of the matrix
*           in columns K to I2.
*
            DO 60 J = K, I2
               SUM = CONJG( T1 )*H( K, J ) + T2*H( K+1, J )
               H( K, J ) = H( K, J ) - SUM
               H( K+1, J ) = H( K+1, J ) - SUM*V2
   60       CONTINUE
*
*           Apply G from the right to transform the columns of the
*           matrix in rows I1 to min(K+2,I).
*
            DO 70 J = I1, MIN( K+2, I )
               SUM = T1*H( J, K ) + T2*H( J, K+1 )
               H( J, K ) = H( J, K ) - SUM
               H( J, K+1 ) = H( J, K+1 ) - SUM*CONJG( V2 )
   70       CONTINUE
*
            IF( WANTZ ) THEN
*
*              Accumulate transformations in the matrix Z
*
               DO 80 J = ILOZ, IHIZ
                  SUM = T1*Z( J, K ) + T2*Z( J, K+1 )
                  Z( J, K ) = Z( J, K ) - SUM
                  Z( J, K+1 ) = Z( J, K+1 ) - SUM*CONJG( V2 )
   80          CONTINUE
            END IF
*
            IF( K.EQ.M .AND. M.GT.L ) THEN
*
*              If the QR step was started at row M > L because two
*              consecutive small subdiagonals were found, then extra
*              scaling must be performed to ensure that H(M,M-1) remains
*              real.
*
               TEMP = ONE - T1
               TEMP = TEMP / ABS( TEMP )
               H( M+1, M ) = H( M+1, M )*CONJG( TEMP )
               IF( M+2.LE.I )
     $            H( M+2, M+1 ) = H( M+2, M+1 )*TEMP
               DO 90 J = M, I
                  IF( J.NE.M+1 ) THEN
                     IF( I2.GT.J )
     $                  CALL CSCAL( I2-J, TEMP, H( J, J+1 ), LDH )
                     CALL CSCAL( J-I1, CONJG( TEMP ), H( I1, J ), 1 )
                     IF( WANTZ ) THEN
                        CALL CSCAL( NZ, CONJG( TEMP ), Z( ILOZ, J ), 1 )
                     END IF
                  END IF
   90          CONTINUE
            END IF
  100    CONTINUE
*
*        Ensure that H(I,I-1) is real.
*
         TEMP = H( I, I-1 )
         IF( AIMAG( TEMP ).NE.RZERO ) THEN
            RTEMP = ABS( TEMP )
            H( I, I-1 ) = RTEMP
            TEMP = TEMP / RTEMP
            IF( I2.GT.I )
     $         CALL CSCAL( I2-I, CONJG( TEMP ), H( I, I+1 ), LDH )
            CALL CSCAL( I-I1, TEMP, H( I1, I ), 1 )
            IF( WANTZ ) THEN
               CALL CSCAL( NZ, TEMP, Z( ILOZ, I ), 1 )
            END IF
         END IF
*
  110 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      INFO = I
      RETURN
*
  120 CONTINUE
*
*     H(I,I-1) is negligible: one eigenvalue has converged.
*
      W( I ) = H( I, I )
*
*     Decrement number of remaining iterations, and return to start of
*     the main loop with new value of I.
*
      ITN = ITN - ITS
      I = L - 1
      GO TO 10
*
  130 CONTINUE
      RETURN
*
*     End of CLAHQR
*
      END



