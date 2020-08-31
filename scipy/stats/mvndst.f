* Note: The test program has been removed and a utility routine mvnun has been
* added.  RTK 2004-08-10
*
* Copyright 2000 by Alan Genz.
* Copyright 2004-2005 by Enthought, Inc.
*
* The subroutine MVNUN is copyrighted by Enthought, Inc.
* The rest of the file is copyrighted by Alan Genz and has kindly been offered
* to the Scipy project under it's BSD-style license.
*
* This file contains a short test program and MVNDST, a subroutine
* for computing multivariate normal distribution function values.
* The file is self contained and should compile without errors on (77) 
* standard Fortran compilers. The test program demonstrates the use of
* MVNDST for computing MVN distribution values for a five dimensional 
* example problem, with three different integration limit combinations.
*
*          Alan Genz
*          Department of Mathematics
*          Washington State University
*          Pullman, WA 99164-3113
*          Email : alangenz@wsu.edu
*
      SUBROUTINE mvnun(d, n, lower, upper, means, covar, maxpts, 
     &                   abseps, releps, value, inform)
*  Parameters
*
*   d       integer, dimensionality of the data
*   n       integer, the number of data points
*   lower   double(2), the lower integration limits
*   upper   double(2), the upper integration limits
*   means   double(n), the mean of each kernel
*   covar   double(2,2), the covariance matrix
*   maxpts  integer, the maximum number of points to evaluate at
*   abseps  double, absolute error tolerance
*   releps  double, relative error tolerance
*   value   double intent(out), integral value
*   inform  integer intent(out), 
*               if inform == 0: error < eps
*               elif inform == 1: error > eps, all maxpts used
      integer n, d, infin(d), maxpts, inform, tmpinf
      double precision lower(d), upper(d), releps, abseps,
     &                 error, value, stdev(d), rho(d*(d-1)/2), 
     &                 covar(d,d),
     &                 nlower(d), nupper(d), means(d,n), tmpval,
     &                 inf
      integer i, j

      inf = 0d0
      inf = 1d0 / inf

      do i=1,d
        stdev(i) = dsqrt(covar(i,i))
        if (upper(i).eq.inf.and.lower(i).eq.-inf) then
           infin(i) = -1
        else if (lower(i).eq.-inf) then
           infin(i) = 0
        else if (upper(i).eq.inf) then
           infin(i) = 1
        else
           infin(i) = 2
        end if
      end do
      do i=1,d
        do j=1,i-1
          rho(j+(i-2)*(i-1)/2) = covar(i,j)/stdev(i)/stdev(j)
        end do
      end do
      value = 0d0

      inform = 0

      do i=1,n
        do j=1,d
          nlower(j) = (lower(j) - means(j,i))/stdev(j)
          nupper(j) = (upper(j) - means(j,i))/stdev(j)
        end do
        call mvndst(d,nlower,nupper,infin,rho,maxpts,abseps,releps,
     &              error,tmpval,tmpinf)
        value = value + tmpval
        if (tmpinf .eq. 1) then
            inform = 1
        end if
      end do

      value = value / n
      
      END 


      SUBROUTINE mvnun_weighted(d, n, lower, upper, means, weights,
     &                          covar, maxpts, abseps, releps, 
     &                           value, inform)
*  Parameters
*
*   d       integer, dimensionality of the data
*   n       integer, the number of data points
*   lower   double(2), the lower integration limits
*   upper   double(2), the upper integration limits
*   means   double(n), the mean of each kernel
*   weights double(n), the weight of each kernel
*   covar   double(2,2), the covariance matrix
*   maxpts  integer, the maximum number of points to evaluate at
*   abseps  double, absolute error tolerance
*   releps  double, relative error tolerance
*   value   double intent(out), integral value
*   inform  integer intent(out), 
*               if inform == 0: error < eps
*               elif inform == 1: error > eps, all maxpts used
      integer n, d, infin(d), maxpts, inform, tmpinf
      double precision lower(d), upper(d), releps, abseps,
     &                 error, value, stdev(d), rho(d*(d-1)/2), 
     &                 covar(d,d),
     &                 nlower(d), nupper(d), means(d,n), tmpval,
     &                 inf, weights(n)
      integer i, j

      inf = 0d0
      inf = 1d0 / inf

      do i=1,d
        stdev(i) = dsqrt(covar(i,i))
        if (upper(i).eq.inf.and.lower(i).eq.-inf) then
           infin(i) = -1
        else if (lower(i).eq.-inf) then
           infin(i) = 0
        else if (upper(i).eq.inf) then
           infin(i) = 1
        else
           infin(i) = 2
        end if
      end do
      do i=1,d
        do j=1,i-1
          rho(j+(i-2)*(i-1)/2) = covar(i,j)/stdev(i)/stdev(j)
        end do
      end do
      value = 0d0

      inform = 0

      do i=1,n
        do j=1,d
          nlower(j) = (lower(j) - means(j,i))/stdev(j)
          nupper(j) = (upper(j) - means(j,i))/stdev(j)
        end do
        call mvndst(d,nlower,nupper,infin,rho,maxpts,abseps,releps,
     &              error,tmpval,tmpinf)
        value = value + tmpval * weights(i)
        if (tmpinf .eq. 1) then
            inform = 1
        end if
      end do
      
      END 

      SUBROUTINE MVNDST( N, LOWER, UPPER, INFIN, CORREL, MAXPTS,
     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM )
*
*     A subroutine for computing multivariate normal probabilities.
*     This subroutine uses an algorithm given in the paper
*     "Numerical Computation of Multivariate Normal Probabilities", in
*     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
*          Alan Genz 
*          Department of Mathematics
*          Washington State University 
*          Pullman, WA 99164-3113
*          Email : AlanGenz@wsu.edu
*
*  Parameters
*
*     N      INTEGER, the number of variables.
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, array of correlation coefficients; the correlation
*            coefficient in row I column J of the correlation matrix
*            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
*            THe correlation matrix must be positive semidefinite.
*     MAXPTS INTEGER, maximum number of function values allowed. This 
*            parameter can be used to limit the time. A sensible 
*            strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function values used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 500 or N < 1.
*
      EXTERNAL MVNDFN
      INTEGER N, INFIN(*), MAXPTS, INFORM, INFIS, IVLS
      DOUBLE PRECISION CORREL(*), LOWER(*), UPPER(*), RELEPS, ABSEPS,
     &                 ERROR, VALUE, E, D, MVNDNT, MVNDFN
      COMMON /DKBLCK/IVLS
      IF ( N .GT. 500 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
      ELSE
         INFORM = MVNDNT(N, CORREL, LOWER, UPPER, INFIN, INFIS, D, E)
         IF ( N-INFIS .EQ. 0 ) THEN
            VALUE = 1
            ERROR = 0
         ELSE IF ( N-INFIS .EQ. 1 ) THEN
            VALUE = E - D
            ERROR = 2D-16
         ELSE
*
*        Call the lattice rule integration subroutine
*
            IVLS = 0
            CALL DKBVRC( N-INFIS-1, IVLS, MAXPTS, MVNDFN, 
     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM )
         ENDIF
      ENDIF
      END
      DOUBLE PRECISION FUNCTION MVNDFN( N, W )
*     
*     Integrand subroutine
*
      INTEGER N, INFIN(*), INFIS, NL
      DOUBLE PRECISION W(*), LOWER(*), UPPER(*), CORREL(*), D, E
      PARAMETER ( NL = 500 )
      DOUBLE PRECISION COV(NL*(NL+1)/2), A(NL), B(NL), Y(NL)
      INTEGER INFI(NL), I, J, IJ, IK, INFA, INFB
      DOUBLE PRECISION SUM, AI, BI, DI, EI, PHINVS, BVNMVN, MVNDNT
      SAVE A, B, INFI, COV
      MVNDFN = 1
      INFA = 0
      INFB = 0
      IK = 1
      IJ = 0
      DO I = 1, N+1
         SUM = 0
         DO J = 1, I-1
            IJ = IJ + 1
            IF ( J .LT. IK ) SUM = SUM + COV(IJ)*Y(J)
         END DO
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( INFA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) - SUM )
            ELSE
               AI = A(I) - SUM 
               INFA = 1
            END IF
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( INFB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) - SUM )
            ELSE
               BI = B(I) - SUM 
               INFB = 1
            END IF
         END IF
         IJ = IJ + 1
         IF ( I .EQ. N+1 .OR. COV(IJ+IK+1) .GT. 0 ) THEN 
            CALL MVNLMS( AI, BI, 2*INFA+INFB-1, DI, EI )
            IF ( DI .GE. EI ) THEN
               MVNDFN = 0
               RETURN
            ELSE
               MVNDFN = MVNDFN*( EI - DI )
               IF ( I .LE. N ) Y(IK) = PHINVS( DI + W(IK)*( EI - DI ) )
               IK = IK + 1
               INFA = 0
               INFB = 0
            END IF
         END IF
      END DO
      RETURN
*
*     Entry point for initialization.
*
      ENTRY MVNDNT( N, CORREL, LOWER, UPPER, INFIN, INFIS, D, E )
      MVNDNT = 0
*
*     Initialization and computation of covariance Cholesky factor.
*
      CALL COVSRT( N, LOWER,UPPER,CORREL,INFIN,Y, INFIS,A,B,COV,INFI )
      IF ( N - INFIS .EQ. 1 ) THEN
         CALL MVNLMS( A(1), B(1), INFI(1), D, E ) 
      ELSE IF ( N - INFIS .EQ. 2 ) THEN
         IF ( ABS( COV(3) ) .GT. 0 ) THEN
            D = SQRT( 1 + COV(2)**2 )
            IF ( INFI(2) .NE. 0 ) A(2) = A(2)/D
            IF ( INFI(2) .NE. 1 ) B(2) = B(2)/D
            E = BVNMVN( A, B, INFI, COV(2)/D )
            D = 0
         ELSE
            IF ( INFI(1) .NE. 0 ) THEN
               IF ( INFI(2) .NE. 0 ) A(1) = MAX( A(1), A(2) )
            ELSE
               IF ( INFI(2) .NE. 0 ) A(1) = A(2) 
            END IF
            IF ( INFI(1) .NE. 1 ) THEN
               IF ( INFI(2) .NE. 1 ) B(1) = MIN( B(1), B(2) ) 
            ELSE
               IF ( INFI(2) .NE. 1 ) B(1) = B(2)
            END IF
            IF ( INFI(1) .NE. INFI(2) ) INFI(1) = 2
            CALL MVNLMS( A(1), B(1), INFI(1), D, E ) 
         END IF
         INFIS = INFIS + 1 
      END IF
      END
      SUBROUTINE MVNLMS( A, B, INFIN, LOWER, UPPER )
      DOUBLE PRECISION A, B, LOWER, UPPER, MVNPHI
      INTEGER INFIN
      LOWER = 0
      UPPER = 1
      IF ( INFIN .GE. 0 ) THEN
         IF ( INFIN .NE. 0 ) LOWER = MVNPHI(A)
         IF ( INFIN .NE. 1 ) UPPER = MVNPHI(B)
      ENDIF
      UPPER = MAX( UPPER, LOWER )
      END      
      SUBROUTINE COVSRT( N, LOWER, UPPER, CORREL, INFIN, Y, 
     &               INFIS,     A,     B,    COV,  INFI )
*
*     Subroutine to sort integration limits and determine Cholesky factor.
*
      INTEGER N, INFI(*), INFIN(*), INFIS
      DOUBLE PRECISION 
     &     A(*), B(*), COV(*), LOWER(*), UPPER(*), CORREL(*), Y(*)
      INTEGER I, J, K, L, M, II, IJ, IL, JMIN
      DOUBLE PRECISION SUMSQ, AJ, BJ, SUM, SQTWPI, EPS, D, E
      DOUBLE PRECISION CVDIAG, AMIN, BMIN, DMIN, EMIN, YL, YU
      PARAMETER ( SQTWPI = 2.506628274631001D0, EPS = 1D-10 )
      IJ = 0
      II = 0
      INFIS = 0
      DO I = 1, N
         A(I) = 0
         B(I) = 0
         INFI(I) = INFIN(I) 
         IF ( INFI(I) .LT. 0 ) THEN
            INFIS = INFIS + 1
         ELSE 
            IF ( INFI(I) .NE. 0 ) A(I) = LOWER(I)
            IF ( INFI(I) .NE. 1 ) B(I) = UPPER(I)
         ENDIF
         DO J = 1, I-1
            IJ = IJ + 1
            II = II + 1
            COV(IJ) = CORREL(II)
         END DO
         IJ = IJ + 1
         COV(IJ) = 1
      END DO
*
*     First move any doubly infinite limits to innermost positions.
*
      IF ( INFIS .LT. N ) THEN
         DO I = N, N-INFIS+1, -1
            IF ( INFI(I) .GE. 0 ) THEN 
               DO J = 1,I-1
                  IF ( INFI(J) .LT. 0 ) THEN
                     CALL RCSWP( J, I, A, B, INFI, N, COV )
                     GO TO 10
                  ENDIF
               END DO
            ENDIF
 10      END DO
*
*     Sort remaining limits and determine Cholesky factor.
*
         II = 0
         DO I = 1, N-INFIS
*
*        Determine the integration limits for variable with minimum
*        expected probability and interchange that variable with Ith.
*
            DMIN = 0
            EMIN = 1
            JMIN = I
            CVDIAG = 0
            IJ = II
            DO J = I, N-INFIS
               IF ( COV(IJ+J) .GT. EPS ) THEN
                  SUMSQ = SQRT( COV(IJ+J) )
                  SUM = 0
                  DO K = 1, I-1
                     SUM = SUM + COV(IJ+K)*Y(K)
                  END DO
                  AJ = ( A(J) - SUM )/SUMSQ
                  BJ = ( B(J) - SUM )/SUMSQ
                  CALL MVNLMS( AJ, BJ, INFI(J), D, E )
                  IF ( EMIN + D .GE. E + DMIN ) THEN
                     JMIN = J
                     AMIN = AJ
                     BMIN = BJ
                     DMIN = D
                     EMIN = E
                     CVDIAG = SUMSQ
                  ENDIF
               ENDIF
               IJ = IJ + J 
            END DO
            IF ( JMIN .GT. I ) CALL RCSWP( I, JMIN, A,B, INFI, N, COV )
            COV(II+I) = CVDIAG
*
*        Compute Ith column of Cholesky factor.
*        Compute expected value for Ith integration variable and
*         scale Ith covariance matrix row and limits.
*
            IF ( CVDIAG .GT. 0 ) THEN
               IL = II + I
               DO L = I+1, N-INFIS               
                  COV(IL+I) = COV(IL+I)/CVDIAG
                  IJ = II + I
                  DO J = I+1, L
                     COV(IL+J) = COV(IL+J) - COV(IL+I)*COV(IJ+I)
                     IJ = IJ + J
                  END DO
                  IL = IL + L
               END DO
               IF ( EMIN .GT. DMIN + EPS ) THEN
                  YL = 0
                  YU = 0
                  IF ( INFI(I) .NE. 0 ) YL = -EXP( -AMIN**2/2 )/SQTWPI
                  IF ( INFI(I) .NE. 1 ) YU = -EXP( -BMIN**2/2 )/SQTWPI
                  Y(I) = ( YU - YL )/( EMIN - DMIN )
               ELSE
                  IF ( INFI(I) .EQ. 0 ) Y(I) = BMIN
                  IF ( INFI(I) .EQ. 1 ) Y(I) = AMIN
                  IF ( INFI(I) .EQ. 2 ) Y(I) = ( AMIN + BMIN )/2
               END IF
               DO J = 1, I
                  II = II + 1
                  COV(II) = COV(II)/CVDIAG
               END DO
               A(I) = A(I)/CVDIAG
               B(I) = B(I)/CVDIAG
            ELSE
               IL = II + I
               DO L = I+1, N-INFIS                
                  COV(IL+I) = 0
                  IL = IL + L
               END DO
*
*        If the covariance matrix diagonal entry is zero, 
*         permute limits and/or rows, if necessary.
*
*
               DO J = I-1, 1, -1
                  IF ( ABS( COV(II+J) ) .GT. EPS ) THEN
                     A(I) = A(I)/COV(II+J)
                     B(I) = B(I)/COV(II+J)
                     IF ( COV(II+J) .LT. 0 ) THEN
                        CALL DKSWAP( A(I), B(I) ) 
                        IF ( INFI(I) .NE. 2 ) INFI(I) = 1 - INFI(I)
                     END IF
                     DO L = 1, J
                        COV(II+L) = COV(II+L)/COV(II+J)
                     END DO
                     DO L = J+1, I-1 
                        IF( COV((L-1)*L/2+J+1) .GT. 0 ) THEN
                           IJ = II
                           DO K = I-1, L, -1 
                              DO M = 1, K
                                 CALL DKSWAP( COV(IJ-K+M), COV(IJ+M) )
                              END DO
                              CALL DKSWAP( A(K), A(K+1) ) 
                              CALL DKSWAP( B(K), B(K+1) ) 
                              M = INFI(K)
                              INFI(K) = INFI(K+1)
                              INFI(K+1) = M
                              IJ = IJ - K 
                           END DO
                           GO TO 20
                        END IF
                     END DO
                     GO TO 20
                  END IF
                  COV(II+J) = 0
               END DO
 20            II = II + I
               Y(I) = 0
            END IF
         END DO
      ENDIF
      END
*
      SUBROUTINE DKSWAP( X, Y )
      DOUBLE PRECISION X, Y, T
      T = X
      X = Y
      Y = T
      END
*
      SUBROUTINE RCSWP( P, Q, A, B, INFIN, N, C )
*
*     Swaps rows and columns P and Q in situ, with P <= Q.
*
      DOUBLE PRECISION A(*), B(*), C(*)
      INTEGER INFIN(*), P, Q, N, I, J, II, JJ
      CALL DKSWAP( A(P), A(Q) )
      CALL DKSWAP( B(P), B(Q) )
      J = INFIN(P)
      INFIN(P) = INFIN(Q)
      INFIN(Q) = J
      JJ = ( P*( P - 1 ) )/2
      II = ( Q*( Q - 1 ) )/2
      CALL DKSWAP( C(JJ+P), C(II+Q) )
      DO J = 1, P-1
         CALL DKSWAP( C(JJ+J), C(II+J) )
      END DO
      JJ = JJ + P
      DO I = P+1, Q-1
         CALL DKSWAP( C(JJ+P), C(II+I) )
         JJ = JJ + I
      END DO
      II = II + Q
      DO I = Q+1, N
         CALL DKSWAP( C(II+P), C(II+Q) )
         II = II + I
      END DO
      END
*
      SUBROUTINE DKBVRC( NDIM, MINVLS, MAXVLS, FUNCTN, ABSEPS, RELEPS,
     &                   ABSERR, FINEST, INFORM )
*
*  Automatic Multidimensional Integration Subroutine
*               
*         AUTHOR: Alan Genz
*                 Department of Mathematics
*                 Washington State University
*                 Pulman, WA 99164-3113
*                 Email: AlanGenz@wsu.edu
*
*         Last Change: 1/15/03
*
*  KRBVRC computes an approximation to the integral
*
*      1  1     1
*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
*      0  0     0
*
*
*  DKBVRC uses randomized Korobov rules for the first 100 variables. 
*  The primary references are
*   "Randomization of Number Theoretic Methods for Multiple Integration"
*    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
*  and 
*   "Optimal Parameters for Multidimensional Integration", 
*    P. Keast, SIAM J Numer Anal, 10, pp.831-838.
*  If there are more than 100 variables, the remaining variables are
*  integrated using the rules described in the reference
*   "On a Number-Theoretical Integration Method"
*   H. Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11.
*   
***************  Parameters ********************************************
****** Input parameters
*  NDIM    Number of variables, must exceed 1, but not exceed 40
*  MINVLS  Integer minimum number of function evaluations allowed.
*          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the
*          routine assumes a previous call has been made with 
*          the same integrand and continues that calculation.
*  MAXVLS  Integer maximum number of function evaluations allowed.
*  FUNCTN  EXTERNALly declared user defined function to be integrated.
*          It must have parameters (NDIM,Z), where Z is a real array
*          of dimension NDIM.
*                                     
*  ABSEPS  Required absolute accuracy.
*  RELEPS  Required relative accuracy.
****** Output parameters
*  MINVLS  Actual number of function evaluations used.
*  ABSERR  Estimated absolute accuracy of FINEST.
*  FINEST  Estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when 
*                     ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST))
*                  and 
*                     INTVLS <= MAXCLS.
*          INFORM = 1 If MAXVLS was too small to obtain the required 
*          accuracy. In this case a value FINEST is returned with 
*          estimated absolute accuracy ABSERR.
************************************************************************
      EXTERNAL FUNCTN
      INTEGER NDIM, MINVLS, MAXVLS, INFORM, NP, PLIM, NLIM, KLIM, KLIMI,
     &        SAMPLS, I, INTVLS, MINSMP
      PARAMETER ( PLIM = 28, NLIM = 1000, KLIM = 100, MINSMP = 8 )
      INTEGER P(PLIM), C(PLIM,KLIM-1) 
      DOUBLE PRECISION FUNCTN, ABSEPS, RELEPS, FINEST, ABSERR, DIFINT, 
     &                 FINVAL, VARSQR, VAREST, VARPRD, VALUE
      DOUBLE PRECISION X(2*NLIM), VK(NLIM), ONE
      PARAMETER ( ONE = 1 )
      SAVE P, C, SAMPLS, NP, VAREST
      INFORM = 1
      INTVLS = 0
      KLIMI = KLIM
      IF ( MINVLS .GE. 0 ) THEN
         FINEST = 0
         VAREST = 0
         SAMPLS = MINSMP 
         DO I = MIN( NDIM, 10), PLIM
            NP = I
            IF ( MINVLS .LT. 2*SAMPLS*P(I) ) GO TO 10
         END DO
         SAMPLS = MAX( MINSMP, MINVLS/( 2*P(NP) ) )
      ENDIF
 10   VK(1) = ONE/P(NP)
      DO I = 2, NDIM
         IF ( I .LE. KLIM ) THEN
            VK(I) = MOD( C(NP, MIN(NDIM-1,KLIM-1))*VK(I-1), ONE )
         ELSE
            VK(I) = INT( P(NP)*2**(DBLE(I-KLIM)/(NDIM-KLIM+1)) ) 
            VK(I) = MOD( VK(I)/P(NP), ONE ) 
         END IF
      END DO
      FINVAL = 0
      VARSQR = 0
      DO I = 1, SAMPLS
         CALL DKSMRC( NDIM, KLIMI, VALUE, P(NP), VK, FUNCTN, X )
         DIFINT = ( VALUE - FINVAL )/I
         FINVAL = FINVAL + DIFINT
         VARSQR = ( I - 2 )*VARSQR/I + DIFINT**2
      END DO
      INTVLS = INTVLS + 2*SAMPLS*P(NP)
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/( 1 + VARPRD )
      IF ( VARSQR .GT. 0 ) VAREST = ( 1 + VARPRD )/VARSQR
      ABSERR = 7*SQRT( VARSQR/( 1 + VARPRD ) )/2
      IF ( ABSERR .GT. MAX( ABSEPS, ABS(FINEST)*RELEPS ) ) THEN
         IF ( NP .LT. PLIM ) THEN
            NP = NP + 1
         ELSE
            SAMPLS = MIN( 3*SAMPLS/2, ( MAXVLS - INTVLS )/( 2*P(NP) ) ) 
            SAMPLS = MAX( MINSMP, SAMPLS )
         ENDIF
         IF ( INTVLS + 2*SAMPLS*P(NP) .LE. MAXVLS ) GO TO 10
      ELSE
         INFORM = 0
      ENDIF
      MINVLS = INTVLS
*
*    Optimal Parameters for Lattice Rules
*
      DATA P( 1),(C( 1,I),I = 1,99)/     31, 12, 2*9, 13, 8*12, 3*3, 12,
     & 2*7, 9*12, 3*3, 12, 2*7, 9*12, 3*3, 12, 2*7, 9*12, 3*3, 12, 2*7,
     & 8*12, 7, 3*3, 3*7, 21*3/
      DATA P( 2),(C( 2,I),I = 1,99)/    47, 13, 11, 17, 10, 6*15,
     & 22, 2*15, 3*6, 2*15, 9, 13, 3*2, 13, 2*11, 10, 9*15, 3*6, 2*15,
     & 9, 13, 3*2, 13, 2*11, 10, 9*15, 3*6, 2*15, 9, 13, 3*2, 13, 2*11,
     & 2*10, 8*15, 6, 2, 3, 2, 3, 12*2/
      DATA P( 3),(C( 3,I),I = 1,99)/    73, 27, 28, 10, 2*11, 20,
     & 2*11, 28, 2*13, 28, 3*13, 16*14, 2*31, 3*5, 31, 13, 6*11, 7*13,
     & 16*14, 2*31, 3*5, 11, 13, 7*11, 2*13, 11, 13, 4*5, 14, 13, 8*5/
      DATA P( 4),(C( 4,I),I = 1,99)/   113, 35, 2*27, 36, 22, 2*29,
     & 20, 45, 3*5, 16*21, 29, 10*17, 12*23, 21, 27, 3*3, 24, 2*27,
     & 17, 3*29, 17, 4*5, 16*21, 3*17, 6, 2*17, 6, 3, 2*6, 5*3/
      DATA P( 5),(C( 5,I),I = 1,99)/   173, 64, 66, 2*28, 2*44, 55,
     & 67, 6*10, 2*38, 5*10, 12*49, 2*38, 31, 2*4, 31, 64, 3*4, 64,
     & 6*45, 19*66, 11, 9*66, 45, 11, 7, 3, 3*2, 27, 5, 2*3, 2*5, 7*2/
      DATA P( 6),(C( 6,I),I = 1,99)/   263, 111, 42, 54, 118, 20,
     & 2*31, 72, 17, 94, 2*14, 11, 3*14, 94, 4*10, 7*14, 3*11, 7*8,
     & 5*18, 113, 2*62, 2*45, 17*113, 2*63, 53, 63, 15*67, 5*51, 12,
     & 51, 12, 51, 5, 2*3, 2*2, 5/
      DATA P( 7),(C( 7,I),I = 1,99)/   397, 163, 154, 83, 43, 82,
     & 92, 150, 59, 2*76, 47, 2*11, 100, 131, 6*116, 9*138, 21*101,
     & 6*116, 5*100, 5*138, 19*101, 8*38, 5*3/
      DATA P( 8),(C( 8,I),I = 1,99)/   593, 246, 189, 242, 102,
     & 2*250, 102, 250, 280, 118, 196, 118, 191, 215, 2*121,
     & 12*49, 34*171, 8*161, 17*14, 6*10, 103, 4*10, 5/
      DATA P( 9),(C( 9,I),I = 1,99)/   907, 347, 402, 322, 418,
     & 215, 220, 3*339, 337, 218, 4*315, 4*167, 361, 201, 11*124,
     & 2*231, 14*90, 4*48, 23*90, 10*243, 9*283, 16, 283, 16, 2*283/
      DATA P(10),(C(10,I),I = 1,99)/  1361, 505, 220, 601, 644,
     & 612, 160, 3*206, 422, 134, 518, 2*134, 518, 652, 382,
     & 206, 158, 441, 179, 441, 56, 2*559, 14*56, 2*101, 56,
     & 8*101, 7*193, 21*101, 17*122, 4*101/
      DATA P(11),(C(11,I),I = 1,99)/  2053, 794, 325, 960, 528,
     & 2*247, 338, 366, 847, 2*753, 236, 2*334, 461, 711, 652,
     & 3*381, 652, 7*381, 226, 7*326, 126, 10*326, 2*195, 19*55,
     & 7*195, 11*132, 13*387/
      DATA P(12),(C(12,I),I = 1,99)/  3079, 1189, 888, 259, 1082, 725,      
     & 811, 636, 965, 2*497, 2*1490, 392, 1291, 2*508, 2*1291, 508,
     & 1291, 2*508, 4*867, 934, 7*867, 9*1284, 4*563, 3*1010, 208,
     & 838, 3*563, 2*759, 564, 2*759, 4*801, 5*759, 8*563, 22*226/
      DATA P(13),(C(13,I),I = 1,99)/  4621, 1763, 1018, 1500, 432,
     & 1332, 2203, 126, 2240, 1719, 1284, 878, 1983, 4*266,
     & 2*747, 2*127, 2074, 127, 2074, 1400, 10*1383, 1400, 7*1383,
     & 507, 4*1073, 5*1990, 9*507, 17*1073, 6*22, 1073, 6*452, 318,
     & 4*301, 2*86, 15/
      DATA P(14),(C(14,I),I = 1,99)/  6947, 2872, 3233, 1534, 2941,
     & 2910, 393, 1796, 919, 446, 2*919, 1117, 7*103, 2311, 3117, 1101,
     & 2*3117, 5*1101, 8*2503, 7*429, 3*1702, 5*184, 34*105, 13*784/
      DATA P(15),(C(15,I),I = 1,99)/ 10427, 4309, 3758, 4034, 1963,
     & 730, 642, 1502, 2246, 3834, 1511, 2*1102, 2*1522, 2*3427,
     & 3928, 2*915, 4*3818, 3*4782, 3818, 4782, 2*3818, 7*1327, 9*1387,
     & 13*2339, 18*3148, 3*1776, 3*3354, 925, 2*3354, 5*925, 8*2133/
      DATA P(16),(C(16,I),I = 1,99)/ 15641, 6610, 6977, 1686, 3819,
     & 2314, 5647, 3953, 3614, 5115, 2*423, 5408, 7426, 2*423,
     & 487, 6227, 2660, 6227, 1221, 3811, 197, 4367, 351,
     & 1281, 1221, 3*351, 7245, 1984, 6*2999, 3995, 4*2063, 1644,
     & 2063, 2077, 3*2512, 4*2077, 19*754, 2*1097, 4*754, 248, 754,
     & 4*1097, 4*222, 754,11*1982/
      DATA P(17),(C(17,I),I = 1,99)/ 23473, 9861, 3647, 4073, 2535,
     & 3430, 9865, 2830, 9328, 4320, 5913, 10365, 8272, 3706, 6186,
     & 3*7806, 8610, 2563, 2*11558, 9421, 1181, 9421, 3*1181, 9421,
     & 2*1181, 2*10574, 5*3534, 3*2898, 3450, 7*2141, 15*7055, 2831,
     & 24*8204, 3*4688, 8*2831/
      DATA P(18),(C(18,I),I = 1,99)/ 35221, 10327, 7582, 7124, 8214,
     & 9600, 10271, 10193, 10800, 9086, 2365, 4409, 13812,
     & 5661, 2*9344, 10362, 2*9344, 8585, 11114, 3*13080, 6949,
     & 3*3436, 13213, 2*6130, 2*8159, 11595, 8159, 3436, 18*7096,
     & 4377, 7096, 5*4377, 2*5410, 32*4377, 2*440, 3*1199/
      DATA P(19),(C(19,I),I = 1,99)/ 52837, 19540, 19926, 11582,
     & 11113, 24585, 8726, 17218, 419, 3*4918, 15701, 17710,
     & 2*4037, 15808, 11401, 19398, 2*25950, 4454, 24987, 11719,
     & 8697, 5*1452, 2*8697, 6436, 21475, 6436, 22913, 6434, 18497,
     & 4*11089, 2*3036, 4*14208, 8*12906, 4*7614, 6*5021, 24*10145,
     & 6*4544, 4*8394/    
      DATA P(20),(C(20,I),I = 1,99)/ 79259, 34566, 9579, 12654,
     & 26856, 37873, 38806, 29501, 17271, 3663, 10763, 18955,
     & 1298, 26560, 2*17132, 2*4753, 8713, 18624, 13082, 6791,
     & 1122, 19363, 34695, 4*18770, 15628, 4*18770, 33766, 6*20837,
     & 5*6545, 14*12138, 5*30483, 19*12138, 9305, 13*11107, 2*9305/
      DATA P(21),(C(21,I),I = 1,99)/118891, 31929, 49367, 10982, 3527,
     & 27066, 13226, 56010, 18911, 40574, 2*20767, 9686, 2*47603, 
     & 2*11736, 41601, 12888, 32948, 30801, 44243, 2*53351, 16016, 
     & 2*35086, 32581, 2*2464, 49554, 2*2464, 2*49554, 2464, 81, 27260, 
     & 10681, 7*2185, 5*18086, 2*17631, 3*18086, 37335, 3*37774, 
     & 13*26401, 12982, 6*40398, 3*3518, 9*37799, 4*4721, 4*7067/
      DATA P(22),(C(22,I),I = 1,99)/178349, 40701, 69087, 77576, 64590, 
     & 39397, 33179, 10858, 38935, 43129, 2*35468, 5279, 2*61518, 27945,
     & 2*70975, 2*86478, 2*20514, 2*73178, 2*43098, 4701,
     & 2*59979, 58556, 69916, 2*15170, 2*4832, 43064, 71685, 4832,
     & 3*15170, 3*27679, 2*60826, 2*6187, 5*4264, 45567, 4*32269,
     & 9*62060, 13*1803, 12*51108, 2*55315, 5*54140, 13134/
      DATA P(23),(C(23,I),I = 1,99)/267523, 103650, 125480, 59978,
     & 46875, 77172, 83021, 126904, 14541, 56299, 43636, 11655,
     & 52680, 88549, 29804, 101894, 113675, 48040, 113675,
     & 34987, 48308, 97926, 5475, 49449, 6850, 2*62545, 9440,
     & 33242, 9440, 33242, 9440, 33242, 9440, 62850, 3*9440,
     & 3*90308, 9*47904, 7*41143, 5*36114, 24997, 14*65162, 7*47650,
     & 7*40586, 4*38725, 5*88329/
      DATA P(24),(C(24,I),I = 1,99)/401287, 165843, 90647, 59925,
     & 189541, 67647, 74795, 68365, 167485, 143918, 74912,
     & 167289, 75517, 8148, 172106, 126159,3*35867, 121694,
     & 52171, 95354, 2*113969, 76304, 2*123709, 144615, 123709,
     & 2*64958, 32377, 2*193002, 25023, 40017, 141605, 2*189165,
     & 141605, 2*189165, 3*141605, 189165, 20*127047, 10*127785,
     & 6*80822, 16*131661, 7114, 131661/
      DATA P(25),(C(25,I),I = 1,99)/601942, 130365, 236711, 110235,
     & 125699, 56483, 93735, 234469, 60549, 1291, 93937,
     & 245291, 196061, 258647, 162489, 176631, 204895, 73353,
     & 172319, 28881, 136787,2*122081, 275993, 64673, 3*211587,
     & 2*282859, 211587, 242821, 3*256865, 122203, 291915, 122203,
     & 2*291915, 122203, 2*25639, 291803, 245397, 284047,
     & 7*245397, 94241, 2*66575, 19*217673, 10*210249, 15*94453/
      DATA P(26),(C(26,I),I = 1,99)/902933, 333459, 375354, 102417,            
     & 383544, 292630, 41147, 374614, 48032, 435453, 281493, 358168, 
     & 114121, 346892, 238990, 317313, 164158, 35497, 2*70530, 434839,  
     & 3*24754, 393656, 2*118711, 148227, 271087, 355831, 91034, 
     & 2*417029, 2*91034, 417029, 91034, 2*299843, 2*413548, 308300,  
     & 3*413548, 3*308300, 413548, 5*308300, 4*15311, 2*176255, 6*23613, 
     & 172210, 4* 204328, 5*121626, 5*200187, 2*121551, 12*248492, 
     & 5*13942/
      DATA P(27), (C(27,I), I = 1,99)/ 1354471, 500884, 566009, 399251,
     & 652979, 355008, 430235, 328722, 670680, 2*405585, 424646, 
     & 2*670180, 641587, 215580, 59048, 633320, 81010, 20789, 2*389250,  
     & 2*638764, 2*389250, 398094, 80846, 2*147776, 296177, 2*398094,  
     & 2*147776, 396313, 3*578233, 19482, 620706, 187095, 620706, 
     & 187095, 126467, 12*241663, 321632, 2*23210, 3*394484, 3*78101, 
     & 19*542095, 3*277743, 12*457259/
      DATA P(28), (C(28,I), I = 1, 99)/ 2031713, 858339, 918142, 501970, 
     & 234813, 460565, 31996, 753018, 256150, 199809, 993599, 245149,      
     & 794183, 121349, 150619, 376952, 2*809123, 804319, 67352, 969594, 
     & 434796, 969594, 804319, 391368, 761041, 754049, 466264, 2*754049,
     & 466264, 2*754049, 282852, 429907, 390017, 276645, 994856, 250142, 
     & 144595, 907454, 689648, 4*687580, 978368, 687580, 552742, 105195, 
     & 942843, 768249, 4*307142, 7*880619, 11*117185, 11*60731,  
     & 4*178309, 8*74373, 3*214965/
*
      END
*
      SUBROUTINE DKSMRC( NDIM, KLIM, SUMKRO, PRIME, VK, FUNCTN, X )
      EXTERNAL FUNCTN
      INTEGER NDIM, NK, KLIM, PRIME, K, J, JP
      DOUBLE PRECISION SUMKRO, VK(*), FUNCTN, X(*), ONE, XT, MVNUNI
      PARAMETER ( ONE = 1 )
      SUMKRO = 0
      NK = MIN( NDIM, KLIM )
      DO J = 1, NK - 1
         JP = J + MVNUNI()*( NK + 1 - J ) 
         XT = VK(J)
         VK(J) = VK(JP)
         VK(JP) = XT
      END DO
      DO J = 1, NDIM
         X(NDIM+J) = MVNUNI()
      END DO
      DO K = 1, PRIME
         DO J = 1, NDIM
            X(J) = ABS( 2*MOD( K*VK(J) + X(NDIM+J), ONE ) - 1 )
         END DO
         SUMKRO = SUMKRO + ( FUNCTN(NDIM,X) - SUMKRO )/( 2*K - 1 )
         DO J = 1, NDIM
            X(J) = 1 - X(J)
         END DO
         SUMKRO = SUMKRO + ( FUNCTN(NDIM,X) - SUMKRO )/( 2*K )
      END DO
      END
*
      DOUBLE PRECISION FUNCTION MVNPHI( Z )
*     
*     Normal distribution probabilities accurate to 1.e-15.
*     Z = no. of standard deviations from the mean.
*     
*     Based upon algorithm 5666 for the error function, from:
*     Hart, J.F. et al, 'Computer Approximations', Wiley 1968
*     
*     Programmer: Alan Miller
*     
*     Latest revision - 30 March 1986
*     
      DOUBLE PRECISION P0, P1, P2, P3, P4, P5, P6, 
     *     Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     *     Z, P, EXPNTL, CUTOFF, ROOTPI, ZABS
      PARAMETER(
     *     P0 = 220.20 68679 12376 1D0,
     *     P1 = 221.21 35961 69931 1D0, 
     *     P2 = 112.07 92914 97870 9D0,
     *     P3 = 33.912 86607 83830 0D0,
     *     P4 = 6.3739 62203 53165 0D0,
     *     P5 = .70038 30644 43688 1D0, 
     *     P6 = .035262 49659 98910 9D0 )
      PARAMETER(
     *     Q0 = 440.41 37358 24752 2D0,
     *     Q1 = 793.82 65125 19948 4D0, 
     *     Q2 = 637.33 36333 78831 1D0,
     *     Q3 = 296.56 42487 79673 7D0, 
     *     Q4 = 86.780 73220 29460 8D0,
     *     Q5 = 16.064 17757 92069 5D0, 
     *     Q6 = 1.7556 67163 18264 2D0,
     *     Q7 = .088388 34764 83184 4D0 )
      PARAMETER( ROOTPI = 2.5066 28274 63100 1D0 )
      PARAMETER( CUTOFF = 7.0710 67811 86547 5D0 )
*     
      ZABS = ABS(Z)
*     
*     |Z| > 37
*     
      IF ( ZABS .GT. 37 ) THEN
         P = 0
      ELSE
*     
*     |Z| <= 37
*     
         EXPNTL = EXP( -ZABS**2/2 )
*     
*     |Z| < CUTOFF = 10/SQRT(2)
*     
         IF ( ZABS .LT. CUTOFF ) THEN
            P = EXPNTL*( (((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS
     *           + P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS
     *           + Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS
     *           + Q0 )
*     
*     |Z| >= CUTOFF.
*     
         ELSE
            P = EXPNTL/( ZABS + 1/( ZABS + 2/( ZABS + 3/( ZABS 
     *                        + 4/( ZABS + 0.65D0 ) ) ) ) )/ROOTPI
         END IF
      END IF
      IF ( Z .GT. 0 ) P = 1 - P
      MVNPHI = P
      END
      DOUBLE PRECISION FUNCTION PHINVS(P)
*
*	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
*
*	Produces the normal deviate Z corresponding to a given lower
*	tail area of P.
*
*	The hash sums below are the sums of the mantissas of the
*	coefficients.   They are included for use in checking
*	transcription.
*
      DOUBLE PRECISION SPLIT1, SPLIT2, CONST1, CONST2, 
     *     A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7, 
     *     C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7, 
     *     E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, 
     *     P, Q, R
      PARAMETER ( SPLIT1 = 0.425, SPLIT2 = 5,
     *            CONST1 = 0.180625D0, CONST2 = 1.6D0 )
*     
*     Coefficients for P close to 0.5
*     
      PARAMETER (
     *     A0 = 3.38713 28727 96366 6080D0,
     *     A1 = 1.33141 66789 17843 7745D+2,
     *     A2 = 1.97159 09503 06551 4427D+3,
     *     A3 = 1.37316 93765 50946 1125D+4,
     *     A4 = 4.59219 53931 54987 1457D+4,
     *     A5 = 6.72657 70927 00870 0853D+4,
     *     A6 = 3.34305 75583 58812 8105D+4,
     *     A7 = 2.50908 09287 30122 6727D+3,
     *     B1 = 4.23133 30701 60091 1252D+1,
     *     B2 = 6.87187 00749 20579 0830D+2,
     *     B3 = 5.39419 60214 24751 1077D+3,
     *     B4 = 2.12137 94301 58659 5867D+4,
     *     B5 = 3.93078 95800 09271 0610D+4,
     *     B6 = 2.87290 85735 72194 2674D+4,
     *     B7 = 5.22649 52788 52854 5610D+3 )
*     HASH SUM AB    55.88319 28806 14901 4439
*     
*     Coefficients for P not close to 0, 0.5 or 1.
*     
      PARAMETER (
     *     C0 = 1.42343 71107 49683 57734D0,
     *     C1 = 4.63033 78461 56545 29590D0,
     *     C2 = 5.76949 72214 60691 40550D0,
     *     C3 = 3.64784 83247 63204 60504D0,
     *     C4 = 1.27045 82524 52368 38258D0,
     *     C5 = 2.41780 72517 74506 11770D-1,
     *     C6 = 2.27238 44989 26918 45833D-2,
     *     C7 = 7.74545 01427 83414 07640D-4,
     *     D1 = 2.05319 16266 37758 82187D0,
     *     D2 = 1.67638 48301 83803 84940D0,
     *     D3 = 6.89767 33498 51000 04550D-1,
     *     D4 = 1.48103 97642 74800 74590D-1,
     *     D5 = 1.51986 66563 61645 71966D-2,
     *     D6 = 5.47593 80849 95344 94600D-4,
     *     D7 = 1.05075 00716 44416 84324D-9 )
*     HASH SUM CD    49.33206 50330 16102 89036
*
*	Coefficients for P near 0 or 1.
*
      PARAMETER (
     *     E0 = 6.65790 46435 01103 77720D0,
     *     E1 = 5.46378 49111 64114 36990D0,
     *     E2 = 1.78482 65399 17291 33580D0,
     *     E3 = 2.96560 57182 85048 91230D-1,
     *     E4 = 2.65321 89526 57612 30930D-2,
     *     E5 = 1.24266 09473 88078 43860D-3,
     *     E6 = 2.71155 55687 43487 57815D-5,
     *     E7 = 2.01033 43992 92288 13265D-7,
     *     F1 = 5.99832 20655 58879 37690D-1,
     *     F2 = 1.36929 88092 27358 05310D-1,
     *     F3 = 1.48753 61290 85061 48525D-2,
     *     F4 = 7.86869 13114 56132 59100D-4,
     *     F5 = 1.84631 83175 10054 68180D-5,
     *     F6 = 1.42151 17583 16445 88870D-7,
     *     F7 = 2.04426 31033 89939 78564D-15 )
*     HASH SUM EF    47.52583 31754 92896 71629
*     
      Q = ( 2*P - 1 )/2
      IF ( ABS(Q) .LE. SPLIT1 ) THEN
         R = CONST1 - Q*Q
         PHINVS = Q*( ( ( ((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     *                  *R + A2 )*R + A1 )*R + A0 )
     *            /( ( ( ((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     *                  *R + B2 )*R + B1 )*R + 1 )
      ELSE
         R = MIN( P, 1 - P )
         IF ( R .GT. 0 ) THEN
            R = SQRT( -LOG(R) )
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               PHINVS = ( ( ( ((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     *                      *R + C2 )*R + C1 )*R + C0 ) 
     *                /( ( ( ((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     *                      *R + D2 )*R + D1 )*R + 1 )
            ELSE
               R = R - SPLIT2
               PHINVS = ( ( ( ((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     *                      *R + E2 )*R + E1 )*R + E0 )
     *                /( ( ( ((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     *                      *R + F2 )*R + F1 )*R + 1 )
            END IF
         ELSE
            PHINVS = 9
         END IF
         IF ( Q .LT. 0 ) PHINVS = - PHINVS
      END IF
      END
      DOUBLE PRECISION FUNCTION BVNMVN( LOWER, UPPER, INFIN, CORREL )
*
*     A function for computing bivariate normal probabilities.
*
*  Parameters
*
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, correlation coefficient.
*
      DOUBLE PRECISION LOWER(*), UPPER(*), CORREL, BVU
      INTEGER INFIN(*)
      IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 2 ) THEN
         BVNMVN =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( UPPER(1), LOWER(2), CORREL )
     +           - BVU ( LOWER(1), UPPER(2), CORREL )
     +           + BVU ( UPPER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 1 ) THEN
         BVNMVN =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( UPPER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 2 ) THEN
         BVNMVN =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( LOWER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 0 ) THEN
         BVNMVN =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - BVU ( -LOWER(1), -UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 2 ) THEN
         BVNMVN =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - BVU ( -UPPER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 0 ) THEN
         BVNMVN =  BVU ( LOWER(1), -UPPER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 1 ) THEN
         BVNMVN =  BVU ( -UPPER(1), LOWER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 1 ) THEN
         BVNMVN =  BVU ( LOWER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 0 ) THEN
         BVNMVN =  BVU ( -UPPER(1), -UPPER(2), CORREL )
      END IF
      END 
      DOUBLE PRECISION FUNCTION BVU( SH, SK, R )
*
*     A function for computing bivariate normal probabilities.
*
*       Yihong Ge
*       Department of Computer Science and Electrical Engineering
*       Washington State University
*       Pullman, WA 99164-2752
*     and
*       Alan Genz
*       Department of Mathematics
*       Washington State University
*       Pullman, WA 99164-3113
*       Email : alangenz@wsu.edu
*
* BVN - calculate the probability that X is larger than SH and Y is
*       larger than SK.
*
* Parameters
*
*   SH  REAL, integration limit
*   SK  REAL, integration limit
*   R   REAL, correlation coefficient
*   LG  INTEGER, number of Gauss Rule Points and Weights
*
      DOUBLE PRECISION BVN, SH, SK, R, ZERO, TWOPI 
      INTEGER I, LG, NG
      PARAMETER ( ZERO = 0, TWOPI = 6.283185307179586D0 ) 
      DOUBLE PRECISION X(10,3), W(10,3), AS, A, B, C, D, RS, XS
      DOUBLE PRECISION MVNPHI, SN, ASR, H, K, BS, HS, HK
      SAVE X, W
*     Gauss Legendre Points and Weights, N =  6
      DATA ( W(I,1), X(I,1), I = 1,3) /
     *  0.1713244923791705D+00,-0.9324695142031522D+00,
     *  0.3607615730481384D+00,-0.6612093864662647D+00,
     *  0.4679139345726904D+00,-0.2386191860831970D+00/
*     Gauss Legendre Points and Weights, N = 12
      DATA ( W(I,2), X(I,2), I = 1,6) /
     *  0.4717533638651177D-01,-0.9815606342467191D+00,
     *  0.1069393259953183D+00,-0.9041172563704750D+00,
     *  0.1600783285433464D+00,-0.7699026741943050D+00,
     *  0.2031674267230659D+00,-0.5873179542866171D+00,
     *  0.2334925365383547D+00,-0.3678314989981802D+00,
     *  0.2491470458134029D+00,-0.1252334085114692D+00/
*     Gauss Legendre Points and Weights, N = 20
      DATA ( W(I,3), X(I,3), I = 1,10) /
     *  0.1761400713915212D-01,-0.9931285991850949D+00,
     *  0.4060142980038694D-01,-0.9639719272779138D+00,
     *  0.6267204833410906D-01,-0.9122344282513259D+00,
     *  0.8327674157670475D-01,-0.8391169718222188D+00,
     *  0.1019301198172404D+00,-0.7463319064601508D+00,
     *  0.1181945319615184D+00,-0.6360536807265150D+00,
     *  0.1316886384491766D+00,-0.5108670019508271D+00,
     *  0.1420961093183821D+00,-0.3737060887154196D+00,
     *  0.1491729864726037D+00,-0.2277858511416451D+00,
     *  0.1527533871307259D+00,-0.7652652113349733D-01/
      IF ( ABS(R) .LT. 0.3 ) THEN
         NG = 1
         LG = 3
      ELSE IF ( ABS(R) .LT. 0.75 ) THEN
         NG = 2
         LG = 6
      ELSE 
         NG = 3
         LG = 10
      ENDIF
      H = SH
      K = SK 
      HK = H*K
      BVN = 0
      IF ( ABS(R) .LT. 0.925 ) THEN
         HS = ( H*H + K*K )/2
         ASR = ASIN(R)
         DO I = 1, LG
            SN = SIN(ASR*( X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
            SN = SIN(ASR*(-X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
         END DO
         BVN = BVN*ASR/(2*TWOPI) + MVNPHI(-H)*MVNPHI(-K) 
      ELSE
         IF ( R .LT. 0 ) THEN
            K = -K
            HK = -HK
         ENDIF
         IF ( ABS(R) .LT. 1 ) THEN
            AS = ( 1 - R )*( 1 + R )
            A = SQRT(AS)
            BS = ( H - K )**2
            C = ( 4 - HK )/8 
            D = ( 12 - HK )/16
            BVN = A*EXP( -(BS/AS + HK)/2 )
     +             *( 1 - C*(BS - AS)*(1 - D*BS/5)/3 + C*D*AS*AS/5 )
            IF ( HK .GT. -160 ) THEN
               B = SQRT(BS)
               BVN = BVN - EXP(-HK/2)*SQRT(TWOPI)*MVNPHI(-B/A)*B
     +                    *( 1 - C*BS*( 1 - D*BS/5 )/3 ) 
            ENDIF
            A = A/2
            DO I = 1, LG
               XS = ( A*(X(I,NG)+1) )**2
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I,NG)*
     +              ( EXP( -BS/(2*XS) - HK/(1+RS) )/RS 
     +              - EXP( -(BS/XS+HK)/2 )*( 1 + C*XS*( 1 + D*XS ) ) )
               XS = AS*(-X(I,NG)+1)**2/4
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I,NG)*EXP( -(BS/XS + HK)/2 )
     +                    *( EXP( -HK*(1-RS)/(2*(1+RS)) )/RS 
     +                       - ( 1 + C*XS*( 1 + D*XS ) ) )
            END DO
            BVN = -BVN/TWOPI
         ENDIF
         IF ( R .GT. 0 ) BVN =  BVN + MVNPHI( -MAX( H, K ) )
         IF ( R .LT. 0 ) BVN = -BVN + MAX( ZERO, MVNPHI(-H)-MVNPHI(-K) )
      ENDIF
      BVU = BVN
      END
      DOUBLE PRECISION FUNCTION MVNUNI()
*
*     Uniform (0,1) random number generator
*
*     Reference:
*     L'Ecuyer, Pierre (1996), 
*     "Combined Multiple Recursive Random Number Generators"
*     Operations Research 44, pp. 816-822.
*
*
      INTEGER A12, A13, A21, A23, P12, P13, P21, P23
      INTEGER Q12, Q13, Q21, Q23, R12, R13, R21, R23
      INTEGER X10, X11, X12, X20, X21, X22, Z, M1, M2, H 
      DOUBLE PRECISION INVMP1
      PARAMETER ( M1 = 2147483647, M2 = 2145483479 )
      PARAMETER ( A12 =   63308, Q12 = 33921, R12 = 12979 )
      PARAMETER ( A13 = -183326, Q13 = 11714, R13 =  2883 )
      PARAMETER ( A21 =   86098, Q21 = 24919, R21 =  7417 )
      PARAMETER ( A23 = -539608, Q23 =  3976, R23 =  2071 )
      PARAMETER ( INVMP1 = 4.656612873077392578125D-10 ) 
*                 INVMP1 = 1/(M1+1)
      SAVE X10, X11, X12, X20, X21, X22
      DATA       X10,      X11,      X12,      X20,      X21,      X22
     &    / 15485857, 17329489, 36312197, 55911127, 75906931, 96210113 /
*
*     Component 1
*
      H = X10/Q13
      P13 = -A13*( X10 - H*Q13 ) - H*R13
      H = X11/Q12
      P12 =  A12*( X11 - H*Q12 ) - H*R12
      IF ( P13 .LT. 0 ) P13 = P13 + M1
      IF ( P12 .LT. 0 ) P12 = P12 + M1
      X10 = X11 
      X11 = X12
      X12 = P12 - P13
      IF ( X12 .LT. 0 ) X12 = X12 + M1
*
*     Component 2
*
      H = X20/Q23
      P23 = -A23*( X20 - H*Q23 ) - H*R23
      H = X22/Q21
      P21 =  A21*( X22 - H*Q21 ) - H*R21
      IF ( P23 .LT. 0 ) P23 = P23 + M2
      IF ( P21 .LT. 0 ) P21 = P21 + M2
      X20 = X21 
      X21 = X22
      X22 = P21 - P23
      IF ( X22 .LT. 0 ) X22 = X22 + M2
*
*     Combination
*
      Z = X12 - X22
      IF ( Z .LE. 0 ) Z = Z + M1
      MVNUNI = Z*INVMP1
      END
