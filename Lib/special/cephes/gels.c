/*
C
C     ..................................................................
C
C        SUBROUTINE GELS
C
C        PURPOSE
C           TO SOLVE A SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS WITH
C           SYMMETRIC COEFFICIENT MATRIX UPPER TRIANGULAR PART OF WHICH
C           IS ASSUMED TO BE STORED COLUMNWISE.
C
C        USAGE
C           CALL GELS(R,A,M,N,EPS,IER,AUX)
C
C        DESCRIPTION OF PARAMETERS
C           R      - M BY N RIGHT HAND SIDE MATRIX.  (DESTROYED)
C                    ON RETURN R CONTAINS THE SOLUTION OF THE EQUATIONS.
C           A      - UPPER TRIANGULAR PART OF THE SYMMETRIC
C                    M BY M COEFFICIENT MATRIX.  (DESTROYED)
C           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.
C           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS RELATIVE
C                    TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS
C                    IER=0  - NO ERROR,
C                    IER=-1 - NO RESULT BECAUSE OF M LESS THAN 1 OR
C                             PIVOT ELEMENT AT ANY ELIMINATION STEP
C                             EQUAL TO 0,
C                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-
C                             CANCE INDICATED AT ELIMINATION STEP K+1,
C                             WHERE PIVOT ELEMENT WAS LESS THAN OR
C                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES
C                             ABSOLUTELY GREATEST MAIN DIAGONAL
C                             ELEMENT OF MATRIX A.
C           AUX    - AN AUXILIARY STORAGE ARRAY WITH DIMENSION M-1.
C
C        REMARKS
C           UPPER TRIANGULAR PART OF MATRIX A IS ASSUMED TO BE STORED
C           COLUMNWISE IN M*(M+1)/2 SUCCESSIVE STORAGE LOCATIONS, RIGHT
C           HAND SIDE MATRIX R COLUMNWISE IN N*M SUCCESSIVE STORAGE
C           LOCATIONS. ON RETURN SOLUTION MATRIX R IS STORED COLUMNWISE
C           TOO.
C           THE PROCEDURE GIVES RESULTS IF THE NUMBER OF EQUATIONS M IS
C           GREATER THAN 0 AND PIVOT ELEMENTS AT ALL ELIMINATION STEPS
C           ARE DIFFERENT FROM 0. HOWEVER WARNING IER=K - IF GIVEN -
C           INDICATES POSSIBLE LOSS OF SIGNIFICANCE. IN CASE OF A WELL
C           SCALED MATRIX A AND APPROPRIATE TOLERANCE EPS, IER=K MAY BE
C           INTERPRETED THAT MATRIX A HAS THE RANK K. NO WARNING IS
C           GIVEN IN CASE M=1.
C           ERROR PARAMETER IER=-1 DOES NOT NECESSARILY MEAN THAT
C           MATRIX A IS SINGULAR, AS ONLY MAIN DIAGONAL ELEMENTS
C           ARE USED AS PIVOT ELEMENTS. POSSIBLY SUBROUTINE GELG (WHICH
C           WORKS WITH TOTAL PIVOTING) WOULD BE ABLE TO FIND A SOLUTION.
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           SOLUTION IS DONE BY MEANS OF GAUSS-ELIMINATION WITH
C           PIVOTING IN MAIN DIAGONAL, IN ORDER TO PRESERVE
C           SYMMETRY IN REMAINING COEFFICIENT MATRICES.
C
C     ..................................................................
C
*/
#define ANSIPROT
#ifndef ANSIPROT
double fabs();
#else
extern double fabs(double);
int gels( double [], double [], int, double, double [] );
#endif

int
gels( A, R, M, EPS, AUX )
double A[],R[];
int M;
double EPS;
double AUX[];
{
int I = 0, J = 0, K, L, IER;
int II, LL, LLD, LR, LT, LST, LLST, LEND;
double tb, piv, tol, pivi;

if( M <= 0 )
	{
fatal:
	IER = -1;
	goto done;
	}
/* SEARCH FOR GREATEST MAIN DIAGONAL ELEMENT */

/*  Diagonal elements are at A(i,i) = 1, 3, 6, 10, ...
 *  A(i,j) = A( i(i-1)/2 + j )
 */
IER = 0;
piv = 0.0;
L = 0;
for( K=1; K<=M; K++ )
	{
	L += K;
	tb = fabs( A[L-1] );
	if( tb > piv )
		{
		piv = tb;
		I = L;
		J = K;
		}
	}
tol = EPS * piv;

/*
C     MAIN DIAGONAL ELEMENT A(I)=A(J,J) IS FIRST PIVOT ELEMENT.
C     PIV CONTAINS THE ABSOLUTE VALUE OF A(I).
*/

/*     START ELIMINATION LOOP */
LST = 0;
LEND = M - 1;
for( K=1; K<=M; K++ )
	{
/*     TEST ON USEFULNESS OF SYMMETRIC ALGORITHM */
	if( piv <= 0.0 )
		goto fatal;
	if( IER == 0 )
		{
		if( piv <= tol )
			{
			IER = K - 1;
			}
		}
	LT = J - K;
	LST += K;

/*  PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R */
	pivi = 1.0 / A[I-1];
	L = K;
	LL = L + LT;
	tb = pivi * R[LL-1];
	R[LL-1] = R[L-1];
	R[L-1] = tb;
/* IS ELIMINATION TERMINATED */
	if( K >= M )
		break;
/*
C     ROW AND COLUMN INTERCHANGE AND PIVOT ROW REDUCTION IN MATRIX A.
C     ELEMENTS OF PIVOT COLUMN ARE SAVED IN AUXILIARY VECTOR AUX.
*/
	LR = LST + (LT*(K+J-1))/2;
	LL = LR;
	L=LST;
	for( II=K; II<=LEND; II++ )
		{
		L += II;
		LL += 1;
		if( L == LR )
			{
			A[LL-1] = A[LST-1];
			tb = A[L-1];
			goto lab13;
			}
		if( L > LR )
			LL = L + LT;

		tb = A[LL-1];
		A[LL-1] = A[L-1];
lab13:
		AUX[II-1] = tb;
		A[L-1] = pivi * tb;
		}
/* SAVE COLUMN INTERCHANGE INFORMATION */
	A[LST-1] = LT;
/* ELEMENT REDUCTION AND SEARCH FOR NEXT PIVOT */
	piv = 0.0;
	LLST = LST;
	LT = 0;
	for( II=K; II<=LEND; II++ )
		{
		pivi = -AUX[II-1];
		LL = LLST;
		LT += 1;
		for( LLD=II; LLD<=LEND; LLD++ )
			{
			LL += LLD;
			L = LL + LT;
			A[L-1] += pivi * A[LL-1];
			}
		LLST += II;
		LR = LLST + LT;
		tb =fabs( A[LR-1] );
		if( tb > piv )
			{
			piv = tb;
			I = LR;
			J = II + 1;
			}
		LR = K;
		LL = LR + LT;
		R[LL-1] += pivi * R[LR-1];
		}
	}
/* END OF ELIMINATION LOOP */

/* BACK SUBSTITUTION AND BACK INTERCHANGE */

if( LEND <= 0 )
	{
	if( LEND < 0 )
		goto fatal;
	goto done;
	}
II = M;
for( I=2; I<=M; I++ )
	{
	LST -= II;
	II -= 1;
	L = A[LST-1] + 0.5;
	J = II;
	tb = R[J-1];
	LL = J;
	K = LST;
	for( LT=II; LT<=LEND; LT++ )
		{
		LL += 1;
		K += LT;
		tb -= A[K-1] * R[LL-1];
		}
	K = J + L;
	R[J-1] = R[K-1];
	R[K-1] = tb;
	}
done:
return( IER );
}
