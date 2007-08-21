/*							mtransp.c
 *
 *	Matrix transpose
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double A[n*n], T[n*n];
 *
 * mtransp( n, A, T );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 * T[r][c] = A[c][r]
 *
 *
 * Transposes the n by n square matrix A and puts the result in T.
 * The output, T, may occupy the same storage as A.
 *
 *
 *
 */


#define ANSIPROT
#ifdef ANSIPROT
void mtransp( int,double*,double* );
#endif

void
mtransp( n, A, T )
int n;
double *A, *T;
{
int i, j, np1;
double *pAc, *pAr, *pTc, *pTr, *pA0, *pT0;
double x;

np1 = n+1;
pA0 = A;
pT0 = T;
for( i=0; i<n-1; i++ ) /* row index */
	{
	pAc = pA0; /* next diagonal element of input */
	pAr = pAc + n; /* next row down underneath the diagonal element */
	pTc = pT0; /* next diagonal element of the output */
	pTr = pTc + n; /* next row underneath */
	*pTc++ = *pAc++; /* copy the diagonal element */
	for( j=i+1; j<n; j++ ) /* column index */
		{
		x = *pAr;
		*pTr = *pAc++;
		*pTc++ = x;
		pAr += n;
		pTr += n;
		}
	pA0 += np1; /* &A[n*i+i] for next i */
	pT0 += np1; /* &T[n*i+i] for next i */
	}
*pT0 = *pA0; /* copy the diagonal element */
}

