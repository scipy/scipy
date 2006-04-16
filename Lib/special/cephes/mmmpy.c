/*							mmmpy.c
 *
 *	Matrix multiply
 *
 *
 *
 * SYNOPSIS:
 *
 * int r, c;
 * double A[r*c], B[c*r], Y[r*r];
 *
 * mmmpy( r, c, A, B, Y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Y = A B
 *              c-1
 *              --
 * Y[i][j]  =   >   A[i][k] B[k][j]
 *              --
 *              k=0
 *
 * Multiplies an r (rows) by c (columns) matrix A on the left
 * by a c (rows) by r (columns) matrix B on the right
 * to produce an r by r matrix Y.
 *
 *
 */


#define ANSIPROT
#ifdef ANSIPROT
void mmmpy( int,int,double*,double*,double* );
#endif

void
mmmpy( r, c, A, B, Y )
int r, c;
double *A, *B, *Y;
{
register double s;
double *pA, *pB, *pY, *pt;
int i, j, k;

pY = Y;
pB = B;
for( i=0; i<r; i++ )
	{
	pA = A;
	for( j=0; j<r; j++ )
		{
		pt = pB;
		s = 0.0;
		for( k=0; k<c; k++ )
			{
			s += *pA++ * *pt;
			pt += r; /* increment to next row underneath */
			}
		*pY++ = s;
		}
	pB += 1;
	}
}

