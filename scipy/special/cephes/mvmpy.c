/*							mvmpy.c
 *
 *	Matrix times vector
 *
 *
 *
 * SYNOPSIS:
 *
 * int r, c;
 * double A[r*c], V[c], Y[r];
 *
 * mvmpy( r, c, A, V, Y );
 *
 *
 *
 * DESCRIPTION:
 *
 *          c-1
 *          --
 * Y[j] =   >   A[j][k] V[k] ,  j = 1, ..., r
 *          --
 *          k=0
 *
 * Multiplies the r (rows) by c (columns) matrix A on the left
 * by column vector V of dimension c on the right
 * to produce a (column) vector Y output of dimension r.
 *
 *
 *
 *
 */


void
mvmpy( r, c, A, V, Y )
int r, c;
double *A, *V, *Y;
{
register double s;
double *pA, *pV, *pY;
int i, j;

pA = A;
pY = Y;
for( i=0; i<r; i++ )
	{
	pV = V;
	s = 0.0;
	for( j=0; j<c; j++ )
		{
		s += *pA++ * *pV++;
		}
	*pY++ = s;
	}
}

