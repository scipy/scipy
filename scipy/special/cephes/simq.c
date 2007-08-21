/*							simq.c
 *
 *	Solution of simultaneous linear equations AX = B
 *	by Gaussian elimination with partial pivoting
 *
 *
 *
 * SYNOPSIS:
 *
 * double A[n*n], B[n], X[n];
 * int n, flag;
 * int IPS[];
 * int simq();
 *
 * ercode = simq( A, B, X, n, flag, IPS );
 *
 *
 *
 * DESCRIPTION:
 *
 * B, X, IPS are vectors of length n.
 * A is an n x n matrix (i.e., a vector of length n*n),
 * stored row-wise: that is, A(i,j) = A[ij],
 * where ij = i*n + j, which is the transpose of the normal
 * column-wise storage.
 *
 * The contents of matrix A are destroyed.
 *
 * Set flag=0 to solve.
 * Set flag=-1 to do a new back substitution for different B vector
 * using the same A matrix previously reduced when flag=0.
 *
 * The routine returns nonzero on error; messages are printed.
 *
 *
 * ACCURACY:
 *
 * Depends on the conditioning (range of eigenvalues) of matrix A.
 *
 *
 * REFERENCE:
 *
 * Computer Solution of Linear Algebraic Systems,
 * by George E. Forsythe and Cleve B. Moler; Prentice-Hall, 1967.
 *
 */

/*							simq	2 */

#include <stdio.h>
#define ANSIPROT
#ifdef ANSIPROT
int simq(double [], double [], double [], int, int, int [] );
#endif

#define fabs(x) ((x) < 0 ? -(x) : (x))

int simq( A, B, X, n, flag, IPS )
double A[], B[], X[];
int n, flag;
int IPS[];
{
int i, j, ij, ip, ipj, ipk, ipn;
int idxpiv, iback;
int k, kp, kp1, kpk, kpn;
int nip, nkp, nm1;
double em, q, rownrm, big, size, pivot, sum;

nm1 = n-1;
if( flag < 0 )
	goto solve;

/*	Initialize IPS and X	*/

ij=0;
for( i=0; i<n; i++ )
	{
	IPS[i] = i;
	rownrm = 0.0;
	for( j=0; j<n; j++ )
		{
		q = fabs( A[ij] );
		if( rownrm < q )
			rownrm = q;
		++ij;
		}
	if( rownrm == 0.0 )
		{
		puts("SIMQ ROWNRM=0");
		return(1);
		}
	X[i] = 1.0/rownrm;
	}

/*							simq	3 */
/*	Gaussian elimination with partial pivoting 	*/

for( k=0; k<nm1; k++ )
	{
	big= 0.0;
	idxpiv = 0;
	for( i=k; i<n; i++ )
		{
		ip = IPS[i];
		ipk = n*ip + k;
		size = fabs( A[ipk] ) * X[ip];
		if( size > big )
			{
			big = size;
			idxpiv = i;
			}
		}

	if( big == 0.0 )
		{
		puts( "SIMQ BIG=0" );
		return(2);
		}
	if( idxpiv != k )
		{
		j = IPS[k];
		IPS[k] = IPS[idxpiv];
		IPS[idxpiv] = j;
		}
	kp = IPS[k];
	kpk = n*kp + k;
	pivot = A[kpk];
	kp1 = k+1;
	for( i=kp1; i<n; i++ )
		{
		ip = IPS[i];
		ipk = n*ip + k;
		em = -A[ipk]/pivot;
		A[ipk] = -em;
		nip = n*ip;
		nkp = n*kp;
		for( j=kp1; j<n; j++ )
			{
			ipj = nip + j;
			A[ipj] = A[ipj] + em * A[nkp + j];
			}
		}
	}
kpn = n * IPS[n-1] + n - 1;	/* last element of IPS[n] th row */
if( A[kpn] == 0.0 )
	{
	puts( "SIMQ A[kpn]=0");
	return(3);
	}

/*							simq 4 */
/*	back substitution	*/

solve:
ip = IPS[0];
X[0] = B[ip];
for( i=1; i<n; i++ )
	{
	ip = IPS[i];
	ipj = n * ip;
	sum = 0.0;
	for( j=0; j<i; j++ )
		{
		sum += A[ipj] * X[j];
		++ipj;
		}
	X[i] = B[ip] - sum;
	}

ipn = n * IPS[n-1] + n - 1;
X[n-1] = X[n-1]/A[ipn];

for( iback=1; iback<n; iback++ )
	{
/* i goes (n-1),...,1	*/
	i = nm1 - iback;
	ip = IPS[i];
	nip = n*ip;
	sum = 0.0;
	for( j=i+1; j<n; j++ )
		sum += A[nip+j] * X[j];
	X[i] = (X[i] - sum)/A[nip+i];
	}
return(0);
}
