/*		Levnsn.c		*/
/* Levinson-Durbin LPC
 *
 * | R0 R1 R2 ... RN-1 |   | A1 |       | -R1 |
 * | R1 R0 R1 ... RN-2 |   | A2 |       | -R2 |
 * | R2 R1 R0 ... RN-3 |   | A3 |   =   | -R3 |
 * |          ...      |   | ...|       | ... |
 * | RN-1 RN-2... R0   |   | AN |       | -RN |
 *
 * Ref: John Makhoul, "Linear Prediction, A Tutorial Review"
 * Proc. IEEE Vol. 63, PP 561-580 April, 1975.
 *
 * R is the input autocorrelation function.  R0 is the zero lag
 * term.  A is the output array of predictor coefficients.  Note
 * that a filter impulse response has a coefficient of 1.0 preceding
 * A1.  E is an array of mean square error for each prediction order
 * 1 to N.  REFL is an output array of the reflection coefficients.
 */

#define abs(x) ( (x) < 0 ? -(x) : (x) )

levnsn( n, r, a, e, refl )
int n;
double r[], a[], e[], refl[];
{
int k, km1, kp1, i, kmi, j;
double ai, akk, err, err1, r0, t, akmi;
double *pa, *pr;

for( i=0; i<n; i++ )
	{
	a[i] = 0.0;
	e[i] = 0.0;
	refl[i] = 0.0;
	}
r0 = r[0];
e[0] = r0;
err = r0;

akk = -r[1]/err;
err = (1.0 - akk*akk) * err;
e[1] = err;
a[1] = akk;
refl[1] = akk;

if( err < 1.0e-2 )
	return;

for( k=2; k<n; k++ )
	{
	t = 0.0;
	pa = &a[1];
	pr = &r[k-1];
	for( j=1; j<k; j++ )
		t += *pa++ * *pr--;
	akk = -( r[k] + t )/err;
	refl[k] = akk;
	km1 = k/2;
	for( j=1; j<=km1; j++ )
		{
		kmi = k-j;
		ai = a[j];
		akmi = a[kmi];
		a[j] = ai + akk*akmi;
		if( i == kmi )
			goto nxtk;
		a[kmi] = akmi + akk*ai;
		}
nxtk:
	a[k] = akk;
	err1 = (1.0 - akk*akk)*err;
	e[k] = err1;
	if( err1 < 0 )
		err1 = -err1;
/*	err1 = abs(err1);*/
/*	if( (err1 < 1.0e-2) || (err1 >= err) )*/
	if( err1 < 1.0e-2 )
		return;
	err = err1;
	}
}
