/*							simpsn.c	*/
/* simpsn.c
 * Numerical integration of function tabulated
 * at equally spaced arguments
 */

/* Coefficients for Cote integration formulas */

/* Note: these numbers were computed using 40-decimal precision. */

#define NCOTE 8

/* 6th order formula */
/*
static double simcon[] =
{
  4.88095238095238095E-2,
  2.57142857142857142857E-1,
  3.2142857142857142857E-2,
  3.2380952380952380952E-1,
};
*/

/* 8th order formula */
static double simcon[] =
{
  3.488536155202821869E-2,
  2.076895943562610229E-1,
 -3.27336860670194003527E-2,
  3.7022927689594356261E-1,
 -1.6014109347442680776E-1,
};

/* 10th order formula */
/*
static double simcon[] =
{
  2.68341483619261397039E-2,
  1.77535941424830313719E-1,
 -8.1043570626903960237E-2,
  4.5494628827962161295E-1,
 -4.3515512265512265512E-1,
  7.1376463043129709796E-1,
};
*/

/*							simpsn.c 2	*/
/* 20th order formula */
/*
static double simcon[] =
{
  1.182527324903160319E-2,
  1.14137717644606974987E-1,
 -2.36478370511426964E-1,
  1.20618689348187566E+0,
 -3.7710317267153304677E+0,
  1.03367982199398011435E+1,
 -2.270881584397951229796E+1,
  4.1828057422193554603E+1,
 -6.4075279490154004651555E+1,
  8.279728347247285172085E+1,
 -9.0005367135242894657916E+1,
};
*/

/*							simpsn.c 3	*/

double simpsn( double [], double );

double simpsn( f, delta )
double f[];	/* tabulated function */
double delta;	/* spacing of arguments */
{
extern double simcon[];
double ans;
int i;


ans = simcon[NCOTE/2] * f[NCOTE/2];
for( i=0; i < NCOTE/2; i++ )
	ans += simcon[i] * ( f[i] + f[NCOTE-i] );

return( ans * delta * NCOTE );
}
