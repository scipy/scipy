#include <stdio.h>

void mysub(int n, double *x, double *y)
{
    return;
}

main()
{
    /* Parameters */    
#define NMAX    100
#define ITS     10000
    
    int      i, j;
    double   alpha, avg, t1, t2, tnotim;
    double   x[NMAX], y[NMAX];
    double   SuperLU_timer_();

    /* Initialize X and Y */
    for (i = 0; i < NMAX; ++i) {
	x[i] = 1.0 / (double)(i+1);
	y[i] = (double)(NMAX - i) / (double)NMAX;
    }
    alpha = 0.315;

    /* Time 1,000,000 DAXPY operations */
    t1 = SuperLU_timer_();
    for (j = 0; j < ITS; ++j) {
	for (i = 0; i < NMAX; ++i)
	    y[i] += alpha * x[i];
	alpha = -alpha;
    }
    t2 = SuperLU_timer_();
    printf("Time for 1,000,000 DAXPY ops  = %10.3g seconds\n", t2-t1);
    if ( t2-t1 > 0. ) 
	printf("DAXPY performance rate        = %10.3g mflops\n", 2./(t2-t1));
    else
	printf("*** Error:  Time for operations was zero\n");
	
    tnotim = t2 - t1;

    /* Time 1,000,000 DAXPY operations with SuperLU_timer_() 
       in the outer loop */
    t1 = SuperLU_timer_();
    for (j = 0; j < ITS; ++j) {
	for (i = 0; i < NMAX; ++i)
	    y[i] += alpha * x[i];
	alpha = -alpha;
	t2 = SuperLU_timer_();
    }

    /* Compute the time in milliseconds used by an average call to 
       SuperLU_timer_(). */
    printf("Including DSECND, time        = %10.3g seconds\n", t2-t1);
    avg = ( (t2 - t1) - tnotim )*1000. / (double)ITS;
    printf("Average time for DSECND       = %10.3g milliseconds\n", avg);

    /* Compute the equivalent number of floating point operations used
       by an average call to DSECND.    */
    if ( tnotim > 0. )
	printf("Equivalent floating point ops = %10.3g ops\n",
	       1000.*avg / tnotim);

    mysub(NMAX, x, y);
    return 0;
}

