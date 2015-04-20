/*                                                     polrt.c
 *
 *     Find roots of a polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * typedef struct
 *     {
 *     double r;
 *     double i;
 *     }cmplx;
 *
 * double xcof[], cof[];
 * int m;
 * cmplx root[];
 *
 * polrt( xcof, cof, m, root )
 *
 *
 *
 * DESCRIPTION:
 *
 * Iterative determination of the roots of a polynomial of
 * degree m whose coefficient vector is xcof[].  The
 * coefficients are arranged in ascending order; i.e., the
 * coefficient of x**m is xcof[m].
 *
 * The array cof[] is working storage the same size as xcof[].
 * root[] is the output array containing the complex roots.
 *
 *
 * ACCURACY:
 *
 * Termination depends on evaluation of the polynomial at
 * the trial values of the roots.  The values of multiple roots
 * or of roots that are nearly equal may have poor relative
 * accuracy after the first root in the neighborhood has been
 * found.
 *
 */

/*                                                     polrt   */
/* Complex roots of real polynomial */
/* number of coefficients is m + 1 ( i.e., m is degree of polynomial) */

#include "mconf.h"
/*
 * typedef struct
 * {
 * double r;
 * double i;
 * }cmplx;
 */

int polrt(xcof, cof, m, root)
double xcof[], cof[];
int m;
cmplx root[];
{
    register double *p, *q;
    int i, j, nsav, n, n1, n2, nroot, iter, retry;
    int final;
    double mag, cofj;
    cmplx x0, x, xsav, dx, t, t1, u, ud;

    xsav.r = 0;
    xsav.i = 0;
    final = 0;
    n = m;
    if (n <= 0)
	return (1);
    if (n > 36)
	return (2);
    if (xcof[m] == 0.0)
	return (4);

    n1 = n;
    n2 = n;
    nroot = 0;
    nsav = n;
    q = &xcof[0];
    p = &cof[n];
    for (j = 0; j <= nsav; j++)
	*p-- = *q++;		/*      cof[ n-j ] = xcof[j]; */

  nxtrut:
    x0.r = 0.00500101;
    x0.i = 0.01000101;
    retry = 0;

  tryagn:
    retry += 1;
    x.r = x0.r;

    x0.r = -10.0 * x0.i;
    x0.i = -10.0 * x.r;

    x.r = x0.r;
    x.i = x0.i;

  finitr:
    iter = 0;

    while (iter < 500) {
	u.r = cof[n];
	if (u.r == 0.0) {	/* this root is zero */
	    x.r = 0;
	    n1 -= 1;
	    n2 -= 1;
	    goto zerrut;
	}
	u.i = 0;
	ud.r = 0;
	ud.i = 0;
	t.r = 1.0;
	t.i = 0;
	p = &cof[n - 1];
	for (i = 0; i < n; i++) {
	    t1.r = x.r * t.r - x.i * t.i;
	    t1.i = x.r * t.i + x.i * t.r;
	    cofj = *p--;	/* evaluate polynomial */
	    u.r += cofj * t1.r;
	    u.i += cofj * t1.i;
	    cofj = cofj * (i + 1);	/* derivative */
	    ud.r += cofj * t.r;
	    ud.i -= cofj * t.i;
	    t.r = t1.r;
	    t.i = t1.i;
	}

	mag = ud.r * ud.r + ud.i * ud.i;
	if (mag == 0.0) {
	    if (!final)
		goto tryagn;
	    x.r = xsav.r;
	    x.i = xsav.i;
	    goto findon;
	}
	dx.r = (u.i * ud.i - u.r * ud.r) / mag;
	x.r += dx.r;
	dx.i = -(u.r * ud.i + u.i * ud.r) / mag;
	x.i += dx.i;
	if ((fabs(dx.i) + fabs(dx.r)) < 1.0e-6)
	    goto lupdon;
	iter += 1;
    }				/* while iter < 500 */

    if (final)
	goto lupdon;
    if (retry < 5)
	goto tryagn;
    return (3);

  lupdon:
    /* Swap original and reduced polynomials */
    q = &xcof[nsav];
    p = &cof[0];
    for (j = 0; j <= n2; j++) {
	cofj = *q;
	*q-- = *p;
	*p++ = cofj;
    }
    i = n;
    n = n1;
    n1 = i;

    if (!final) {
	final = 1;
	if (fabs(x.i / x.r) < 1.0e-4)
	    x.i = 0.0;
	xsav.r = x.r;
	xsav.i = x.i;
	goto finitr;		/* do final iteration on original polynomial */
    }

  findon:
    final = 0;
    if (fabs(x.i / x.r) >= 1.0e-5) {
	cofj = x.r + x.r;
	mag = x.r * x.r + x.i * x.i;
	n -= 2;
    }
    else {			/* root is real */
      zerrut:
	x.i = 0;
	cofj = x.r;
	mag = 0;
	n -= 1;
    }
    /* divide working polynomial cof(z) by z - x */
    p = &cof[1];
    *p += cofj * *(p - 1);
    for (j = 1; j < n; j++) {
	*(p + 1) += cofj * *p - mag * *(p - 1);
	p++;
    }

  setrut:
    root[nroot].r = x.r;
    root[nroot].i = x.i;
    nroot += 1;
    if (mag != 0.0) {
	x.i = -x.i;
	mag = 0;
	goto setrut;		/* fill in the complex conjugate root */
    }
    if (n > 0)
	goto nxtrut;
    return (0);
}
