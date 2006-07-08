#include "Python.h"
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define NO_IMPORT_ARRAY
#include "numpy/arrayobject.h"

void compute_root_from_lambda(double, double *, double *);

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */    
#endif

#define CONJ(a) ((a))
#define ABSQ(a) ( (a*CONJ(a)))

void D_IIR_order1(double,double,double*,double*,int,int,int); 
void D_IIR_order2(double,double,double,double*,double*,int,int,int);
void D_IIR_order2_cascade(double,double,double,double,double*,double*,int,int,int);
int D_IIR_forback1(double,double,double*,double*,int,int,int,double);
void D_FIR_mirror_symmetric(double*,double*,int,double*,int,int,int);
int D_separable_2Dconvolve_mirror(double*,double*,int,int,double*,double*,int,int,npy_intp*,npy_intp*);
int D_IIR_forback2(double,double,double*,double*,int,int,int,double); 
int D_cubic_spline2D(double*,double*,int,int,double,npy_intp*,npy_intp*,double);
int D_quadratic_spline2D(double*,double*,int,int,double,npy_intp*,npy_intp*,double);

/* Implement the following difference equation */
/* y[n] = a1 * x[n] + a2 * y[n-1]  */
/* with a given starting value loaded into the array */

void 
D_IIR_order1 (a1, a2, x, y, N, stridex, stridey) 
     double a1; 
     double a2; 
     double *x; 
     double *y; 
     int N, stridex, stridey; 
{ 
    double *yvec = y+stridey; 
    double *xvec = x+stridex; 
    int n; 

    for (n=1; n < N; n++) { 
	*yvec = *xvec * a1 + *(yvec-stridey) * a2; 
	yvec += stridey; 
	xvec += stridex; 
    } 
}


/* Implement the following difference equation */
/* y[n] = a1 * x[n] + a2 * y[n-1]  + a3 * y[n-2] */
/* with two starting values loaded into the array */
void 
D_IIR_order2 (a1, a2, a3, x, y, N, stridex, stridey) 
     double a1; 
     double a2; 
     double a3; 
     double *x; 
     double *y; 
     int N, stridex, stridey; 
{ 
    double *yvec = y+2*stridey; 
    double *xvec = x+2*stridex; 
    int n; 

    for (n=2; n < N; n++) { 
	*yvec = *xvec * a1 + *(yvec-stridey) * a2 + *(yvec-2*stridey) * a3; 
	yvec += stridey; 
	xvec += stridex;
    } 
}

/* Implement a second order IIR difference equation using a cascade
   of first order sections.  Suppose the transfer function is 
                  cs   
   H(z) =   -------------------
            (1-z1/z) ( 1-z2/z)                    

   then the following pair is implemented with one starting value loaded in
   the output array and the starting value for the intermediate array
   passed in as yp0.

   y1[n] = x[n] + z1 y1[n-1]
   yp[n] = cs y1[n] + z2 yp[n-1]

*/

void 
D_IIR_order2_cascade (cs, z1, z2, y1_0, x, yp, N, stridex, stridey) 
     double cs; 
     double z1; 
     double z2; 
     double y1_0; 
     double *x; 
     double *yp; 
     int N, stridex, stridey; 
{ 
    double *yvec = yp+stridey; 
    double *xvec = x+stridex; 
    int n; 

    for (n=1; n < N; n++) { 
	y1_0 = *xvec + y1_0 * z1; 
	*yvec = cs * y1_0 + *(yvec-stridey) * z2; 
	yvec += stridey; 
	xvec += stridex;  
    } 
}


/* Implement a smoothing IIR filter with mirror-symmetric boundary conditions
   using a cascade of first-order sections.  The second section uses a 
   reversed sequence.  This implements the following transfer function:
                    c0
   H(z) = ---------------------------
           (1-z1/z) (1 - z1 z)

   with the following difference equations:

   yp[n] = x[n] + z1 yp[n-1]  
     with starting condition: 
   yp[0] = x[0] + Sum(z1^(k+1) x[k],k=0..Infinity)

   and 

   y[n] = z1 y[n+1] + c0 yp[n]
     with starting condition:
   y[N-1] = z1 / (z1-1) yp[N-1] 

   The resulting signal will have mirror symmetric boundary conditions as well.

   If memory could not be allocated for the temporary vector yp, the 
   function returns -1 otherwise it returns 0.
   
   z1 should be less than 1;
   
*/

int 
D_IIR_forback1 (c0, z1, x, y, N, stridex, stridey, precision)
     double c0; 
     double z1; 
     double *x; 
     double *y; 
     int N, stridex, stridey; 
     double precision; 
{ 
    double *yp = NULL; 
    double *xptr = x;
    double yp0;
    double powz1;  
    double diff;
    double err;
    int k;

    if (ABSQ(z1) >= 1.0) return -2; /* z1 not less than 1 */

    /* Initialize memory for loop */ 
    if ((yp = malloc(N*sizeof(double)))==NULL) return -1; 

   /* Fix starting value assuming mirror-symmetric boundary conditions. */
    yp0 = x[0];
    powz1 = 1.0;
    k = 0;
    precision *= precision;
    do {
	yp[0] = yp0;
	powz1 *= z1;
	yp0 += powz1 * (*xptr);
	diff = powz1;
	err = ABSQ(diff);
	xptr += stridex;
	k++;
    } while((err > precision) && (k < N));
    if (k >= N) return -3;     /* sum did not converge */ 
    yp[0] = yp0;

    D_IIR_order1(1.0, z1, x, yp, N, stridex, 1); 

    *(y + (N-1)*stridey) = -c0 / (z1 - 1.0) * yp[N-1];

    D_IIR_order1(c0, z1, yp+N-1, y+(N-1)*stridey, N, -1, -stridey);

    free(yp);
    return 0;
}


/* h must be odd length */
/* strides in units of sizeof(double) bytes */
void 
D_FIR_mirror_symmetric (in, out, N, h, Nh, instride, outstride)
     double *in;
     double *out;
     int N, Nh;
     double *h;
     int instride, outstride;
{    
    int n, k;
    int Nhdiv2 = Nh >> 1;
    double *outptr;
    double *inptr;
    double *hptr;

    /* first part boundary conditions */
    outptr = out;
    for (n=0; n < Nhdiv2; n++) {
	*outptr = 0.0;
	hptr = h;
	inptr = in + (n+Nhdiv2)*instride;
	for (k=-Nhdiv2; k <= n; k++) {
	    *outptr += *hptr++ * *inptr;
	    inptr -= instride;
	}
	inptr += instride;
	for (k=n+1; k <= Nhdiv2; k++) {
	    *outptr += *hptr++ * *inptr;
	    inptr += instride;
	}	
	outptr += outstride;
    }
    
    /* middle section */
    outptr = out + Nhdiv2*outstride;
    for (n=Nhdiv2; n < N-Nhdiv2; n++) {
	*outptr = 0.0;
	hptr = h;
	inptr = in + (n+Nhdiv2)*instride;
	for (k=-Nhdiv2; k <= Nhdiv2; k++) {
	    *outptr += *hptr++ * *inptr;
	    inptr -= instride;
	}
	outptr += outstride;
    }

    /* end boundary conditions */
    outptr = out + (N-Nhdiv2)*outstride;
    for (n=N-Nhdiv2; n < N; n++) {
	*outptr = 0.0;
	hptr = h;
	inptr = in + (2*N-1-n-Nhdiv2)*instride;
	for (k=-Nhdiv2; k <= n-N; k++) {
	    *outptr += *hptr++ * *inptr;
	    inptr += instride;
	}
	inptr -= instride;
	for (k=n+1-N; k <= Nhdiv2; k++) {
	    *outptr += *hptr++ * *inptr;
	    inptr -= instride;
	}
	outptr += outstride;
    }
   
}

int
D_separable_2Dconvolve_mirror(in, out, M, N, hr, hc, Nhr, 
				   Nhc, instrides, outstrides)
     double *in;
     double *out;
     int M, N;
     double *hr, *hc;
     int Nhr, Nhc;
     npy_intp *instrides, *outstrides;
{
    int m, n;
    double *tmpmem;
    double *inptr=0, *outptr=0;
    
    tmpmem = malloc(M*N*sizeof(double));
    if (tmpmem == NULL) return -1;

    if (Nhr > 0) {
	/* filter across rows */
	inptr = in;
	outptr = tmpmem;    
	for (m = 0; m < M; m++) {
	    D_FIR_mirror_symmetric (inptr, outptr, N, hr, Nhr, instrides[1], 1);
	    inptr += instrides[0];
	    outptr += N;
	}
    }
    else 
	memmove(tmpmem, inptr, M*N*sizeof(double));
	
    if (Nhc > 0) {
	/* filter down columns */
	inptr = tmpmem;
	outptr = out;
	for (n = 0; n < N; n++) {
	    D_FIR_mirror_symmetric (inptr, outptr, M, hc, Nhc, N, outstrides[0]);
	    outptr += outstrides[1];
	    inptr += 1;
	}
    }
    else
	memmove(outptr, tmpmem, M*N*sizeof(double));

    free(tmpmem);
    return 0;
}


static double D_hc(int,double,double,double);
static double D_hs(int,double,double,double);

double
D_hc(k, cs, r, omega)
     int k;
     double cs;
     double r, omega;
{
    if (k < 0) return 0.0;
    if (omega == 0.0) 
	return cs * pow(r, (double )k) * (k+1); 
    else if (omega == M_PI) 
	return cs * pow(r, (double )k) * (k+1) * (1 - 2*(k % 2));
    return cs * pow(r, (double) k) * sin(omega * (k+1)) / sin(omega);
}

double
D_hs(k, cs, rsq, omega)
     int k;
     double cs;
     double rsq, omega;
{
    double cssq;
    double c0;
    double gamma, rsupk;

    cssq = cs * cs;
    k = abs(k);
    rsupk = pow(rsq, ((double ) k) / 2.0);
    if (omega == 0.0) {
	c0 = (1+rsq)/ ((1-rsq)*(1-rsq)*(1-rsq)) * cssq;
	gamma = (1-rsq) / (1+rsq);
	return c0 * rsupk * (1 + gamma * k);
    }
    if (omega == M_PI) {
	c0 = (1+rsq)/ ((1-rsq)*(1-rsq)*(1-rsq)) * cssq;
	gamma = (1-rsq) / (1+rsq) * (1 - 2 * (k % 2));
	return c0 * rsupk * (1 + gamma * k);
    }
    c0 = cssq * (1.0+rsq)/(1.0-rsq) / (1-2*rsq*cos(2*omega) + rsq*rsq);
    gamma = (1.0 - rsq)/ (1.0+rsq) / tan(omega);
    return c0 * rsupk * (cos(omega*k) + gamma * sin(omega * k));
}


/* Implement a smoothing IIR filter with mirror-symmetric boundary conditions
   using a cascade of second-order sections.  The second section uses a 
   reversed sequence.  This implements the following transfer function:

                          cs^2  
   H(z) = --------------------------------------
          (1 - a2/z - a3/z^2) (1 - a2 z -a3 z^2 )

   where a2 = (2 r cos omega)
         a3 = - r^2 
	 cs = 1 - 2 r cos omega + r^2 

   with the following difference equations:

   yp[n] = cs*x[n] - b1 yp[n-1] - b2 yp[n-2] 
     with starting conditions: 
   yp[0] = hc[0] x[0] + Sum(hc[k+1]*x[k],k=0..Infinity)
   yp[1] = hc[0] x[1] + hc[1] x[0] + Sum(hc[k+2] x[k], k=0..Infinity)

   and 

   y[n] = cs*yp[n] - b1 y[n+1] -b2 y[n+2]
     with starting conditions:
   y[N-1] = Sum((hs[k] + hs[k+1])x[N-1-k],k=0..Infinity)
   y[N-2] = Sum((hs[k-1] + hs[k+2])x[N-1-k],k=0..Infinity)

   The resulting signal will have mirror symmetric boundary conditions as well.

   If memory could not be allocated for the temporary vector yp, the 
   function returns -1 otherwise it returns 0.
   
   z1 should be less than 1;
   
*/

int 
D_IIR_forback2 (r, omega, x, y, N, stridex, stridey, precision) 
     double r,omega;
     double *x;
     double *y; 
     int N, stridex, stridey; 
     double precision; 
{ 
    double cs;
    double *yp = NULL; 
    double *yptr;
    double *xptr;
    double yp0;
    double yp1;
    double rsq;
    double diff;
    double err; 
    double a2, a3;
    int k;

    if (r >= 1.0) return -2; /* z1 not less than 1 */

    /* Initialize memory for loop */ 
    if ((yp = malloc(N*sizeof(double)))==NULL) return -1; 
    
    rsq = r * r;
    a2 = 2 * r * cos(omega);
    a3 = -rsq;
    cs = 1 - 2 * r * cos(omega) + rsq;

   /* Fix starting values assuming mirror-symmetric boundary conditions. */
    yp0 = D_hc(0, cs, r, omega) * x[0];
    k = 0;
    precision *= precision;
    xptr = x;
    do {
	yp[0] = yp0;
	diff = D_hc(k+1, cs, r, omega); 
  	yp0 += diff * (*xptr); 
	err = diff * diff;
	xptr += stridex;
	k++;
    } while((err > precision) && (k < N));
    if (k >= N) {free(yp); return -3;}     /* sum did not converge */ 
    yp[0] = yp0;

    yp1 = D_hc(0, cs, r, omega) * (*(x+stridex));
    yp1 += D_hc(1, cs, r, omega) * x[0];
    k = 0;
    xptr = x;
    do {
	yp[1] = yp1;
	diff = D_hc(k+2, cs, r, omega);
	yp1 += diff * (*xptr);
	err = diff * diff;
	xptr += stridex;
	k++;
    } while((err > precision) && (k < N));
    if (k >= N) {free(yp); return -3;}     /* sum did not converge */ 
    yp[1] = yp1;

    D_IIR_order2(cs, a2, a3, x, yp, N, stridex, 1); 

   /* Fix starting values assuming mirror-symmetric boundary conditions. */
    yp0 = 0.0;
    k = 0;
    yptr = y + (N-1)*stridey;
    xptr = x + (N-1)*stridex;
    do {
	*yptr = yp0;
	diff = (D_hs(k, cs, rsq, omega) + D_hs(k+1, cs, rsq, omega));
	yp0 += diff * (*xptr);
	err = diff * diff;
	xptr -= stridex;
	k++;
    } while((err > precision) && (k < N));
    if (k >= N) {free(yp); return -3;}     /* sum did not converge */ 
    *yptr = yp0;

    yp1 = 0.0;
    k = 0;
    yptr -= stridey;        /* Initialize in next-to-last slot in output array */
    xptr = x + (N-1)*stridex;
    do {
	*yptr = yp1;
	diff = (D_hs(k-1, cs, rsq, omega) + D_hs(k+2, cs, rsq, omega));
	yp1 += diff * (*xptr);
	err = diff * diff;
	xptr -= stridex;
	k++;
    } while((err > precision) && (k < N));
    if (k >= N) {free(yp); return -3;}     /* sum did not converge */ 
    *yptr = yp1;

    D_IIR_order2(cs, a2, a3, yp+N-1, yptr+stridey, N, -1, -stridey);

    free(yp);
    return 0;
}

/* Find the cubic spline coefficients of an image 
   image is M rows by N columns stored rowise in memory (vary column number
          first). It will be replaced with the spline coefficients. 
   lambda is a smoothing parameter (lambda = 100 approximately corresponds
          to a cutoff frequency of 0.1*(sample freq))
   strides is an integer array [rowstride, colstride]
          telling how much memory in units of sizeof(double) bytes to skip
	  to get to the next element.
*/

/* to get the (smoothed) image back mirror-symmetric convolve with a length 
   three separable FIR filter [1.0, 4.0, 1.0]/ 6.0
*/

int 
D_cubic_spline2D(image, coeffs, M, N, lambda, strides, cstrides, precision)
     double *image;
     double *coeffs;
     int M, N;
     double lambda;
     npy_intp *strides, *cstrides; 
     double precision;
{    
    double r, omega;
    double *inptr;
    double *coptr;
    double *tmpmem;
    double *tptr;
    int m,n, retval=0;

    tmpmem = malloc(N*M*sizeof(double));
    if (tmpmem == NULL) return -1;

    if (lambda <= 1.0 / 144.0) { 
	/* normal cubic spline */	
	r = -2 + sqrt(3.0);

	/* Loop over rows */
	inptr = image;
	tptr = tmpmem;    
	for (m = 0; m < M; m++) { 
	    retval = D_IIR_forback1 (-r*6.0, r, inptr, tptr, N, strides[1], 1, precision);
	    if (retval < 0) break;
	    inptr += strides[0];   
	    tptr += N;
	}

	if (retval >=0) {
	    /* Loop over columns */
	    tptr = tmpmem;
	    coptr = coeffs;
	    for (n = 0; n < N; n++) {
		retval = D_IIR_forback1 (-r*6.0, r, tptr, coptr, M, N, cstrides[0], precision);
		if (retval < 0) break;
		coptr += cstrides[1];
		tptr += 1;
	    }
	}
	free(tmpmem);
	return retval;
    }

    /* Smoothing spline */

    /* Compute r and omega from lambda */
    compute_root_from_lambda(lambda, &r, &omega);

    /* Loop over rows */
    inptr = image;
    tptr = tmpmem;    
    for (m = 0; m < M; m++) {
	retval = D_IIR_forback2 (r, omega, inptr, tptr, N, strides[1], 
				       1, precision);
	if (retval < 0) break;
	inptr += strides[0];   
	tptr += N;
    }
    /* Loop over columns */
    tptr = tmpmem;
    coptr = coeffs;
    for (n = 0; n < N; n++) {
	retval = D_IIR_forback2 (r, omega, tptr, coptr, M, N, 
				       cstrides[0], precision);
	if (retval < 0) break;
	coptr += cstrides[1];
	tptr += 1;
    }

    free(tmpmem);
    return retval;
}

/* Find the quadratic spline coefficients of an image 
   image is M rows by N columns stored rowise in memory (vary column number
          first). It will be replaced with the spline coefficients. 
   lambda is a smoothing parameter (lambda = 100 approximately corresponds
          to a cutoff frequency of 0.1*(sample freq))
	  must be zero for now.
   strides is an integer array [rowstride, colstride]
          telling how much memory in units of sizeof(double) bytes to skip
	  to get to the next element.
*/

/* to get the (smoothed) image back mirror-symmetric convolve with a length 
   three separable FIR filter [1.0, 6.0, 1.0]/ 8.0
*/

int 
D_quadratic_spline2D(image, coeffs, M, N, lambda, strides, cstrides, precision)
     double *image;
     double *coeffs;
     int M, N;
     double lambda;
     npy_intp *strides, *cstrides; 
     double precision;
{    
    double r;
    double *inptr;
    double *coptr;
    double *tmpmem;
    double *tptr;
    int m,n, retval=0;

    tmpmem = malloc(N*M*sizeof(double));
    if (tmpmem == NULL) return -1;

    if (lambda > 0) return -2;
    /* normal quadratic spline */	
    r = -3 + 2*sqrt(2.0);
    
    /* Loop over rows */
    inptr = image;
    tptr = tmpmem;    
    for (m = 0; m < M; m++) { 
      retval = D_IIR_forback1 (-r*8.0, r, inptr, tptr, N, strides[1], 1, precision);
      if (retval < 0) break;
      inptr += strides[0];   
      tptr += N;
    }
    
    if (retval >=0) {
    /* Loop over columns */
      tptr = tmpmem;
      coptr = coeffs;
      for (n = 0; n < N; n++) {
	retval = D_IIR_forback1 (-r*8.0, r, tptr, coptr, M, N, cstrides[0], precision);
	if (retval < 0) break;
	coptr += cstrides[1];
	tptr += 1;
      }
    }
    free(tmpmem);
    return retval;
}





