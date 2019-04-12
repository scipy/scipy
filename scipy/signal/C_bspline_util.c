#include "Python.h"
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define NO_IMPORT_ARRAY
#include "numpy/arrayobject.h"

void compute_root_from_lambda(double, double *, double *);





#define CONJ(a) (~(a))
#define ABSQ(a) (__real__ (a*CONJ(a)))
#ifdef __GNUC__

/* Implement the following difference equation */
/* y[n] = a1 * x[n] + a2 * y[n-1]  */
/* with a given starting value loaded into the array */

void C_IIR_order1 (__complex__ float,__complex__ float,__complex__ float*,__complex__ float*,int,int,int);
void C_IIR_order2 (__complex__ float,__complex__ float,__complex__ float,__complex__ float*,__complex__ float*,int,int,int);
void C_IIR_order2_cascade (__complex__ float,__complex__ float,__complex__ float,__complex__ float,__complex__ float*,__complex__ float*,int,int,int);
int C_IIR_forback1(__complex__ float,__complex__ float,__complex__ float*,__complex__ float*,int,int,int,float);
void C_FIR_mirror_symmetric(__complex__ float*,__complex__ float*,int,__complex__ float*,int,int,int);
int C_separable_2Dconvolve_mirror(__complex__ float*,__complex__ float*,int,int,__complex__ float*,__complex__ float*,int,int,npy_intp*,npy_intp*);

void 
C_IIR_order1(__complex__ float a1, __complex__ float a2,
             __complex__ float *x, __complex__ float *y,
             int N, int stridex, int stridey)
{ 
    __complex__ float *yvec = y+stridey; 
    __complex__ float *xvec = x+stridex; 
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
C_IIR_order2(__complex__ float a1, __complex__ float a2, __complex__ float a3,
             __complex__ float *x, __complex__ float *y,
             int N, int stridex, int stridey)
{ 
    __complex__ float *yvec = y+2*stridey; 
    __complex__ float *xvec = x+2*stridex; 
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
C_IIR_order2_cascade(__complex__ float cs,
                     __complex__ float z1, __complex__ float z2,
                     __complex__ float y1_0,
                     __complex__ float *x, __complex__ float *yp,
                     int N, int stridex, int stridey)
{ 
    __complex__ float *yvec = yp+stridey; 
    __complex__ float *xvec = x+stridex; 
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
C_IIR_forback1(__complex__ float c0, __complex__ float z1,
               __complex__ float *x, __complex__ float *y,
               int N, int stridex, int stridey, float precision)
{
    __complex__ float *yp = NULL; 
    __complex__ float *xptr = x;
    __complex__ float yp0;
    __complex__ float powz1;  
    __complex__ float diff;
    float err;
    int k;

    if (ABSQ(z1) >= 1.0) return -2; /* z1 not less than 1 */

    /* Initialize memory for loop */ 
    if ((yp = malloc(N*sizeof(__complex__ float)))==NULL) return -1; 

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

    C_IIR_order1(1.0, z1, x, yp, N, stridex, 1); 

    *(y + (N-1)*stridey) = -c0 / (z1 - 1.0) * yp[N-1];

    C_IIR_order1(c0, z1, yp+N-1, y+(N-1)*stridey, N, -1, -stridey);

    free(yp);
    return 0;
}


/* h must be odd length */
/* strides in units of sizeof(__complex__ float) bytes */
void 
C_FIR_mirror_symmetric(__complex__ float *in, __complex__ float *out,
                       int N, __complex__ float *h, int Nh,
                       int instride, int outstride)
{    
    int n, k;
    int Nhdiv2 = Nh >> 1;
    __complex__ float *outptr;
    __complex__ float *inptr;
    __complex__ float *hptr;

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
C_separable_2Dconvolve_mirror(__complex__ float *in, __complex__ float *out,
                              int M, int N,
                              __complex__ float *hr, __complex__ float *hc,
                              int Nhr, int Nhc,
                              npy_intp *instrides, npy_intp *outstrides)
{
    int m, n;
    __complex__ float *tmpmem;
    __complex__ float *inptr=NULL, *outptr=NULL;
    
    tmpmem = malloc(M*N*sizeof(__complex__ float));
    if (tmpmem == NULL) return -1;

    if (Nhr > 0) {
	/* filter across rows */
	inptr = in;
	outptr = tmpmem;    
	for (m = 0; m < M; m++) {
	    C_FIR_mirror_symmetric (inptr, outptr, N, hr, Nhr, instrides[1], 1);
	    inptr += instrides[0];
	    outptr += N;
	}
    }
    else 
	memmove(tmpmem, inptr, M*N*sizeof(__complex__ float));
	
    if (Nhc > 0) {
	/* filter down columns */
	inptr = tmpmem;
	outptr = out;
	for (n = 0; n < N; n++) {
	    C_FIR_mirror_symmetric (inptr, outptr, M, hc, Nhc, N, outstrides[0]);
	    outptr += outstrides[1];
	    inptr += 1;
	}
    }
    else
	memmove(outptr, tmpmem, M*N*sizeof(__complex__ float));

    free(tmpmem);
    return 0;
}
#endif
