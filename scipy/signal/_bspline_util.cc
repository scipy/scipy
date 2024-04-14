#include "Python.h"
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "_splinemodule.h"
#include "_bspline_util.h"

#define NO_IMPORT_ARRAY
#include "numpy/arrayobject.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif



#ifdef __GNUC__
int C_SYM_IIR1_initial(__complex__ float z1, __complex__ float *x, __complex__ float *yp0,
                             int M, int N, float precision)
{
    return _sym_iir1_initial(z1, x, yp0, M, N, precision);

}

int Z_SYM_IIR1_initial(__complex__ double z1, __complex__ double *x, __complex__ double *yp0,
                             int M, int N, double precision)
{
    return _sym_iir1_initial(z1, x, yp0, M, N, precision);
}
#endif

int S_SYM_IIR1_initial(float z1, float *x, float *yp0,
                             int M, int N, float precision)
{
    return _sym_iir1_initial(z1, x, yp0, M, N, precision);
}

int D_SYM_IIR1_initial(double z1, double *x, double *yp0,
                             int M, int N, double precision)
{
    return _sym_iir1_initial(z1, x, yp0, M, N, precision);
}




/**
Compute the starting initial conditions for the system
cs / (1 - a2 * z^-1 - a3 * z^-2) against signals x
**/
int S_SYM_IIR2_initial_fwd(double r, double omega,
        float *x, float *yp, int M, int N, float precision)
{
    return _sym_iir2_initial_fwd(r, omega, x, yp, M, N, precision);
}

int D_SYM_IIR2_initial_fwd(double r, double omega,
        double *x, double *yp, int M, int N, double precision) {
    return _sym_iir2_initial_fwd(r, omega, x, yp, M, N, precision);
}

/*
Compute the starting initial conditions for the system (ran in backwards)
cs / (1 - a2 * z - a3 * z^2) against signal x
*/

int S_SYM_IIR2_initial_bwd(double r, double omega,
        float *x, float *yp, int M, int N, float precision)
{
    return _sym_iir2_initial_bwd(r, omega, x, yp, M, N, precision);

}

int D_SYM_IIR2_initial_bwd(double r, double omega,
        double *x, double *yp, int M, int N, double precision)
{
    return _sym_iir2_initial_bwd(r, omega, x, yp, M, N, precision);
}



// _FIR_mirros_symmetric instantiations //


/* h must be odd length */
/* strides in units of sizeof(DATA TYPE) bytes */

void S_FIR_mirror_symmetric(float *in, float *out, int N,
                            float *h, int Nh, int instride, int outstride)
{
    return _fir_mirror_symmetric(in, out, N, h, Nh, instride, outstride);
}

void D_FIR_mirror_symmetric(double *in, double *out, int N,
                            double *h, int Nh, int instride, int outstride)
{
    return _fir_mirror_symmetric(in, out, N, h, Nh, instride, outstride);
}


#ifdef __GNUC__
void C_FIR_mirror_symmetric(__complex__ float *in, __complex__ float *out, int N,
                            __complex__ float *h, int Nh, int instride, int outstride)
{
    return _fir_mirror_symmetric(in, out, N, h, Nh, instride, outstride);
}

void Z_FIR_mirror_symmetric(__complex__ double *in, __complex__ double *out, int N,
                            __complex__ double *h, int Nh, int instride, int outstride)
{
    return _fir_mirror_symmetric(in, out, N, h, Nh, instride, outstride);
}
#endif



// _separable_2Dconvolve_mirror instantiations

int S_separable_2Dconvolve_mirror(float *in, float *out,
                              int M, int N, float *hr, float *hc,
                              int Nhr, int Nhc,
                              npy_intp *instrides, npy_intp *outstrides)
{
    return _separable_2Dconvolve_mirror(in, out, M, N, hr, hc, Nhr, Nhc, instrides, outstrides);
}

int D_separable_2Dconvolve_mirror(double *in, double *out,
                              int M, int N, double *hr, double *hc,
                              int Nhr, int Nhc,
                              npy_intp *instrides, npy_intp *outstrides)
{
    return _separable_2Dconvolve_mirror(in, out, M, N, hr, hc, Nhr, Nhc, instrides, outstrides);
}


#ifdef __GNUC__
int C_separable_2Dconvolve_mirror(__complex__ float *in, __complex__ float *out,
                              int M, int N, __complex__ float *hr, __complex__ float *hc,
                              int Nhr, int Nhc,
                              npy_intp *instrides, npy_intp *outstrides)
{
    return _separable_2Dconvolve_mirror(in, out, M, N, hr, hc, Nhr, Nhc, instrides, outstrides);
}

int Z_separable_2Dconvolve_mirror(__complex__ double *in, __complex__ double *out,
                              int M, int N, __complex__ double *hr, __complex__ double *hc,
                              int Nhr, int Nhc,
                              npy_intp *instrides, npy_intp *outstrides)
{
    return _separable_2Dconvolve_mirror(in, out, M, N, hr, hc, Nhr, Nhc, instrides, outstrides);
}
#endif



/**
Approximate the steady state of a two-pole filter in polar form for a
step input.
**/
float S_hc(int k, float cs, double r, double omega)
{
    return _hc(k, cs, r, omega);
}


double D_hc(int k, double cs, double r, double omega)
{
    return _hc(k, cs, r, omega);
}


/**
Approximate the steady state of a two-pole filer in polar form ran in backwards
for a step input.
**/

float S_hs(int k, float cs, double rsq, double omega)
{
    return _hs(k, cs, rsq, omega);
}

double D_hs(int k, double cs, double rsq, double omega)
{
    return _hs(k, cs, rsq, omega);
}


