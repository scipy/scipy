#pragma once
#include "Python.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define PYERR(message)                                                                                                 \
    do {                                                                                                               \
        PyErr_SetString(PyExc_ValueError, message);                                                                    \
        goto fail;                                                                                                     \
    } while (0)


/**
* Compute the starting initial conditions for the system
*
*    1 / (1 - z1 * z^-1) against signal x.
*
* The initial conditions are defined as x[0] + z1 * \sum{k = 0}^{n - 1} x[k] * z1^k
* this sum will be aggregated until its precision is below a certain threshold.
*
* Arguments
* ----------------
*  z1: double or complex
*      Exponential decay parameter in the transfer function
*  x: float* or double* or complex*
*     2D strided pointer signal of size (M, N). When M > 1, multiple signals
*     will be processed independently.
*  yp0: float* or double* or complex*
*       Output state condition pointer of size (M, 1)
*  M: int
*      Number of signals to compute initial conditions for.
*  N: int
*      Length of the signals.
*  precision: double* or float*
*      Precision up to which the initial conditions will be computed.
**/
template <typename T, typename C>
int _sym_iir1_initial(C z1, C *x, C *yp0, int M, int N, T precision) {
    // XXX: remove templating on C,T : C === T or std::complex<T>
    C powz1, diff;
    T err;
    int k;

    if (std::abs(z1) >= 1.0)
        return -2; /* z1 not less than 1 */

   /* Fix starting value assuming mirror-symmetric boundary conditions. */
    for(int i = 0; i < M; i++) {
        yp0[i] = x[N * i];
    }

    powz1 = 1.0;
    k = 0;
    precision *= precision;
    do {
        powz1 *= z1;
        for(int i = 0; i < M; i++) {
            yp0[i] += powz1 * x[N * i + k];
        }
        diff = powz1;
        err = std::abs(diff);
        err *= err;
        k++;
    } while((err > precision) && (k < N));

    if (k >= N){
        /* sum did not converge */
        return -3;
    }

    return 0;
}


/**
Approximate the steady state of a two-pole filer in polar form ran in backwards
for a step input.
**/
template<typename T>
T _hs(int k, T cs, double rsq, double omega) {
    T cssq;
    T c0;
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


/**
Approximate the steady state of a two-pole filter in polar form for a
step input.
**/
template<typename T>
T _hc(int k, T cs, double r, double omega) {
    if (k < 0)
        return 0.0;
    if (omega == 0.0)
        return cs * pow(r, (double )k) * (k+1);
    else if (omega == M_PI)
        return cs * pow(r, (double )k) * (k+1) * (1 - 2*(k % 2));
    return cs * pow(r, (double) k) * sin(omega * (k+1)) / sin(omega);
}



/**
Compute the starting initial conditions for the system (ran in backwards)
cs / (1 - a2 * z - a3 * z^2) against signal x, where:

a2 = 2 * r * cos(omega)
a3 = - r ** 2
cs = 1 - 2 * r * cos(omega) + r ** 2

Arguments
---------
x: double* or float*
    2D strided pointer signal of size (M, N). When M > 1, multiple signals
    will be processed independently.
yp: double* or float*
    Output state condition pointer of size (M, 2), yp[:, 0] will contain the
    initial conditions y[n + 1], whereas yp[:, 1] will contain the initial
    condition y[n + 2].
M: int
    Number of signals to compute initial conditions for.
N: int
    Length of the signals.
precision: double* or float*
    Precision up to which the initial conditions will be computed.
**/
template<typename T>
int _sym_iir2_initial_bwd(double r, double omega, T *x, T *yp, int M, int N, T precision) {
    double rsq = r * r;
    T cs = 1 - 2 * r * cos(omega) + rsq;

    // Fix starting values assuming mirror-symmetric boundary conditions.
    int k = 0;

    T err;
    T diff;

    do {
        diff = (_hs(k, cs, rsq, omega) + _hs(k+1, cs, rsq, omega));
        for(int i = 0; i < M; i++) {
            // Compute initial condition y[n + 1]
            yp[2 * i] += diff * x[N * i + N - 1 - k];
        }
        err = diff * diff;
        k++;
    } while((err > precision) && (k < N));

    if (k >= N) {
        return -3;
    } // sum did not converge

    k = 0;
    do {
        diff = (_hs(k-1, cs, rsq, omega) + _hs(k+2, cs, rsq, omega));
        for(int i = 0; i < M; i++) {
            // Compute initial condition y[n + 2]
            yp[2 * i + 1] += diff * x[N * i + N - 1 - k];
        }
        err = diff * diff;
        k++;
    } while((err > precision) && (k < N));

    if (k >= N) {
        return -3;
    } // sum did not converge

    return 0;
}



/**
Compute the starting initial conditions for the system
cs / (1 - a2 * z^-1 - a3 * z^-2) against signals x, where:

a2 = 2 * r * cos(omega)
a3 = - r ** 2
cs = 1 - 2 * r * cos(omega) + r ** 2

Arguments
---------
x: double* or float*
    2D strided pointer signal of size (M, N). When M > 1, multiple signals
    will be processed independently.
yp: double* or float*
    Strided output state condition pointer of size (M, 2).
    yp[:, 0] will contain the initial conditions y[n - 1],
    whereas yp[:, 1] will contain the initial conditions y[n - 2].
M: int
    Number of signals to compute initial conditions for.
N: int
    Length of the signals.
precision: double* or float*
    Precision up to which the initial conditions will be computed.
**/
template<typename T>
int _sym_iir2_initial_fwd(double r, double omega, T *x, T *yp, int M, int N, T precision) {
    /* Fix starting values assuming mirror-symmetric boundary conditions. */
    T cs = 1 - 2 * r * cos(omega) + r * r;

    for(int i = 0; i < M; i++) {
        // Compute starting condition y[n - 1]
        yp[2 * i] = _hc(0, cs, r, omega) * x[N * i];
    }

    int k = 0;
    precision *= precision;

    T err;
    T diff;

    do {
        diff = _hc(k+1, cs, r, omega);
        for(int i = 0; i < M; i++) {
            // Keep computing starting condition y[n - 1]
            yp[2 * i] += diff * x[N * i + k];
        }
        err = diff * diff;
        k++;
    } while((err > precision) && (k < N));

    if (k >= N) {
        return -3;
    } /* sum did not converge */

    for(int i = 0; i < M; i++) {
        // Compute starting condition y[n - 2]
        yp[2 * i + 1] = _hc(0, cs, r, omega) * x[N * i + 1];
        yp[2 * i + 1] += _hc(1, cs, r, omega) * x[N * i];
    }

    k = 0;
    do {
        diff = _hc(k+2, cs, r, omega);
        for(int i = 0; i < M; i++) {
            // Keep computing starting condition y[n - 2]
            yp[2 * i + 1] += diff * x[N * i + k];
        }
        err = diff * diff;
        k++;
    } while((err > precision) && (k < N));

    if (k >= N) {
        return -3;
    } /* sum did not converge */
    return 0;
}



template <typename T>
void _fir_mirror_symmetric(T *in, T *out, int N, T *h, int Nh, int instride, int outstride) {
    int n, k;
    int Nhdiv2 = Nh >> 1;
    T *outptr;
    T *inptr;
    T *hptr;

    /* first part boundary conditions */
    outptr = out;
    for (n=0; n < Nhdiv2; n++) {
        *outptr = 0.0;
        hptr = h;
        inptr = in + (n + Nhdiv2)*instride;
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
        inptr = in + (n + Nhdiv2)*instride;
        for (k=-Nhdiv2; k <= Nhdiv2; k++) {
            *outptr += *hptr++ * *inptr;
            inptr -= instride;
        }
        outptr += outstride;
    }

    /* end boundary conditions */
    outptr = out + (N - Nhdiv2)*outstride;
    for (n=N-Nhdiv2; n < N; n++) {
        *outptr = 0.0;
        hptr = h;
        inptr = in + (2*N - 1 - n - Nhdiv2)*instride;
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



template<typename T>
int _separable_2Dconvolve_mirror(T *in, T *out, int M, int N, T *hr, T *hc, int Nhr, int Nhc, Py_ssize_t *instrides,
                                 Py_ssize_t *outstrides) {
    int m, n;
    T *tmpmem;
    T *inptr = NULL;
    T *outptr = NULL;

    tmpmem = (T *)malloc(M*N*sizeof(T));
    if (tmpmem == NULL) {
        return -1;
    }

    if (Nhr > 0) {
        /* filter across rows */
        inptr = in;
        outptr = tmpmem;
        for (m = 0; m < M; m++) {
            _fir_mirror_symmetric (inptr, outptr, N, hr, Nhr, instrides[1], 1);
            inptr += instrides[0];
            outptr += N;
        }
    } else
        memmove(tmpmem, in, M*N*sizeof(T));

    if (Nhc > 0) {
        /* filter down columns */
        inptr = tmpmem;
        outptr = out;
        for (n = 0; n < N; n++) {
            _fir_mirror_symmetric (inptr, outptr, M, hc, Nhc, N, outstrides[0]);
            outptr += outstrides[1];
            inptr += 1;
        }
    } else
        memmove(out, tmpmem, M*N*sizeof(T));

    free(tmpmem);
    return 0;
}
