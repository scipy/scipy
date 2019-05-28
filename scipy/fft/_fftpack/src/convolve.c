/*
  Generic functions for computing 1D convolutions of periodic sequences.

  Supported FFT libraries:
    DJBFFT  - optional, used for power-of-two length arrays
    FFTW    - optional
    FFTPACK - used if any of the above libraries is not available

  Author: Pearu Peterson, September 2002
 */

#include "fftpack.h"

/**************** FFTPACK ZFFT **********************/
extern void F_FUNC(dfftf, DFFTF) (int *, double *, double *);
extern void F_FUNC(dfftb, DFFTB) (int *, double *, double *);
extern void F_FUNC(dffti, DFFTI) (int *, double *);
GEN_CACHE(dfftpack, (int n)
          , double *wsave;, (caches_dfftpack[i].n == n)
          , caches_dfftpack[id].wsave =
          (double *) malloc(sizeof(double) * (2 * n + 15));
          F_FUNC(dffti, DFFTI) (&n, caches_dfftpack[id].wsave);,
          free(caches_dfftpack[id].wsave);, 20)

extern void destroy_convolve_cache(void)
{
    destroy_dfftpack_cache();
}

/**************** convolve **********************/
extern void
convolve(int n, double *inout, double *omega, int swap_real_imag)
{
    int i;
    double *wsave = NULL;

    i = get_cache_id_dfftpack(n);
    wsave = caches_dfftpack[i].wsave;
    F_FUNC(dfftf, DFFTF) (&n, inout, wsave);
    if (swap_real_imag) {
        double c;
        int n1 = n - 1;
        inout[0] *= omega[0];
        if (!(n % 2))
            inout[n - 1] *= omega[n - 1];
        for (i = 1; i < n1; i += 2) {
            c = inout[i] * omega[i];
            inout[i] = inout[i + 1] * omega[i + 1];
            inout[i + 1] = c;
        }
    } else
        for (i = 0; i < n; ++i)
            inout[i] *= omega[i];
    F_FUNC(dfftb, DFFTB) (&n, inout, wsave);
}

/**************** convolve **********************/
extern void
convolve_z(int n, double *inout, double *omega_real, double *omega_imag)
{
    int i;
    double *wsave = NULL;
    i = get_cache_id_dfftpack(n);
    wsave = caches_dfftpack[i].wsave;
    F_FUNC(dfftf, DFFTF) (&n, inout, wsave);
    {
        double c;
        int n1 = n - 1;
        inout[0] *= (omega_real[0] + omega_imag[0]);
        if (!(n % 2))
            inout[n - 1] *= (omega_real[n - 1] + omega_imag[n - 1]);
        for (i = 1; i < n1; i += 2) {
            c = inout[i] * omega_imag[i];
            inout[i] *= omega_real[i];
            inout[i] += inout[i + 1] * omega_imag[i + 1];
            inout[i + 1] *= omega_real[i + 1];
            inout[i + 1] += c;
        }
    }
    F_FUNC(dfftb, DFFTB) (&n, inout, wsave);
}

extern void
init_convolution_kernel(int n, double *omega, int d,
                        double (*kernel_func) (int), int zero_nyquist)
{
    /*
       omega[k] = pow(sqrt(-1),d) * kernel_func(k)
       omega[0] = kernel_func(0)
       conjugate(omega[-k]) == omega[k]
     */
    int j, k, l = (n % 2 ? n : n - 1);
    omega[0] = (*kernel_func) (0) / n;
    switch (d % 4) {
        case 0:
            for (k = j = 1; j < l; j += 2, ++k)
                omega[j] = omega[j + 1] = (*kernel_func) (k) / n;
            if (!(n % 2))
                omega[n - 1] =
                    (zero_nyquist ? 0.0 : (*kernel_func) (k) / n);
            break;
        case 1:;
        case -3:
            for (k = j = 1; j < l; j += 2, ++k) {
                omega[j] = (*kernel_func) (k) / n;
                omega[j + 1] = -omega[j];
            }
            if (!(n % 2))
                omega[n - 1] =
                    (zero_nyquist ? 0.0 : (*kernel_func) (k) / n);
            break;
        case 2:;
        case -2:
            for (k = j = 1; j < l; j += 2, ++k)
                omega[j] = omega[j + 1] = -(*kernel_func) (k) / n;
            if (!(n % 2))
                omega[n - 1] =
                    (zero_nyquist ? 0.0 : -(*kernel_func) (k) / n);
            break;
        case 3:;
        case -1:
            for (k = j = 1; j < l; j += 2, ++k) {
                omega[j] = -(*kernel_func) (k) / n;
                omega[j + 1] = -omega[j];
            }
            if (!(n % 2))
                omega[n - 1] =
                    (zero_nyquist ? 0.0 : -(*kernel_func) (k) / n);
            break;
    }
}
