/*
  Interface to FFTLog.
  Author: Dieter Werthmüller, January 2017
 */

#include "fftpack.h"

extern void F_FUNC(dffti, DFFTI) (int *, double *);
extern void F_FUNC(fftl, FFTL) (int *, double *, double *, int *, double *,
                                double *);
extern void F_FUNC(fhti, FHTI) (int *, double *, double *, double *, double *,
                                double *);

/* Cache fftl_w is the same as cache drfft in drfft.c, and could potentially be
 * combined. */
GEN_CACHE(fftl_w, (int n)
      , double *wsave;
      , (caches_fftl_w[i].n == n)
      , caches_fftl_w[id].wsave =
      (double *) malloc(sizeof(double) * (2 * n + 15));
      F_FUNC(dffti, DFFTI) (&n, caches_fftl_w[id].wsave);
      , free(caches_fftl_w[id].wsave);
      , 10)

GEN_CACHE(fftl_x
      , (int n, double mu, double q, double dlnr, double kr, int size)
      , double *xsave; double mu; double q; double dlnr; double kr;
      , ((caches_fftl_x[i].n == n) && (caches_fftl_x[i].mu == mu) &&
         (caches_fftl_x[i].q == q) && (caches_fftl_x[i].dlnr == dlnr) &&
         (caches_fftl_x[i].kr == kr))
      , caches_fftl_x[id].xsave = (double *)
        malloc(sizeof(double) * size);
        F_FUNC(fhti, FHTI) (&n, &mu, &q, &dlnr, &kr, caches_fftl_x[id].xsave);
      , free(caches_fftl_x[id].xsave);
      , 10)


void drfftl(double *inout, int n, double mu, double q, double dlnr, double kr,
            double rk, int direction)
{
    int size;
    double *ptr = inout;
    double *xsave = NULL;
    double *wsave = NULL;

    if ( q != 0) {
        size = 3 * (n / 2) + 4;
    }
    else {
        size = 2 * (n / 2) + 3;
    }

    wsave = caches_fftl_w[get_cache_id_fftl_w(n)].wsave;
    xsave = caches_fftl_x[get_cache_id_fftl_x(n, mu, q, dlnr, kr, size)].xsave;

    F_FUNC(fftl, FFTL)(&n, ptr, &rk, &direction, wsave, xsave);
}
