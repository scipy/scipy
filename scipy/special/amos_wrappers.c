#include "amos_wrappers.h"
#define CADDR(z) ((double *) (&(z))), (&(((double *) (&(z)))[1]))

#ifndef CMPLX
#define CMPLX(x, y) ((double complex)((double)(x) + I * (double)(y)))
#endif /* CMPLX */

int ierr_to_sferr(int nz, int ierr) {
  /* Return sf_error equivalents for ierr values */

  if (nz != 0) return SF_ERROR_UNDERFLOW;

  switch (ierr) {
  case 1:
    return SF_ERROR_DOMAIN;
  case 2:
    return SF_ERROR_OVERFLOW;
  case 3:
    return SF_ERROR_LOSS;
  case 4:
    return SF_ERROR_NO_RESULT;
  case 5:   /* Algorithm termination condition not met */
    return SF_ERROR_NO_RESULT;
  }
  return -1;
}


void set_nan_if_no_computation_done(npy_cdouble *v, int ierr) {
  if (v != NULL && (ierr == 1 || ierr == 2 || ierr == 4 || ierr == 5)) {
    npymath_csetreal(v, NAN);
    npymath_csetimag(v, NAN);
  }
}

double sin_pi(double x)
{
    if (floor(x) == x && fabs(x) < 1e14) {
        /* Return 0 when at exact zero, as long as the floating point number is
         * small enough to distinguish integer points from other points.
         */
        return 0;
    }
    return sin(M_PI * x);
}

static double cos_pi(double x)
{
    double x05 = x + 0.5;
    if (floor(x05) == x05 && fabs(x) < 1e14) {
        /* Return 0 when at exact zero, as long as the floating point number is
         * small enough to distinguish integer points from other points.
         */
        return 0;
    }
    return cos(M_PI * x);
}

static npy_cdouble
rotate(npy_cdouble z, double v)
{
    npy_cdouble w;
    double c = cos_pi(v);
    double s = sin_pi(v);
    npymath_csetreal(&w, npymath_creal(z)*c - npymath_cimag(z)*s);
    npymath_csetimag(&w, npymath_creal(z)*s + npymath_cimag(z)*c);
    return w;
}

static npy_cdouble
rotate_jy(npy_cdouble j, npy_cdouble y, double v)
{
    npy_cdouble w;
    double c = cos_pi(v);
    double s = sin_pi(v);
    npymath_csetreal(&w, npymath_creal(j) * c - npymath_creal(y) * s);
    npymath_csetimag(&w, npymath_cimag(j) * c - npymath_cimag(y) * s);
    return w;
}

static int
reflect_jy(npy_cdouble *jy, double v)
{
    /* NB: Y_v may be huge near negative integers -- so handle exact
     *     integers carefully
     */
    int i;
    if (v != floor(v))
        return 0;

    i = v - 16384.0 * floor(v / 16384.0);
    if (i & 1) {
        npymath_csetreal(jy, -npymath_creal(*jy));
        npymath_csetimag(jy, -npymath_cimag(*jy));
    }
    return 1;
}

static int
reflect_i(npy_cdouble *ik, double v)
{
    if (v != floor(v))
        return 0;
    return 1; /* I is symmetric for integer v */
}

static npy_cdouble
rotate_i(npy_cdouble i, npy_cdouble k, double v)
{
    npy_cdouble w;
    double s = sin(v * M_PI)*(2.0/M_PI);
    npymath_csetreal(&w, npymath_creal(i) + s*npymath_creal(k));
    npymath_csetimag(&w, npymath_cimag(i) + s*npymath_cimag(k));
    return w;
}

int cephes_airy(double x, double *ai, double *aip, double *bi, double *bip);

int airy_wrap(double x, double *ai, double *aip, double *bi, double *bip)
{
    npy_cdouble z, zai, zaip, zbi, zbip;

    /* For small arguments, use Cephes as it's slightly faster.
     * For large arguments, use AMOS as it's more accurate.
     */
    if (x < -10 || x > 10) {
        npymath_csetreal(&z, x);
        npymath_csetimag(&z, 0);
        cairy_wrap(z, &zai, &zaip, &zbi, &zbip);
        *ai  = npymath_creal(zai);
        *aip = npymath_creal(zaip);
        *bi  = npymath_creal(zbi);
        *bip = npymath_creal(zbip);
    }
    else {
        cephes_airy(x, ai, aip, bi, bip);
    }
    return 0;
}

int cairy_wrap(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
  int id = 0;
  int ierr = 0;
  int kode = 1;
  int nz;
  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex res;
  
  npymath_csetreal(ai, NAN);
  npymath_csetimag(ai, NAN);
  npymath_csetreal(bi, NAN);
  npymath_csetimag(bi, NAN);
  npymath_csetreal(aip, NAN);
  npymath_csetimag(aip, NAN);
  npymath_csetreal(bip, NAN);
  npymath_csetimag(bip, NAN);

  res = amos_airy(z99, id, kode, &nz, &ierr);
  npymath_csetreal(ai, creal(res));
  npymath_csetimag(ai, cimag(res));
  DO_SFERR("airy:", ai);

  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  npymath_csetreal(bi, creal(res));
  npymath_csetimag(bi, cimag(res));
  DO_SFERR("airy:", bi);

  id = 1;
  res = amos_airy(z99, id, kode, &nz, &ierr);
  npymath_csetreal(aip, creal(res));
  npymath_csetimag(aip, cimag(res));
  DO_SFERR("airy:", aip);

  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  npymath_csetreal(bip, creal(res));
  npymath_csetimag(bip, cimag(res));
  DO_SFERR("airy:", bip);
  return 0;
}

int cairy_wrap_e(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
  int id = 0;
  int kode = 2;        /* Exponential scaling */
  int nz, ierr;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex res;

  npymath_csetreal(ai, NAN);
  npymath_csetimag(ai, NAN);
  npymath_csetreal(bi, NAN);
  npymath_csetimag(bi, NAN);
  npymath_csetreal(aip, NAN);
  npymath_csetimag(aip, NAN);
  npymath_csetreal(bip, NAN);
  npymath_csetimag(bip, NAN);

  res = amos_airy(z99, id, kode, &nz, &ierr);
  npymath_csetreal(ai, creal(res));
  npymath_csetimag(ai, cimag(res));
  DO_SFERR("airye:", ai);

  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  npymath_csetreal(bi, creal(res));
  npymath_csetimag(bi, cimag(res));
  DO_SFERR("airye:", bi);

  id = 1;
  res = amos_airy(z99, id, kode, &nz, &ierr);
  npymath_csetreal(aip, creal(res));
  npymath_csetimag(aip, cimag(res));
  DO_SFERR("airye:", aip);

  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  npymath_csetreal(bip, creal(res));
  npymath_csetimag(bip, cimag(res));
  DO_SFERR("airye:", bip);
  return 0;
}

int cairy_wrap_e_real(double z, double *ai, double *aip, double *bi, double *bip) {
  int id = 0;
  int kode = 2;        /* Exponential scaling */
  int nz, ierr;
  npy_cdouble cai, caip, cbi, cbip;
  
  double complex z99 = z;
  double complex res;

  npymath_csetreal(&cai, NAN);
  npymath_csetimag(&cai, NAN);
  npymath_csetreal(&cbi, NAN);
  npymath_csetimag(&cbi, NAN);
  npymath_csetreal(&caip, NAN);
  npymath_csetimag(&caip, NAN);
  npymath_csetreal(&cbip, NAN);
  npymath_csetimag(&cbip, NAN);


  if (z < 0) {
      *ai = NAN;
  } else {
      res = amos_airy(z99, id, kode, &nz, &ierr);
      npymath_csetreal(&cai, creal(res));
      npymath_csetimag(&cai, cimag(res));
      DO_SFERR("airye:", &cai);
      *ai = npymath_creal(cai);
  }
  
  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  npymath_csetreal(&cbi, creal(res));
  npymath_csetimag(&cbi, cimag(res));
  DO_SFERR("airye:", &cbi);
  *bi = npymath_creal(cbi);

  id = 1;
  if (z < 0) {
      *aip = NAN;
  } else {
      res = amos_airy(z99, id, kode, &nz, &ierr);
      npymath_csetreal(&caip, creal(res));
      npymath_csetimag(&caip, cimag(res));
      DO_SFERR("airye:", &caip);
      *aip = npymath_creal(caip);
  }

  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  npymath_csetreal(&cbip, creal(res));
  npymath_csetimag(&cbip, cimag(res));
  DO_SFERR("airye:", &cbip);
  *bip = npymath_creal(cbip);
  return 0;
}

npy_cdouble cbesi_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int sign = 1;
  int nz, ierr;
  npy_cdouble cy, cy_k;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy99[1] = { NAN };
  double complex cy_k99[1] = { NAN };

  npymath_csetreal(&cy, NAN);
  npymath_csetimag(&cy, NAN);
  npymath_csetreal(&cy_k, NAN);
  npymath_csetimag(&cy_k, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besi(z99, v, kode, n, cy99, &ierr);
  npymath_csetreal(&cy, creal(cy99[0]));
  npymath_csetimag(&cy, cimag(cy99[0]));
  DO_SFERR("iv:", &cy);
  if (ierr == 2) {
    /* overflow */
    if (npymath_cimag(z) == 0 && (npymath_creal(z) >= 0 || v == floor(v))) {
        if (npymath_creal(z) < 0 && v/2 != floor(v/2))
            npymath_csetreal(&cy, -INFINITY);
        else
            npymath_csetreal(&cy, INFINITY);
        npymath_csetimag(&cy, 0);
    } else {
        cy = cbesi_wrap_e(v*sign, z);
        npymath_csetreal(&cy, npymath_creal(cy) * INFINITY);
        npymath_csetimag(&cy, npymath_cimag(cy) * INFINITY);
    }
  }

  if (sign == -1) {
    if (!reflect_i(&cy, v)) {
      nz = amos_besk(z99, v, kode, n, cy_k99, &ierr);
      npymath_csetreal(&cy_k, creal(cy_k99[0]));
      npymath_csetimag(&cy_k, cimag(cy_k99[0]));
      DO_SFERR("iv(kv):", &cy_k);
      cy = rotate_i(cy, cy_k, v);
    }
  }

  return cy;
}

npy_cdouble cbesi_wrap_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int sign = 1;
  int nz, ierr;
  npy_cdouble cy, cy_k;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy99[1] = { NAN };
  double complex cy_k99[1] = { NAN };

  npymath_csetreal(&cy, NAN);
  npymath_csetimag(&cy, NAN);
  npymath_csetreal(&cy_k, NAN);
  npymath_csetimag(&cy_k, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besi(z99, v, kode, n, cy99, &ierr);
  npymath_csetreal(&cy, creal(cy99[0]));
  npymath_csetimag(&cy, cimag(cy99[0]));
  DO_SFERR("ive:", &cy);

  if (sign == -1) {
    if (!reflect_i(&cy, v)) {
      nz = amos_besk(z99, v, kode, n, cy_k99, &ierr);
      npymath_csetreal(&cy_k, creal(cy_k99[0]));
      npymath_csetimag(&cy_k, cimag(cy_k99[0]));
      DO_SFERR("ive(kv):", &cy_k);
      /* adjust scaling to match zbesi */
      cy_k = rotate(cy_k, -npymath_cimag(z)/M_PI);
      if (npymath_creal(z) > 0) {
          npymath_csetreal(&cy_k, npymath_creal(cy_k) * exp(-2*npymath_creal(z)));
          npymath_csetimag(&cy_k, npymath_cimag(cy_k) * exp(-2*npymath_creal(z)));
      }
      /* v -> -v */
      cy = rotate_i(cy, cy_k, v);
    }
  }

  return cy;
}

double cbesi_wrap_e_real(double v, double z) {
  npy_cdouble cy, w;
  if (v != floor(v) && z < 0) {
    return NAN;
  } else {
    npymath_csetreal(&w, z);
    npymath_csetimag(&w, 0);
    cy = cbesi_wrap_e(v, w);
    return npymath_creal(cy);
  }
}

npy_cdouble cbesj_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_j, cy_y;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy_j99[1] = { NAN };
  double complex cy_y99[1] = { NAN };

  npymath_csetreal(&cy_j, NAN);
  npymath_csetimag(&cy_j, NAN);
  npymath_csetreal(&cy_y, NAN);
  npymath_csetimag(&cy_y, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy_j;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besj(z99, v, kode, n, cy_j99, &ierr);
  npymath_csetreal(&cy_j, creal(cy_j99[0]));
  npymath_csetimag(&cy_j, cimag(cy_j99[0]));
  DO_SFERR("jv:", &cy_j);
  if (ierr == 2) {
    /* overflow */
    cy_j = cbesj_wrap_e(v, z);
    npymath_csetreal(&cy_j, npymath_creal(cy_j) * INFINITY);
    npymath_csetimag(&cy_j, npymath_cimag(cy_j) * INFINITY);
  }

  if (sign == -1) {
    if (!reflect_jy(&cy_j, v)) {
      nz = amos_besy(z99, v, kode, n, cy_y99, &ierr);
      npymath_csetreal(&cy_y, creal(cy_y99[0]));
      npymath_csetimag(&cy_y, cimag(cy_y99[0]));
      DO_SFERR("jv(yv):", &cy_y);
      cy_j = rotate_jy(cy_j, cy_y, v);
    }
  }
  return cy_j;
}

double cephes_jv(double v, double x);

double cbesj_wrap_real(double v, double x)
{
    npy_cdouble z, r;

    if (x < 0 && v != (int)v) {
        sf_error("yv", SF_ERROR_DOMAIN, NULL);
        return NAN;
    }

    npymath_csetreal(&z, x);
    npymath_csetimag(&z, 0);
    r = cbesj_wrap(v, z);
    if (npymath_creal(r) != npymath_creal(r)) {
        /* AMOS returned NaN, possibly due to overflow */
        return cephes_jv(v, x);
    }
    return npymath_creal(r);
}

npy_cdouble cbesj_wrap_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_j, cy_y;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy_j99[1] = { NAN };
  double complex cy_y99[1] = { NAN };

  npymath_csetreal(&cy_j, NAN);
  npymath_csetimag(&cy_j, NAN);
  npymath_csetreal(&cy_y, NAN);
  npymath_csetimag(&cy_y, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy_j;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besj(z99, v, kode, n, cy_j99, &ierr);
  npymath_csetreal(&cy_j, creal(cy_j99[0]));
  npymath_csetimag(&cy_j, cimag(cy_j99[0]));
  DO_SFERR("jve:", &cy_j);
  if (sign == -1) {
    if (!reflect_jy(&cy_j, v)) {
      nz = amos_besy(z99, v, kode, n, cy_y99, &ierr);
      npymath_csetreal(&cy_y, creal(cy_y99[0]));
      npymath_csetimag(&cy_y, cimag(cy_y99[0]));
      DO_SFERR("jve(yve):", &cy_y);
      cy_j = rotate_jy(cy_j, cy_y, v);
    }
  }
  return cy_j;
}

double cbesj_wrap_e_real(double v, double z) {
  npy_cdouble cy, w;
  if (v != floor(v) && z < 0) {
    return NAN;
  } else {
    npymath_csetreal(&w, z);
    npymath_csetimag(&w, 0);
    cy = cbesj_wrap_e(v, w);
    return npymath_creal(cy);
  }
}

npy_cdouble cbesy_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_y, cy_j;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy_j99[1] = { NAN };
  double complex cy_y99[1] = { NAN };

  npymath_csetreal(&cy_j, NAN);
  npymath_csetimag(&cy_j, NAN);
  npymath_csetreal(&cy_y, NAN);
  npymath_csetimag(&cy_y, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy_y;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }

  if (npymath_creal(z) == 0 && npymath_cimag(z) == 0) {
      /* overflow */
      npymath_csetreal(&cy_y, -INFINITY);
      npymath_csetimag(&cy_y, 0);
      sf_error("yv", SF_ERROR_OVERFLOW, NULL);
  }
  else {
      nz = amos_besy(z99, v, kode, n, cy_y99, &ierr);
      npymath_csetreal(&cy_y, creal(cy_y99[0]));
      npymath_csetimag(&cy_y, cimag(cy_y99[0]));
      DO_SFERR("yv:", &cy_y);
      if (ierr == 2) {
          if (npymath_creal(z) >= 0 && npymath_cimag(z) == 0) {
              /* overflow */
              npymath_csetreal(&cy_y, -INFINITY);
              npymath_csetimag(&cy_y, 0);
          }
      }
  }

  if (sign == -1) {
    if (!reflect_jy(&cy_y, v)) {
      nz = amos_besj(z99, v, kode, n, cy_j99, &ierr);
      npymath_csetreal(&cy_j, creal(cy_j99[0]));
      npymath_csetimag(&cy_j, cimag(cy_j99[0]));
      // F_FUNC(zbesj,ZBESJ)(CADDR(z), &v,  &kode, &n, CADDR(cy_j), &nz, &ierr);
      DO_SFERR("yv(jv):", &cy_j);
      cy_y = rotate_jy(cy_y, cy_j, -v);
    }
  }
  return cy_y;
}

double cephes_yv(double v, double x);

double cbesy_wrap_real(double v, double x)
{
    npy_cdouble z, r;

    if (x < 0.0) {
        sf_error("yv", SF_ERROR_DOMAIN, NULL);
        return NAN;
    }

    npymath_csetreal(&z, x);
    npymath_csetimag(&z, 0);
    r = cbesy_wrap(v, z);
    if (npymath_creal(r) != npymath_creal(r)) {
        /* AMOS returned NaN, possibly due to overflow */
        return cephes_yv(v, x);
    }
    return npymath_creal(r);
}

npy_cdouble cbesy_wrap_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_y, cy_j;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy_j99[1] = { NAN };
  double complex cy_y99[1] = { NAN };

  npymath_csetreal(&cy_j, NAN);
  npymath_csetimag(&cy_j, NAN);
  npymath_csetreal(&cy_y, NAN);
  npymath_csetimag(&cy_y, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy_y;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besy(z99, v, kode, n, cy_y99, &ierr);
  npymath_csetreal(&cy_y, creal(cy_y99[0]));
  npymath_csetimag(&cy_y, cimag(cy_y99[0]));
  DO_SFERR("yve:", &cy_y);
  if (ierr == 2) {
    if (npymath_creal(z) >= 0 && npymath_cimag(z) == 0) {
      /* overflow */
      npymath_csetreal(&cy_y, INFINITY);
      npymath_csetimag(&cy_y, 0);
    }
  }

  if (sign == -1) {
    if (!reflect_jy(&cy_y, v)) {
      nz = amos_besj(z99, v, kode, n, cy_j99, &ierr);
      npymath_csetreal(&cy_j, creal(cy_j99[0]));
      npymath_csetimag(&cy_j, cimag(cy_j99[0]));
      DO_SFERR("yv(jv):", &cy_j);
      cy_y = rotate_jy(cy_y, cy_j, -v);
    }
  }
  return cy_y;
}

double cbesy_wrap_e_real(double v, double z) {
  npy_cdouble cy, w;
  if (z < 0) {
    return NAN;
  } else {
    npymath_csetreal(&w, z);
    npymath_csetimag(&w, 0);
    cy = cbesy_wrap_e(v, w);
    return npymath_creal(cy);
  }
}

npy_cdouble cbesk_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  npy_cdouble cy;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy99[1] = { NAN };

  npymath_csetreal(&cy, NAN);
  npymath_csetimag(&cy, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    /* K_v == K_{-v} even for non-integer v */
    v = -v;
  }
  nz = amos_besk(z99, v, kode, n, cy99, &ierr);
  npymath_csetreal(&cy, creal(cy99[0]));
  npymath_csetimag(&cy, cimag(cy99[0]));
  DO_SFERR("kv:", &cy);
  if (ierr == 2) {
    if (npymath_creal(z) >= 0 && npymath_cimag(z) == 0) {
      /* overflow */
      npymath_csetreal(&cy, INFINITY);
      npymath_csetimag(&cy, 0);
    }
  }

  return cy;
}

npy_cdouble cbesk_wrap_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  npy_cdouble cy;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy99[1] = { NAN };

  npymath_csetreal(&cy, NAN);
  npymath_csetimag(&cy, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    /* K_v == K_{-v} even for non-integer v */
    v = -v;
  }
  nz = amos_besk(z99, v, kode, n, cy99, &ierr);
  npymath_csetreal(&cy, creal(cy99[0]));
  npymath_csetimag(&cy, cimag(cy99[0]));
  DO_SFERR("kve:", &cy);
  if (ierr == 2) {
    if (npymath_creal(z) >= 0 && npymath_cimag(z) == 0) {
      /* overflow */
      npymath_csetreal(&cy, INFINITY);
      npymath_csetimag(&cy, 0);
    }
  }

  return cy;
}

double cbesk_wrap_real( double v, double z) {
  npy_cdouble cy, w;
  if (z < 0) {
    return NAN;
  }
  else if (z == 0) {
    return INFINITY;
  }
  else if (z > 710 * (1 + fabs(v))) {
      /* Underflow. See uniform expansion https://dlmf.nist.gov/10.41
       * This condition is not a strict bound (it can underflow earlier),
       * rather, we are here working around a restriction in AMOS.
       */
      return 0;
  }
  else {
    npymath_csetreal(&w, z);
    npymath_csetimag(&w, 0);
    cy = cbesk_wrap(v, w);
    return npymath_creal(cy);
  }
}

double cbesk_wrap_real_int(int n, double z)
{
    return cbesk_wrap_real(n, z);
}

double cbesk_wrap_e_real( double v, double z) {
  npy_cdouble cy, w;
  if (z < 0) {
    return NAN;
  }
  else if (z == 0) {
    return INFINITY;
  }
  else {
    npymath_csetreal(&w, z);
    npymath_csetimag(&w, 0);
    cy = cbesk_wrap_e(v, w);
    return npymath_creal(cy);
  }
}

npy_cdouble cbesh_wrap1( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int m = 1;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy99[1] = { NAN };

  npymath_csetreal(&cy, NAN);
  npymath_csetimag(&cy, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besh(z99, v, kode, m, n, cy99, &ierr);
  npymath_csetreal(&cy, creal(cy99[0]));
  npymath_csetimag(&cy, cimag(cy99[0]));
  DO_SFERR("hankel1:", &cy);
  if (sign == -1) {
    cy = rotate(cy, v);
  }
  return cy;
}

npy_cdouble cbesh_wrap1_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int m = 1;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy99[1] = { NAN };

  npymath_csetreal(&cy, NAN);
  npymath_csetimag(&cy, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besh(z99, v, kode, m, n, cy99, &ierr);
  npymath_csetreal(&cy, creal(cy99[0]));
  npymath_csetimag(&cy, cimag(cy99[0]));
  DO_SFERR("hankel1e:", &cy);
  if (sign == -1) {
    cy = rotate(cy, v);
  }
  return cy;
}

npy_cdouble cbesh_wrap2( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int m = 2;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy99[1] = { NAN };

  npymath_csetreal(&cy, NAN);
  npymath_csetimag(&cy, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besh(z99, v, kode, m, n, cy99, &ierr);
  npymath_csetreal(&cy, creal(cy99[0]));
  npymath_csetimag(&cy, cimag(cy99[0]));
  DO_SFERR("hankel2:", &cy);
  if (sign == -1) {
    cy = rotate(cy, -v);
  }
  return cy;
}

npy_cdouble cbesh_wrap2_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int m = 2;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy;

  double complex z99 = CMPLX(npymath_creal(z), npymath_cimag(z));
  double complex cy99[1] = { NAN };

  npymath_csetreal(&cy, NAN);
  npymath_csetimag(&cy, NAN);

  if (isnan(v) || isnan(npymath_creal(z)) || isnan(npymath_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besh(z99, v, kode, m, n, cy99, &ierr);
  npymath_csetreal(&cy, creal(cy99[0]));
  npymath_csetimag(&cy, cimag(cy99[0]));
  DO_SFERR("hankel2e:", &cy);
  if (sign == -1) {
    cy = rotate(cy, -v);
  }
  return cy;
}
