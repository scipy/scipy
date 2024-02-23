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
    NPY_CSETREAL(v, NAN);
    NPY_CSETIMAG(v, NAN);
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
    NPY_CSETREAL(&w, npy_creal(z)*c - npy_cimag(z)*s);
    NPY_CSETIMAG(&w, npy_creal(z)*s + npy_cimag(z)*c);
    return w;
}

static npy_cdouble
rotate_jy(npy_cdouble j, npy_cdouble y, double v)
{
    npy_cdouble w;
    double c = cos_pi(v);
    double s = sin_pi(v);
    NPY_CSETREAL(&w, npy_creal(j) * c - npy_creal(y) * s);
    NPY_CSETIMAG(&w, npy_cimag(j) * c - npy_cimag(y) * s);
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
        NPY_CSETREAL(jy, -npy_creal(*jy));
        NPY_CSETIMAG(jy, -npy_cimag(*jy));
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
    NPY_CSETREAL(&w, npy_creal(i) + s*npy_creal(k));
    NPY_CSETIMAG(&w, npy_cimag(i) + s*npy_cimag(k));
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
        NPY_CSETREAL(&z, x);
        NPY_CSETIMAG(&z, 0);
        cairy_wrap(z, &zai, &zaip, &zbi, &zbip);
        *ai  = npy_creal(zai);
        *aip = npy_creal(zaip);
        *bi  = npy_creal(zbi);
        *bip = npy_creal(zbip);
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
  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex res;
  
  NPY_CSETREAL(ai, NAN);
  NPY_CSETIMAG(ai, NAN);
  NPY_CSETREAL(bi, NAN);
  NPY_CSETIMAG(bi, NAN);
  NPY_CSETREAL(aip, NAN);
  NPY_CSETIMAG(aip, NAN);
  NPY_CSETREAL(bip, NAN);
  NPY_CSETIMAG(bip, NAN);

  res = amos_airy(z99, id, kode, &nz, &ierr);
  NPY_CSETREAL(ai, creal(res));
  NPY_CSETIMAG(ai, cimag(res));
  DO_SFERR("airy:", ai);

  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  NPY_CSETREAL(bi, creal(res));
  NPY_CSETIMAG(bi, cimag(res));
  DO_SFERR("airy:", bi);

  id = 1;
  res = amos_airy(z99, id, kode, &nz, &ierr);
  NPY_CSETREAL(aip, creal(res));
  NPY_CSETIMAG(aip, cimag(res));
  DO_SFERR("airy:", aip);

  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  NPY_CSETREAL(bip, creal(res));
  NPY_CSETIMAG(bip, cimag(res));
  DO_SFERR("airy:", bip);
  return 0;
}

int cairy_wrap_e(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
  int id = 0;
  int kode = 2;        /* Exponential scaling */
  int nz, ierr;

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex res;

  NPY_CSETREAL(ai, NAN);
  NPY_CSETIMAG(ai, NAN);
  NPY_CSETREAL(bi, NAN);
  NPY_CSETIMAG(bi, NAN);
  NPY_CSETREAL(aip, NAN);
  NPY_CSETIMAG(aip, NAN);
  NPY_CSETREAL(bip, NAN);
  NPY_CSETIMAG(bip, NAN);

  res = amos_airy(z99, id, kode, &nz, &ierr);
  NPY_CSETREAL(ai, creal(res));
  NPY_CSETIMAG(ai, cimag(res));
  DO_SFERR("airye:", ai);

  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  NPY_CSETREAL(bi, creal(res));
  NPY_CSETIMAG(bi, cimag(res));
  DO_SFERR("airye:", bi);

  id = 1;
  res = amos_airy(z99, id, kode, &nz, &ierr);
  NPY_CSETREAL(aip, creal(res));
  NPY_CSETIMAG(aip, cimag(res));
  DO_SFERR("airye:", aip);

  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  NPY_CSETREAL(bip, creal(res));
  NPY_CSETIMAG(bip, cimag(res));
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

  NPY_CSETREAL(&cai, NAN);
  NPY_CSETIMAG(&cai, NAN);
  NPY_CSETREAL(&cbi, NAN);
  NPY_CSETIMAG(&cbi, NAN);
  NPY_CSETREAL(&caip, NAN);
  NPY_CSETIMAG(&caip, NAN);
  NPY_CSETREAL(&cbip, NAN);
  NPY_CSETIMAG(&cbip, NAN);


  if (z < 0) {
      *ai = NAN;
  } else {
      res = amos_airy(z99, id, kode, &nz, &ierr);
      NPY_CSETREAL(&cai, creal(res));
      NPY_CSETIMAG(&cai, cimag(res));
      DO_SFERR("airye:", &cai);
      *ai = npy_creal(cai);
  }
  
  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  NPY_CSETREAL(&cbi, creal(res));
  NPY_CSETIMAG(&cbi, cimag(res));
  DO_SFERR("airye:", &cbi);
  *bi = npy_creal(cbi);

  id = 1;
  if (z < 0) {
      *aip = NAN;
  } else {
      res = amos_airy(z99, id, kode, &nz, &ierr);
      NPY_CSETREAL(&caip, creal(res));
      NPY_CSETIMAG(&caip, cimag(res));
      DO_SFERR("airye:", &caip);
      *aip = npy_creal(caip);
  }

  nz = 0;
  res = amos_biry(z99, id, kode, &ierr);
  NPY_CSETREAL(&cbip, creal(res));
  NPY_CSETIMAG(&cbip, cimag(res));
  DO_SFERR("airye:", &cbip);
  *bip = npy_creal(cbip);
  return 0;
}

npy_cdouble cbesi_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int sign = 1;
  int nz, ierr;
  npy_cdouble cy, cy_k;

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy99[1] = { NAN };
  double complex cy_k99[1] = { NAN };

  NPY_CSETREAL(&cy, NAN);
  NPY_CSETIMAG(&cy, NAN);
  NPY_CSETREAL(&cy_k, NAN);
  NPY_CSETIMAG(&cy_k, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besi(z99, v, kode, n, cy99, &ierr);
  NPY_CSETREAL(&cy, creal(cy99[0]));
  NPY_CSETIMAG(&cy, cimag(cy99[0]));
  DO_SFERR("iv:", &cy);
  if (ierr == 2) {
    /* overflow */
    if (npy_cimag(z) == 0 && (npy_creal(z) >= 0 || v == floor(v))) {
        if (npy_creal(z) < 0 && v/2 != floor(v/2))
            NPY_CSETREAL(&cy, -INFINITY);
        else
            NPY_CSETREAL(&cy, INFINITY);
        NPY_CSETIMAG(&cy, 0);
    } else {
        cy = cbesi_wrap_e(v*sign, z);
        NPY_CSETREAL(&cy, npy_creal(cy) * INFINITY);
        NPY_CSETIMAG(&cy, npy_cimag(cy) * INFINITY);
    }
  }

  if (sign == -1) {
    if (!reflect_i(&cy, v)) {
      nz = amos_besk(z99, v, kode, n, cy_k99, &ierr);
      NPY_CSETREAL(&cy_k, creal(cy_k99[0]));
      NPY_CSETIMAG(&cy_k, cimag(cy_k99[0]));
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

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy99[1] = { NAN };
  double complex cy_k99[1] = { NAN };

  NPY_CSETREAL(&cy, NAN);
  NPY_CSETIMAG(&cy, NAN);
  NPY_CSETREAL(&cy_k, NAN);
  NPY_CSETIMAG(&cy_k, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besi(z99, v, kode, n, cy99, &ierr);
  NPY_CSETREAL(&cy, creal(cy99[0]));
  NPY_CSETIMAG(&cy, cimag(cy99[0]));
  DO_SFERR("ive:", &cy);

  if (sign == -1) {
    if (!reflect_i(&cy, v)) {
      nz = amos_besk(z99, v, kode, n, cy_k99, &ierr);
      NPY_CSETREAL(&cy_k, creal(cy_k99[0]));
      NPY_CSETIMAG(&cy_k, cimag(cy_k99[0]));
      DO_SFERR("ive(kv):", &cy_k);
      /* adjust scaling to match zbesi */
      cy_k = rotate(cy_k, -npy_cimag(z)/M_PI);
      if (npy_creal(z) > 0) {
          NPY_CSETREAL(&cy_k, npy_creal(cy_k) * exp(-2*npy_creal(z)));
          NPY_CSETIMAG(&cy_k, npy_cimag(cy_k) * exp(-2*npy_creal(z)));
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
    NPY_CSETREAL(&w, z);
    NPY_CSETIMAG(&w, 0);
    cy = cbesi_wrap_e(v, w);
    return npy_creal(cy);
  }
}

npy_cdouble cbesj_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_j, cy_y;

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy_j99[1] = { NAN };
  double complex cy_y99[1] = { NAN };

  NPY_CSETREAL(&cy_j, NAN);
  NPY_CSETIMAG(&cy_j, NAN);
  NPY_CSETREAL(&cy_y, NAN);
  NPY_CSETIMAG(&cy_y, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy_j;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besj(z99, v, kode, n, cy_j99, &ierr);
  NPY_CSETREAL(&cy_j, creal(cy_j99[0]));
  NPY_CSETIMAG(&cy_j, cimag(cy_j99[0]));
  DO_SFERR("jv:", &cy_j);
  if (ierr == 2) {
    /* overflow */
    cy_j = cbesj_wrap_e(v, z);
    NPY_CSETREAL(&cy_j, npy_creal(cy_j) * INFINITY);
    NPY_CSETIMAG(&cy_j, npy_cimag(cy_j) * INFINITY);
  }

  if (sign == -1) {
    if (!reflect_jy(&cy_j, v)) {
      nz = amos_besy(z99, v, kode, n, cy_y99, &ierr);
      NPY_CSETREAL(&cy_y, creal(cy_y99[0]));
      NPY_CSETIMAG(&cy_y, cimag(cy_y99[0]));
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

    NPY_CSETREAL(&z, x);
    NPY_CSETIMAG(&z, 0);
    r = cbesj_wrap(v, z);
    if (npy_creal(r) != npy_creal(r)) {
        /* AMOS returned NaN, possibly due to overflow */
        return cephes_jv(v, x);
    }
    return npy_creal(r);
}

npy_cdouble cbesj_wrap_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_j, cy_y;

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy_j99[1] = { NAN };
  double complex cy_y99[1] = { NAN };

  NPY_CSETREAL(&cy_j, NAN);
  NPY_CSETIMAG(&cy_j, NAN);
  NPY_CSETREAL(&cy_y, NAN);
  NPY_CSETIMAG(&cy_y, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy_j;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besj(z99, v, kode, n, cy_j99, &ierr);
  NPY_CSETREAL(&cy_j, creal(cy_j99[0]));
  NPY_CSETIMAG(&cy_j, cimag(cy_j99[0]));
  DO_SFERR("jve:", &cy_j);
  if (sign == -1) {
    if (!reflect_jy(&cy_j, v)) {
      nz = amos_besy(z99, v, kode, n, cy_y99, &ierr);
      NPY_CSETREAL(&cy_y, creal(cy_y99[0]));
      NPY_CSETIMAG(&cy_y, cimag(cy_y99[0]));
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
    NPY_CSETREAL(&w, z);
    NPY_CSETIMAG(&w, 0);
    cy = cbesj_wrap_e(v, w);
    return npy_creal(cy);
  }
}

npy_cdouble cbesy_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_y, cy_j;

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy_j99[1] = { NAN };
  double complex cy_y99[1] = { NAN };

  NPY_CSETREAL(&cy_j, NAN);
  NPY_CSETIMAG(&cy_j, NAN);
  NPY_CSETREAL(&cy_y, NAN);
  NPY_CSETIMAG(&cy_y, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy_y;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }

  if (npy_creal(z) == 0 && npy_cimag(z) == 0) {
      /* overflow */
      NPY_CSETREAL(&cy_y, -INFINITY);
      NPY_CSETIMAG(&cy_y, 0);
      sf_error("yv", SF_ERROR_OVERFLOW, NULL);
  }
  else {
      nz = amos_besy(z99, v, kode, n, cy_y99, &ierr);
      NPY_CSETREAL(&cy_y, creal(cy_y99[0]));
      NPY_CSETIMAG(&cy_y, cimag(cy_y99[0]));
      DO_SFERR("yv:", &cy_y);
      if (ierr == 2) {
          if (npy_creal(z) >= 0 && npy_cimag(z) == 0) {
              /* overflow */
              NPY_CSETREAL(&cy_y, -INFINITY);
              NPY_CSETIMAG(&cy_y, 0);
          }
      }
  }

  if (sign == -1) {
    if (!reflect_jy(&cy_y, v)) {
      nz = amos_besj(z99, v, kode, n, cy_j99, &ierr);
      NPY_CSETREAL(&cy_j, creal(cy_j99[0]));
      NPY_CSETIMAG(&cy_j, cimag(cy_j99[0]));
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

    NPY_CSETREAL(&z, x);
    NPY_CSETIMAG(&z, 0);
    r = cbesy_wrap(v, z);
    if (npy_creal(r) != npy_creal(r)) {
        /* AMOS returned NaN, possibly due to overflow */
        return cephes_yv(v, x);
    }
    return npy_creal(r);
}

npy_cdouble cbesy_wrap_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_y, cy_j;

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy_j99[1] = { NAN };
  double complex cy_y99[1] = { NAN };

  NPY_CSETREAL(&cy_j, NAN);
  NPY_CSETIMAG(&cy_j, NAN);
  NPY_CSETREAL(&cy_y, NAN);
  NPY_CSETIMAG(&cy_y, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy_y;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besy(z99, v, kode, n, cy_y99, &ierr);
  NPY_CSETREAL(&cy_y, creal(cy_y99[0]));
  NPY_CSETIMAG(&cy_y, cimag(cy_y99[0]));
  DO_SFERR("yve:", &cy_y);
  if (ierr == 2) {
    if (npy_creal(z) >= 0 && npy_cimag(z) == 0) {
      /* overflow */
      NPY_CSETREAL(&cy_y, INFINITY);
      NPY_CSETIMAG(&cy_y, 0);
    }
  }

  if (sign == -1) {
    if (!reflect_jy(&cy_y, v)) {
      nz = amos_besj(z99, v, kode, n, cy_j99, &ierr);
      NPY_CSETREAL(&cy_j, creal(cy_j99[0]));
      NPY_CSETIMAG(&cy_j, cimag(cy_j99[0]));
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
    NPY_CSETREAL(&w, z);
    NPY_CSETIMAG(&w, 0);
    cy = cbesy_wrap_e(v, w);
    return npy_creal(cy);
  }
}

npy_cdouble cbesk_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  npy_cdouble cy;

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy99[1] = { NAN };

  NPY_CSETREAL(&cy, NAN);
  NPY_CSETIMAG(&cy, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    /* K_v == K_{-v} even for non-integer v */
    v = -v;
  }
  nz = amos_besk(z99, v, kode, n, cy99, &ierr);
  NPY_CSETREAL(&cy, creal(cy99[0]));
  NPY_CSETIMAG(&cy, cimag(cy99[0]));
  DO_SFERR("kv:", &cy);
  if (ierr == 2) {
    if (npy_creal(z) >= 0 && npy_cimag(z) == 0) {
      /* overflow */
      NPY_CSETREAL(&cy, INFINITY);
      NPY_CSETIMAG(&cy, 0);
    }
  }

  return cy;
}

npy_cdouble cbesk_wrap_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  npy_cdouble cy;

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy99[1] = { NAN };

  NPY_CSETREAL(&cy, NAN);
  NPY_CSETIMAG(&cy, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    /* K_v == K_{-v} even for non-integer v */
    v = -v;
  }
  nz = amos_besk(z99, v, kode, n, cy99, &ierr);
  NPY_CSETREAL(&cy, creal(cy99[0]));
  NPY_CSETIMAG(&cy, cimag(cy99[0]));
  DO_SFERR("kve:", &cy);
  if (ierr == 2) {
    if (npy_creal(z) >= 0 && npy_cimag(z) == 0) {
      /* overflow */
      NPY_CSETREAL(&cy, INFINITY);
      NPY_CSETIMAG(&cy, 0);
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
    NPY_CSETREAL(&w, z);
    NPY_CSETIMAG(&w, 0);
    cy = cbesk_wrap(v, w);
    return npy_creal(cy);
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
    NPY_CSETREAL(&w, z);
    NPY_CSETIMAG(&w, 0);
    cy = cbesk_wrap_e(v, w);
    return npy_creal(cy);
  }
}

npy_cdouble cbesh_wrap1( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int m = 1;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy;

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy99[1] = { NAN };

  NPY_CSETREAL(&cy, NAN);
  NPY_CSETIMAG(&cy, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besh(z99, v, kode, m, n, cy99, &ierr);
  NPY_CSETREAL(&cy, creal(cy99[0]));
  NPY_CSETIMAG(&cy, cimag(cy99[0]));
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

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy99[1] = { NAN };

  NPY_CSETREAL(&cy, NAN);
  NPY_CSETIMAG(&cy, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besh(z99, v, kode, m, n, cy99, &ierr);
  NPY_CSETREAL(&cy, creal(cy99[0]));
  NPY_CSETIMAG(&cy, cimag(cy99[0]));
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

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy99[1] = { NAN };

  NPY_CSETREAL(&cy, NAN);
  NPY_CSETIMAG(&cy, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besh(z99, v, kode, m, n, cy99, &ierr);
  NPY_CSETREAL(&cy, creal(cy99[0]));
  NPY_CSETIMAG(&cy, cimag(cy99[0]));
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

  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex cy99[1] = { NAN };

  NPY_CSETREAL(&cy, NAN);
  NPY_CSETIMAG(&cy, NAN);

  if (isnan(v) || isnan(npy_creal(z)) || isnan(npy_cimag(z))) {
    return cy;
  }
  if (v < 0) {
    v = -v;
    sign = -1;
  }
  nz = amos_besh(z99, v, kode, m, n, cy99, &ierr);
  NPY_CSETREAL(&cy, creal(cy99[0]));
  NPY_CSETIMAG(&cy, cimag(cy99[0]));
  DO_SFERR("hankel2e:", &cy);
  if (sign == -1) {
    cy = rotate(cy, -v);
  }
  return cy;
}
