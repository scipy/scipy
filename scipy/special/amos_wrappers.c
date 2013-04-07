/* This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 */

#include "amos_wrappers.h"

#define CADDR(z) (double *)(&((z).real)), (double*)(&((z).imag))
#define F2C_CST(z) (double *)&((z)->real), (double *)&((z)->imag)

#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif

extern int F_FUNC(zairy,ZAIRY)
     (double*, double*, int*, int*, double*, double*, int*, int*);
extern int F_FUNC(zbiry,ZBIRY)
     (double*, double*, int*, int*, double*, double*, int*);
extern int F_FUNC(zbesi,ZBESI)
     (double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern int F_FUNC(zbesj,ZBESJ)
     (double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern int F_FUNC(zbesk,ZBESK)
     (double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern int F_FUNC(zbesy,ZBESY)
     (double*, double*, double*, int*, int*, double*, double*, int*, double*, double*, int*);
extern int F_FUNC(zbesh,ZBESH)
     (double*, double*, double*, int*, int*, int*, double*, double*, int*, int*);

/* This must be linked with fortran
 */

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
    v->real = NPY_NAN;
    v->imag = NPY_NAN;
  }
}

static npy_cdouble
rotate(npy_cdouble z, double v)
{
    npy_cdouble w;
    double c = cos(v * NPY_PI);
    double s = sin(v * NPY_PI);
    w.real = z.real*c - z.imag*s;
    w.imag = z.real*s + z.imag*c;
    return w;
}

static npy_cdouble
rotate_jy(npy_cdouble j, npy_cdouble y, double v)
{
    npy_cdouble w;
    double c = cos(v * NPY_PI);
    double s = sin(v * NPY_PI);
    w.real = j.real * c - y.real * s;
    w.imag = j.imag * c - y.imag * s;
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
        jy->real = -jy->real;
        jy->imag = -jy->imag;
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
    double s = sin(v * NPY_PI)*(2.0/NPY_PI);
    w.real = i.real + s*k.real;
    w.imag = i.imag + s*k.imag;
    return w;
}

int cairy_wrap(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
  int id = 0;
  int ierr = 0;
  int kode = 1;
  int nz;

  F_FUNC(zairy,ZAIRY)(CADDR(z), &id, &kode, F2C_CST(ai), &nz, &ierr);
  DO_SFERR("airy:", ai);
  nz = 0;
  F_FUNC(zbiry,ZBIRY)(CADDR(z), &id, &kode, F2C_CST(bi), &ierr);
  DO_SFERR("airy:", bi);
  
  id = 1;
  F_FUNC(zairy,ZAIRY)(CADDR(z), &id, &kode, F2C_CST(aip), &nz, &ierr);
  DO_SFERR("airy:", aip);
  nz = 0;
  F_FUNC(zbiry,ZBIRY)(CADDR(z), &id, &kode, F2C_CST(bip), &ierr);
  DO_SFERR("airy:", bip);
  return 0;
}

int cairy_wrap_e(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
  int id = 0;
  int kode = 2;        /* Exponential scaling */
  int nz, ierr;

  F_FUNC(zairy,ZAIRY)(CADDR(z), &id, &kode, F2C_CST(ai), &nz, &ierr);
  DO_SFERR("airye:", ai);
  nz = 0;
  F_FUNC(zbiry,ZBIRY)(CADDR(z), &id, &kode, F2C_CST(bi), &ierr);
  DO_SFERR("airye:", bi);
  
  id = 1;
  F_FUNC(zairy,ZAIRY)(CADDR(z), &id, &kode, F2C_CST(aip), &nz, &ierr);
  DO_SFERR("airye:", aip);
  nz = 0;
  F_FUNC(zbiry,ZBIRY)(CADDR(z), &id, &kode, F2C_CST(bip), &ierr);
  DO_SFERR("airye:", bip);
  return 0;
}

int cairy_wrap_e_real(double z, double *ai, double *aip, double *bi, double *bip) {
  int id = 0;
  int kode = 2;        /* Exponential scaling */
  int nz, ierr;
  npy_cdouble cz, cai, caip, cbi, cbip;

  cz.real = z;
  cz.imag = 0;

  if (z < 0) {
      *ai = NPY_NAN;
  } else {
      F_FUNC(zairy,ZAIRY)(CADDR(cz), &id, &kode, CADDR(cai), &nz, &ierr);
      DO_SFERR("airye:", &cai);
      *ai = cai.real;
  }
  nz = 0;
  F_FUNC(zbiry,ZBIRY)(CADDR(cz), &id, &kode, CADDR(cbi), &ierr);
  DO_SFERR("airye:", &cbi);
  *bi = cbi.real;
  
  id = 1;
  if (z < 0) {
      *aip = NPY_NAN;
  } else {
      F_FUNC(zairy,ZAIRY)(CADDR(cz), &id, &kode, CADDR(caip), &nz, &ierr);
      DO_SFERR("airye:", &caip);
      *aip = caip.real;
  }
  nz = 0;
  F_FUNC(zbiry,ZBIRY)(CADDR(cz), &id, &kode, CADDR(cbip), &ierr);
  DO_SFERR("airye:", &cbip);
  *bip = cbip.real;
  return 0;
}

npy_cdouble cbesi_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int sign = 1;
  int nz, ierr;
  npy_cdouble cy, cy_k;

  if (v < 0) {
    v = -v;
    sign = -1;
  }
  F_FUNC(zbesi,ZBESI)(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_SFERR("iv:", &cy);
  if (ierr == 2) {
    /* overflow */
    if (z.imag == 0 && (z.real >= 0 || v == floor(v))) {
        if (z.real < 0 && v/2 != floor(v/2))
            cy.real = -NPY_INFINITY;
        else
            cy.real = NPY_INFINITY;
        cy.imag = 0;
    } else {
        cy = cbesi_wrap_e(v*sign, z);
        cy.real *= NPY_INFINITY;
        cy.imag *= NPY_INFINITY;
    }
  }

  if (sign == -1) {
    if (!reflect_i(&cy, v)) {
      F_FUNC(zbesk,ZBESK)(CADDR(z), &v,  &kode, &n, CADDR(cy_k), &nz, &ierr);
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

  if (v < 0) {
    v = -v;
    sign = -1;
  }
  F_FUNC(zbesi,ZBESI)(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_SFERR("ive:", &cy);

  if (sign == -1) {
    if (!reflect_i(&cy, v)) {
      F_FUNC(zbesk,ZBESK)(CADDR(z), &v,  &kode, &n, CADDR(cy_k), &nz, &ierr);
      DO_SFERR("ive(kv):", &cy_k);
      /* adjust scaling to match zbesi */
      cy_k = rotate(cy_k, -z.imag/NPY_PI);
      if (z.real > 0) {
          cy_k.real *= exp(-2*z.real);
          cy_k.imag *= exp(-2*z.real);
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
    return NPY_NAN;
  } else {
    w.real = z;
    w.imag = 0;
    cy = cbesi_wrap_e(v, w);
    return cy.real;
  }
}
  
npy_cdouble cbesj_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_j, cy_y, cwork;

  if (v < 0) {
    v = -v;
    sign = -1;
  }
  F_FUNC(zbesj,ZBESJ)(CADDR(z), &v,  &kode, &n, CADDR(cy_j), &nz, &ierr);
  DO_SFERR("jv:", &cy_j);
  if (ierr == 2) {
    /* overflow */
    cy_j = cbesj_wrap_e(v, z);
    cy_j.real *= NPY_INFINITY;
    cy_j.imag *= NPY_INFINITY;
  }

  if (sign == -1) {
    if (!reflect_jy(&cy_j, v)) {
      F_FUNC(zbesy,ZBESY)(CADDR(z), &v,  &kode, &n, CADDR(cy_y), &nz, CADDR(cwork), &ierr);
      DO_SFERR("jv(yv):", &cy_y);
      cy_j = rotate_jy(cy_j, cy_y, v);
    }
  }
  return cy_j;
}

npy_cdouble cbesj_wrap_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_j, cy_y, cwork;

  if (v < 0) {
    v = -v;
    sign = -1;
  }
  F_FUNC(zbesj,ZBESJ)(CADDR(z), &v, &kode, &n, CADDR(cy_j), &nz, &ierr);
  DO_SFERR("jve:", &cy_j);
  if (sign == -1) {
    if (!reflect_jy(&cy_j, v)) {
      F_FUNC(zbesy,ZBESY)(CADDR(z), &v,  &kode, &n, CADDR(cy_y), &nz, CADDR(cwork), &ierr);
      DO_SFERR("jve(yve):", &cy_y);
      cy_j = rotate_jy(cy_j, cy_y, v);
    }
  }
  return cy_j;
}

double cbesj_wrap_e_real(double v, double z) {
  npy_cdouble cy, w;
  if (v != floor(v) && z < 0) {
    return NPY_NAN;
  } else {
    w.real = z;
    w.imag = 0;
    cy = cbesj_wrap_e(v, w);
    return cy.real;
  }
}
  
npy_cdouble cbesy_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_y, cy_j, cwork;

  if (v < 0) {
    v = -v;
    sign = -1;
  }
  F_FUNC(zbesy,ZBESY)(CADDR(z), &v,  &kode, &n, CADDR(cy_y), &nz, CADDR(cwork), &ierr);
  DO_SFERR("yv:", &cy_y);
  if (ierr == 2) {
    if (z.real >= 0 && z.imag == 0) {
      /* overflow */
      cy_y.real = NPY_INFINITY;
      cy_y.imag = 0;
    }
  }

  if (sign == -1) {
    if (!reflect_jy(&cy_y, v)) {
      F_FUNC(zbesj,ZBESJ)(CADDR(z), &v,  &kode, &n, CADDR(cy_j), &nz, &ierr);
      DO_SFERR("yv(jv):", &cy_j);
      cy_y = rotate_jy(cy_y, cy_j, -v);
    }
  }
  return cy_y;
}

npy_cdouble cbesy_wrap_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy_y, cy_j, cwork;

  if (v < 0) {
    v = -v;
    sign = -1;
  }
  F_FUNC(zbesy,ZBESY)(CADDR(z), &v, &kode, &n, CADDR(cy_y), &nz, CADDR(cwork), &ierr);
  DO_SFERR("yve:", &cy_y);
  if (ierr == 2) {
    if (z.real >= 0 && z.imag == 0) {
      /* overflow */
      cy_y.real = NPY_INFINITY;
      cy_y.imag = 0;
    }
  }

  if (sign == -1) {
    if (!reflect_jy(&cy_y, v)) {
      F_FUNC(zbesj,ZBESJ)(CADDR(z), &v,  &kode, &n, CADDR(cy_j), &nz, &ierr);
      DO_SFERR("yv(jv):", &cy_j);
      cy_y = rotate_jy(cy_y, cy_j, -v);
    }
  }
  return cy_y;
}

double cbesy_wrap_e_real(double v, double z) {
  npy_cdouble cy, w;
  if (z < 0) {
    return NPY_NAN;
  } else {
    w.real = z;
    w.imag = 0;
    cy = cbesy_wrap_e(v, w);
    return cy.real;
  }
}
  
npy_cdouble cbesk_wrap( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  npy_cdouble cy;

  if (v < 0) {
    /* K_v == K_{-v} even for non-integer v */
    v = -v;
  }
  F_FUNC(zbesk,ZBESK)(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_SFERR("kv:", &cy);
  if (ierr == 2) {
    if (z.real >= 0 && z.imag == 0) {
      /* overflow */
      cy.real = NPY_INFINITY;
      cy.imag = 0;
    }
  }

  return cy;
}

npy_cdouble cbesk_wrap_e( double v, npy_cdouble z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  npy_cdouble cy;

  if (v < 0) {
    /* K_v == K_{-v} even for non-integer v */
    v = -v;
  }
  F_FUNC(zbesk,ZBESK)(CADDR(z), &v, &kode, &n, CADDR(cy), &nz, &ierr);
  DO_SFERR("kve:", &cy);
  if (ierr == 2) {
    if (z.real >= 0 && z.imag == 0) {
      /* overflow */
      cy.real = NPY_INFINITY;
      cy.imag = 0;
    }
  }

  return cy;
}
  
double cbesk_wrap_real( double v, double z) {
  npy_cdouble cy, w;
  if (z < 0) {
    return NPY_NAN;
  }
  else if (z > 710 * (1 + fabs(v))) {
      /* Underflow. See uniform expansion http://dlmf.nist.gov/10.41
       * This condition is not a strict bound (it can underflow earlier),
       * rather, we are here working around a restriction in AMOS.
       */
      return 0;
  }
  else {
    w.real = z;
    w.imag = 0;
    cy = cbesk_wrap(v, w);
    return cy.real;
  }
}

double cbesk_wrap_real_int(int n, double z)
{
    return cbesk_wrap_real(n, z);
}

double cbesk_wrap_e_real( double v, double z) {
  npy_cdouble cy, w;
  if (z < 0) {
    return NPY_NAN;
  } else {
    w.real = z;
    w.imag = 0;
    cy = cbesk_wrap_e(v, w);
    return cy.real;
  }
}
  
npy_cdouble cbesh_wrap1( double v, npy_cdouble z) {
  int n = 1;
  int kode = 1;
  int m = 1;
  int nz, ierr;
  int sign = 1;
  npy_cdouble cy;

  if (v < 0) {
    v = -v;
    sign = -1;
  }
  F_FUNC(zbesh,ZBESH)(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
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

  if (v < 0) {
    v = -v;
    sign = -1;
  }
  F_FUNC(zbesh,ZBESH)(CADDR(z), &v, &kode, &m, &n, CADDR(cy), &nz, &ierr);
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

  if (v < 0) {
    v = -v;
    sign = -1;
  }
  F_FUNC(zbesh,ZBESH)(CADDR(z), &v, &kode, &m, &n, CADDR(cy), &nz, &ierr);
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

  if (v < 0) {
    v = -v;
    sign = -1;
  }
  F_FUNC(zbesh,ZBESH)(CADDR(z), &v, &kode, &m, &n, CADDR(cy), &nz, &ierr);
  DO_SFERR("hankel2e:", &cy);
  if (sign == -1) {
    cy = rotate(cy, -v);
  }
  return cy;
}
