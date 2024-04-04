#pragma once

#include "amos/amos.h"
#include "config.h"
#include "error.h"

#include "cephes/airy.h"
#include "cephes/jv.h"
#include "cephes/yv.h"

inline int cephes_airy(float xf, float *aif, float *aipf, float *bif, float *bipf) {
    double ai;
    double aip;
    double bi;
    double bip;
    int res = cephes_airy(xf, &ai, &aip, &bi, &bip);

    *aif = ai;
    *aipf = aip;
    *bif = bi;
    *bipf = bip;
    return res;
}

namespace special {

inline int ierr_to_sferr(int nz, int ierr) {
    /* Return sf_error equivalents for ierr values */

    if (nz != 0)
        return SF_ERROR_UNDERFLOW;

    switch (ierr) {
    case 1:
        return SF_ERROR_DOMAIN;
    case 2:
        return SF_ERROR_OVERFLOW;
    case 3:
        return SF_ERROR_LOSS;
    case 4:
        return SF_ERROR_NO_RESULT;
    case 5: /* Algorithm termination condition not met */
        return SF_ERROR_NO_RESULT;
    }
    return -1;
}

template <typename T>
void set_nan_if_no_computation_done(std::complex<T> *v, int ierr) {
    if (v != NULL && (ierr == 1 || ierr == 2 || ierr == 4 || ierr == 5)) {
        v->real(NAN);
        v->imag(NAN);
    }
}

template <typename T>
void do_sferr(const char *name, std::complex<T> *ai, int nz, int ierr) {
    if (nz != 0 || ierr != 0) {
        set_error(name, (sf_error_t) ierr_to_sferr(nz, ierr), NULL);
        set_nan_if_no_computation_done(ai, ierr);
    }
}

inline double sin_pi(double x) {
    if (floor(x) == x && fabs(x) < 1e14) {
        /* Return 0 when at exact zero, as long as the floating point number is
         * small enough to distinguish integer points from other points.
         */
        return 0;
    }
    return sin(M_PI * x);
}

inline double cos_pi(double x) {
    double x05 = x + 0.5;
    if (floor(x05) == x05 && fabs(x) < 1e14) {
        /* Return 0 when at exact zero, as long as the floating point number is
         * small enough to distinguish integer points from other points.
         */
        return 0;
    }
    return cos(M_PI * x);
}

inline std::complex<double> rotate(std::complex<double> z, double v) {
    std::complex<double> w;
    double c = cos_pi(v);
    double s = sin_pi(v);
    w.real(z.real() * c - z.imag() * s);
    w.imag(z.real() * s + z.imag() * c);

    return w;
}

inline std::complex<double> rotate_jy(std::complex<double> j, std::complex<double> y, double v) {
    std::complex<double> w;
    double c = cos_pi(v);
    double s = sin_pi(v);
    w.real(j.real() * c - y.real() * s);
    w.imag(j.imag() * c - y.imag() * s);
    return w;
}

inline int reflect_jy(std::complex<double> *jy, double v) {
    /* NB: Y_v may be huge near negative integers -- so handle exact
     *     integers carefully
     */
    int i;
    if (v != floor(v))
        return 0;

    i = v - 16384.0 * floor(v / 16384.0);
    if (i & 1) {
        jy->real(-jy->real());
        jy->imag(-jy->imag());
    }
    return 1;
}

inline int reflect_i(std::complex<double> *ik, double v) {
    if (v != floor(v))
        return 0;
    return 1; /* I is symmetric for integer v */
}

inline std::complex<double> rotate_i(std::complex<double> i, std::complex<double> k, double v) {
    std::complex<double> w;
    double s = sin(v * M_PI) * (2.0 / M_PI);
    w.real(i.real() + s * k.real());
    w.imag(i.imag() + s * k.imag());
    return w;
}

<<<<<<< HEAD
template <typename T>
void cairy_wrap(std::complex<T> z, std::complex<T> *ai, std::complex<T> *aip, std::complex<T> *bi,
                std::complex<T> *bip) {
=======
inline int cairy_wrap(std::complex<double> z, std::complex<double> *ai, std::complex<double> *aip, std::complex<double> *bi,
               std::complex<double> *bip) {
>>>>>>> a20d809ca6 (cephes translation is nearing completion)
    int id = 0;
    int ierr = 0;
    int kode = 1;
    int nz;
    std::complex<double> z99 = z;
    std::complex<double> res;

    ai->real(NAN);
    ai->imag(NAN);
    bi->real(NAN);
    bi->imag(NAN);
    aip->real(NAN);
    aip->imag(NAN);
    bip->real(NAN);
    bip->imag(NAN);

    res = amos::airy(z99, id, kode, &nz, &ierr);
    ai->real(res.real());
    ai->imag(res.imag());
    do_sferr("airy:", ai, nz, ierr);

    nz = 0;
    res = amos::biry(z99, id, kode, &ierr);
    bi->real(res.real());
    bi->imag(res.imag());
    do_sferr("airy:", bi, nz, ierr);

    id = 1;
    res = amos::airy(z99, id, kode, &nz, &ierr);
    aip->real(res.real());
    aip->imag(res.imag());
    do_sferr("airy:", aip, nz, ierr);

    nz = 0;
    res = amos::biry(z99, id, kode, &ierr);
    bip->real(res.real());
    bip->imag(res.imag());
    do_sferr("airy:", bip, nz, ierr);
    return 0;
}

<<<<<<< HEAD
template <typename T>
void cairy_wrap_e(std::complex<T> z, std::complex<T> *ai, std::complex<T> *aip, std::complex<T> *bi,
                  std::complex<T> *bip) {
=======
inline int cairy_wrap_e(std::complex<double> z, std::complex<double> *ai, std::complex<double> *aip, std::complex<double> *bi,
                 std::complex<double> *bip) {
>>>>>>> a20d809ca6 (cephes translation is nearing completion)
    int id = 0;
    int kode = 2; /* Exponential scaling */
    int nz, ierr;

    std::complex<double> z99 = z;
    std::complex<double> res;

    ai->real(NAN);
    ai->imag(NAN);
    bi->real(NAN);
    bi->imag(NAN);
    aip->real(NAN);
    aip->imag(NAN);
    bip->real(NAN);
    bip->imag(NAN);

    res = amos::airy(z99, id, kode, &nz, &ierr);
    ai->real(res.real());
    ai->imag(std::imag(res));
    do_sferr("airye:", ai, nz, ierr);

    nz = 0;
    res = amos::biry(z99, id, kode, &ierr);
    bi->real(res.real());
    bi->imag(std::imag(res));
    do_sferr("airye:", bi, nz, ierr);

    id = 1;
    res = amos::airy(z99, id, kode, &nz, &ierr);
    aip->real(res.real());
    aip->imag(std::imag(res));
    do_sferr("airye:", aip, nz, ierr);

    nz = 0;
    res = amos::biry(z99, id, kode, &ierr);
    bip->real(res.real());
    bip->imag(std::imag(res));
    do_sferr("airye:", bip, nz, ierr);
    return 0;
}

<<<<<<< HEAD
template <typename T>
inline void cairy_wrap_e_real(T z, T *ai, T *aip, T *bi, T *bip) {
=======
inline int cairy_wrap_e_real(double z, double *ai, double *aip, double *bi, double *bip) {
>>>>>>> a20d809ca6 (cephes translation is nearing completion)
    int id = 0;
    int kode = 2; /* Exponential scaling */
    int nz, ierr;
    std::complex<T> cai, caip, cbi, cbip;

    std::complex<double> z99 = z;
    std::complex<double> res;

    cai.real(NAN);
    cai.imag(NAN);
    cbi.real(NAN);
    cbi.imag(NAN);
    caip.real(NAN);
    caip.imag(NAN);
    cbip.real(NAN);
    cbip.real(NAN);

    if (z < 0) {
        *ai = NAN;
    } else {
        res = amos::airy(z99, id, kode, &nz, &ierr);
        cai.real(res.real());
        cai.imag(res.imag());
        do_sferr("airye:", &cai, nz, ierr);
        *ai = cai.real();
    }

    nz = 0;
    res = amos::biry(z99, id, kode, &ierr);
    cbi.real(res.real());
    cbi.imag(res.imag());
    do_sferr("airye:", &cbi, nz, ierr);
    *bi = std::real(cbi);

    id = 1;
    if (z < 0) {
        *aip = NAN;
    } else {
        res = amos::airy(z99, id, kode, &nz, &ierr);
        caip.real(res.real());
        caip.imag(res.imag());
        do_sferr("airye:", &caip, nz, ierr);
        *aip = std::real(caip);
    }

    nz = 0;
    res = amos::biry(z99, id, kode, &ierr);
    cbip.real(res.real());
    cbip.imag(res.imag());
    do_sferr("airye:", &cbip, nz, ierr);
    *bip = cbip.real();
    return 0;
}

<<<<<<< HEAD
template <typename T>
void airy_wrap(T x, T *ai, T *aip, T *bi, T *bip) {
    std::complex<T> z, zai, zaip, zbi, zbip;
=======
inline int airy_wrap(double x, double *ai, double *aip, double *bi, double *bip) {
    std::complex<double> z, zai, zaip, zbi, zbip;
>>>>>>> a20d809ca6 (cephes translation is nearing completion)

    /* For small arguments, use Cephes as it's slightly faster.
     * For large arguments, use AMOS as it's more accurate.
     */
    if (x < -10 || x > 10) {
        z.real(x);
        z.imag(0);
        cairy_wrap(z, &zai, &zaip, &zbi, &zbip);
        *ai = zai.real();
        *aip = zaip.real();
        *bi = zbi.real();
        *bip = zbip.real();
    } else {
        cephes::airy(x, ai, aip, bi, bip);
    }
    return 0;
}

inline std::complex<double> cbesi_wrap_e(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int sign = 1;
    int nz, ierr;
    std::complex<double> cy, cy_k;

    std::complex<double> z99 = z;
    std::complex<double> cy99[1] = {NAN};
    std::complex<double> cy_k99[1] = {NAN};

    cy.real(NAN);
    cy.imag(NAN);
    cy_k.real(NAN);
    cy_k.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besi(z99, v, kode, n, cy99, &ierr);
    cy.real(cy99[0].real());
    cy.imag(std::imag(cy99[0]));
    do_sferr("ive:", &cy, nz, ierr);

    if (sign == -1) {
        if (!reflect_i(&cy, v)) {
            nz = amos::besk(z99, v, kode, n, cy_k99, &ierr);
            cy_k.real(cy_k99[0].real());
            cy_k.imag(std::imag(cy_k99[0]));
            do_sferr("ive(kv):", &cy_k, nz, ierr);
            /* adjust scaling to match zbesi */
            cy_k = rotate(cy_k, -z.imag() / M_PI);
            if (z.real() > 0) {
                cy_k.real(cy_k.real() * exp(-2 * z.real()));
                cy_k.imag(cy_k.imag() * exp(-2 * z.real()));
            }
            /* v -> -v */
            cy = rotate_i(cy, cy_k, v);
        }
    }

    return cy;
}

inline double cbesi_wrap_e_real(double v, double z) {
    std::complex<double> cy, w;
    if (v != floor(v) && z < 0) {
        return NAN;
    } else {
        w.real(z);
        w.imag(0);
        cy = cbesi_wrap_e(v, w);
        return cy.real();
    }
}

inline float cbesi_wrap_e_real(float v, float z) {
    return cbesi_wrap_e_real(static_cast<double>(v), static_cast<double>(z));
}

inline std::complex<float> cbesi_wrap_e(float v, std::complex<float> z) {
    std::complex<double> res = cbesi_wrap_e(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline std::complex<double> cbesi_wrap(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int sign = 1;
    int nz, ierr;
    std::complex<double> cy, cy_k;

    std::complex<double> z99 = z;
    std::complex<double> cy99[1] = {NAN};
    std::complex<double> cy_k99[1] = {NAN};

    cy.real(NAN);
    cy.imag(NAN);
    cy_k.real(NAN);
    cy_k.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besi(z99, v, kode, n, cy99, &ierr);
    cy.real(cy99[0].real());
    cy.imag(cy99[0].imag());
    do_sferr("iv:", &cy, nz, ierr);
    if (ierr == 2) {
        /* overflow */
        if (z.imag() == 0 && (z.real() >= 0 || v == floor(v))) {
            if (z.real() < 0 && v / 2 != floor(v / 2))
                cy.real(-INFINITY);
            else
                cy.real(INFINITY);
            cy.imag(0);
        } else {
            cy = cbesi_wrap_e(v * sign, z);
            cy.real(cy.real() * INFINITY);
            cy.imag(cy.imag() * INFINITY);
        }
    }

    if (sign == -1) {
        if (!reflect_i(&cy, v)) {
            nz = amos::besk(z99, v, kode, n, cy_k99, &ierr);
            cy_k.real(cy_k99[0].real());
            cy_k.imag(cy_k99[0].imag());
            do_sferr("iv(kv):", &cy_k, nz, ierr);
            cy = rotate_i(cy, cy_k, v);
        }
    }

    return cy;
}

inline float cbesi_wrap_real(float v, float z) { return cephes_iv(static_cast<double>(v), static_cast<double>(z)); }

inline std::complex<float> cbesi_wrap(float v, std::complex<float> z) {
    std::complex<double> res = cbesi_wrap(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline std::complex<double> cbesj_wrap_e(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy_j, cy_y;

    std::complex<double> z99 = z;
    std::complex<double> cy_j99[1] = {NAN};
    std::complex<double> cy_y99[1] = {NAN};

    cy_j.real(NAN);
    cy_j.imag(NAN);
    cy_y.real(NAN);
    cy_y.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy_j;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besj(z99, v, kode, n, cy_j99, &ierr);
    cy_j.real(cy_j99[0].real());
    cy_j.imag(cy_j99[0].imag());
    do_sferr("jve:", &cy_j, nz, ierr);
    if (sign == -1) {
        if (!reflect_jy(&cy_j, v)) {
            nz = amos::besy(z99, v, kode, n, cy_y99, &ierr);
            cy_y.real(cy_y99[0].real());
            cy_y.imag(cy_y99[0].imag());
            do_sferr("jve(yve):", &cy_y, nz, ierr);
            cy_j = rotate_jy(cy_j, cy_y, v);
        }
    }
    return cy_j;
}

inline std::complex<double> cbesj_wrap(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy_j, cy_y;

    std::complex<double> z99 = z;
    std::complex<double> cy_j99[1] = {NAN};
    std::complex<double> cy_y99[1] = {NAN};

    cy_j.real(NAN);
    cy_j.imag(NAN);
    cy_y.real(NAN);
    cy_y.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy_j;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besj(z99, v, kode, n, cy_j99, &ierr);
    cy_j.real(cy_j99[0].real());
    cy_j.imag(cy_j99[0].imag());
    do_sferr("jv:", &cy_j, nz, ierr);
    if (ierr == 2) {
        /* overflow */
        cy_j = cbesj_wrap_e(v, z);
        cy_j.real(cy_j.real() * INFINITY);
        cy_j.imag(cy_j.imag() * INFINITY);
    }

    if (sign == -1) {
        if (!reflect_jy(&cy_j, v)) {
            nz = amos::besy(z99, v, kode, n, cy_y99, &ierr);
            cy_y.real(cy_y99[0].real());
            cy_y.imag(std::imag(cy_y99[0]));
            do_sferr("jv(yv):", &cy_y, nz, ierr);
            cy_j = rotate_jy(cy_j, cy_y, v);
        }
    }
    return cy_j;
}

inline double cbesj_wrap_e_real(double v, double z) {
    std::complex<double> cy, w;
    if (v != floor(v) && z < 0) {
        return NAN;
    } else {
        w.real(z);
        w.imag(0);
        cy = cbesj_wrap_e(v, w);
        return cy.real();
    }
}

inline float cbesj_wrap_e_real(float v, float z) {
    return cbesj_wrap_e_real(static_cast<double>(v), static_cast<double>(z));
}

inline std::complex<float> cbesj_wrap_e(float v, std::complex<float> z) {
    std::complex<double> res = cbesj_wrap_e(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline double cbesj_wrap_real(double v, double x) {
    std::complex<double> z, r;

    if (x < 0 && v != (int) v) {
        set_error("yv", SF_ERROR_DOMAIN, NULL);
        return NAN;
    }

    z.real(x);
    z.imag(0);
    r = cbesj_wrap(v, z);
    if (r.real() != r.real()) {
        /* AMOS returned NaN, possibly due to overflow */
        return cephes::jv(v, x);
    }
    return r.real();
}

inline float cbesj_wrap_real(float v, float z) {
    return cbesj_wrap_real(static_cast<double>(v), static_cast<double>(z));
}

inline std::complex<float> cbesj_wrap(float v, std::complex<float> z) {
    std::complex<double> res = cbesj_wrap(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline std::complex<double> cbesy_wrap(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy_y, cy_j;

    std::complex<double> z99 = z;
    std::complex<double> cy_j99[1] = {NAN};
    std::complex<double> cy_y99[1] = {NAN};

    cy_j.real(NAN);
    cy_j.imag(NAN);
    cy_y.real(NAN);
    cy_y.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy_y;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }

    if (z.real() == 0 && z.imag() == 0) {
        /* overflow */
        cy_y.real(-INFINITY);
        cy_y.imag(0);
        set_error("yv", SF_ERROR_OVERFLOW, NULL);
    } else {
        nz = amos::besy(z99, v, kode, n, cy_y99, &ierr);
        cy_y.real(cy_y99[0].real());
        cy_y.imag(std::imag(cy_y99[0]));
        do_sferr("yv:", &cy_y, nz, ierr);
        if (ierr == 2) {
            if (z.real() >= 0 && z.imag() == 0) {
                /* overflow */
                cy_y.real(-INFINITY);
                cy_y.imag(0);
            }
        }
    }

    if (sign == -1) {
        if (!reflect_jy(&cy_y, v)) {
            nz = amos::besj(z99, v, kode, n, cy_j99, &ierr);
            cy_j.real(cy_j99[0].real());
            cy_j.imag(std::imag(cy_j99[0]));
            // F_FUNC(zbesj,ZBESJ)(CADDR(z), &v,  &kode, &n, CADDR(cy_j), &nz, &ierr);
            do_sferr("yv(jv):", &cy_j, nz, ierr);
            cy_y = rotate_jy(cy_y, cy_j, -v);
        }
    }
    return cy_y;
}

inline double cbesy_wrap_real(double v, double x) {
    std::complex<double> z, r;

    if (x < 0.0) {
        set_error("yv", SF_ERROR_DOMAIN, NULL);
        return NAN;
    }

    z.real(x);
    z.imag(0);
    r = cbesy_wrap(v, z);
    if (r.real() != r.real()) {
        /* AMOS returned NaN, possibly due to overflow */
        return cephes::yv(v, x);
    }
    return r.real();
}

inline float cbesy_wrap_real(float v, float z) {
    return cbesy_wrap_real(static_cast<double>(v), static_cast<double>(z));
}

inline std::complex<float> cbesy_wrap(float v, std::complex<float> z) {
    std::complex<double> res = cbesy_wrap(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline std::complex<double> cbesy_wrap_e(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy_y, cy_j;

    std::complex<double> z99 = z;
    std::complex<double> cy_j99[1] = {NAN};
    std::complex<double> cy_y99[1] = {NAN};

    cy_j.real(NAN);
    cy_j.imag(NAN);
    cy_y.real(NAN);
    cy_y.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy_y;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besy(z99, v, kode, n, cy_y99, &ierr);
    cy_y.real(cy_y99[0].real());
    cy_y.imag(cy_y99[0].imag());
    do_sferr("yve:", &cy_y, nz, ierr);
    if (ierr == 2) {
        if (z.real() >= 0 && z.imag() == 0) {
            /* overflow */
            cy_y.real(INFINITY);
            cy_y.imag(0);
        }
    }

    if (sign == -1) {
        if (!reflect_jy(&cy_y, v)) {
            nz = amos::besj(z99, v, kode, n, cy_j99, &ierr);
            cy_j.real(cy_j99[0].real());
            cy_j.imag(cy_j99[0].imag());
            do_sferr("yv(jv):", &cy_j, nz, ierr);
            cy_y = rotate_jy(cy_y, cy_j, -v);
        }
    }
    return cy_y;
}

inline std::complex<float> cbesy_wrap_e(float v, std::complex<float> z) {
    std::complex<double> res = cbesy_wrap_e(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline double cbesy_wrap_e_real(double v, double z) {
    std::complex<double> cy, w;
    if (z < 0) {
        return NAN;
    } else {
        w.real(z);
        w.imag(0);
        cy = cbesy_wrap_e(v, w);
        return cy.real();
    }
}

inline float cbesy_wrap_e_real(float v, float z) {
    return cbesy_wrap_e_real(static_cast<double>(v), static_cast<double>(z));
}

inline std::complex<double> cbesk_wrap(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int nz, ierr;
    std::complex<double> cy;

    std::complex<double> z99 = z;
    std::complex<double> cy99[1] = {NAN};

    cy.real(NAN);
    cy.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        /* K_v == K_{-v} even for non-integer v */
        v = -v;
    }
    nz = amos::besk(z99, v, kode, n, cy99, &ierr);
    cy.real(std::real(cy99[0]));
    cy.imag(std::imag(cy99[0]));
    do_sferr("kv:", &cy, nz, ierr);
    if (ierr == 2) {
        if (z.real() >= 0 && z.imag() == 0) {
            /* overflow */
            cy.real(INFINITY);
            cy.imag(0);
        }
    }

    return cy;
}

inline std::complex<double> cbesk_wrap_e(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int nz, ierr;
    std::complex<double> cy;

    std::complex<double> z99 = z;
    std::complex<double> cy99[1] = {NAN};

    cy.real(NAN);
    cy.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        /* K_v == K_{-v} even for non-integer v */
        v = -v;
    }
    nz = amos::besk(z99, v, kode, n, cy99, &ierr);
    cy.real(std::real(cy99[0]));
    cy.imag(std::imag(cy99[0]));
    do_sferr("kve:", &cy, nz, ierr);
    if (ierr == 2) {
        if (z.real() >= 0 && z.imag() == 0) {
            /* overflow */
            cy.real(INFINITY);
            cy.imag(0);
        }
    }

    return cy;
}

inline double cbesk_wrap_real(double v, double z) {
    std::complex<double> cy, w;
    if (z < 0) {
        return NAN;
    } else if (z == 0) {
        return INFINITY;
    } else if (z > 710 * (1 + fabs(v))) {
        /* Underflow. See uniform expansion https://dlmf.nist.gov/10.41
         * This condition is not a strict bound (it can underflow earlier),
         * rather, we are here working around a restriction in AMOS.
         */
        return 0;
    } else {
        w.real(z);
        w.imag(0);
        cy = cbesk_wrap(v, w);
        return cy.real();
    }
}

inline double cbesk_wrap_real_int(int n, double z) { return cbesk_wrap_real(n, z); }

inline float cbesk_wrap_real(float v, float z) {
    return cbesk_wrap_real(static_cast<double>(v), static_cast<double>(z));
}

inline std::complex<float> cbesk_wrap(float v, std::complex<float> z) {
    std::complex<double> res = cbesk_wrap(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline double cbesk_wrap_e_real(double v, double z) {
    std::complex<double> cy, w;
    if (z < 0) {
        return NAN;
    } else if (z == 0) {
        return INFINITY;
    } else {
        w.real(z);
        w.imag(0);
        cy = cbesk_wrap_e(v, w);
        return cy.real();
    }
}

inline float cbesk_wrap_e_real(float v, float z) {
    return cbesk_wrap_e_real(static_cast<double>(v), static_cast<double>(z));
}

inline std::complex<float> cbesk_wrap_e(float v, std::complex<float> z) {
    std::complex<double> res = cbesk_wrap_e(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline std::complex<double> cbesh_wrap1(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int m = 1;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy;

    std::complex<double> z99 = z;
    std::complex<double> cy99[1] = {NAN};

    cy.real(NAN);
    cy.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besh(z99, v, kode, m, n, cy99, &ierr);
    cy.real(std::real(cy99[0]));
    cy.imag(std::imag(cy99[0]));
    do_sferr("hankel1:", &cy, nz, ierr);
    if (sign == -1) {
        cy = rotate(cy, v);
    }
    return cy;
}

inline std::complex<double> cbesh_wrap1_e(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int m = 1;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy;

    std::complex<double> z99 = z;
    std::complex<double> cy99[1] = {NAN};

    cy.real(NAN);
    cy.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besh(z99, v, kode, m, n, cy99, &ierr);
    cy.real(std::real(cy99[0]));
    cy.imag(std::imag(cy99[0]));
    do_sferr("hankel1e:", &cy, nz, ierr);
    if (sign == -1) {
        cy = rotate(cy, v);
    }
    return cy;
}

inline std::complex<float> cbesh_wrap1(float v, std::complex<float> z) {
    std::complex<double> res = cbesh_wrap1(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline std::complex<float> cbesh_wrap1_e(float v, std::complex<float> z) {
    std::complex<double> res = cbesh_wrap1_e(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline std::complex<double> cbesh_wrap2(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int m = 2;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy;

    std::complex<double> z99 = z;
    std::complex<double> cy99[1] = {NAN};

    cy.real(NAN);
    cy.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besh(z99, v, kode, m, n, cy99, &ierr);
    cy.real(std::real(cy99[0]));
    cy.imag(std::imag(cy99[0]));
    do_sferr("hankel2:", &cy, nz, ierr);
    if (sign == -1) {
        cy = rotate(cy, -v);
    }
    return cy;
}

inline std::complex<double> cbesh_wrap2_e(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int m = 2;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy;

    std::complex<double> z99 = z;
    std::complex<double> cy99[1] = {NAN};

    cy.real(NAN);
    cy.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besh(z99, v, kode, m, n, cy99, &ierr);
    cy.real(cy99[0].real());
    cy.imag(cy99[0].imag());
    do_sferr("hankel2e:", &cy, nz, ierr);
    if (sign == -1) {
        cy = rotate(cy, -v);
    }
    return cy;
}

inline std::complex<float> cbesh_wrap2(float v, std::complex<float> z) {
    std::complex<double> res = cbesh_wrap2(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

inline std::complex<float> cbesh_wrap2_e(float v, std::complex<float> z) {
    std::complex<double> res = cbesh_wrap2_e(
        static_cast<double>(v), std::complex(static_cast<double>(std::real(z)), static_cast<double>(std::imag(z))));
    return {static_cast<float>(std::real(res)), static_cast<float>(std::imag(res))};
}

} // namespace special
