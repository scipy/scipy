#pragma once

#include "amos/amos.h"
#include "config.h"
#include "error.h"

extern "C" int cephes_airy(double x, double *ai, double *aip, double *bi, double *bip);

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

template <typename T>
void airy(std::complex<T> z, std::complex<T> &ai, std::complex<T> &aip, std::complex<T> &bi, std::complex<T> &bip) {
    int id = 0;
    int ierr = 0;
    int kode = 1;
    int nz;

    ai.real(NAN);
    ai.imag(NAN);
    bi.real(NAN);
    bi.imag(NAN);
    aip.real(NAN);
    aip.imag(NAN);
    bip.real(NAN);
    bip.imag(NAN);

    ai = amos::airy(z, id, kode, &nz, &ierr);
    do_sferr("airy:", &ai, nz, ierr);

    nz = 0;
    bi = amos::biry(z, id, kode, &ierr);
    do_sferr("airy:", &bi, nz, ierr);

    id = 1;
    aip = amos::airy(z, id, kode, &nz, &ierr);
    do_sferr("airy:", &aip, nz, ierr);

    nz = 0;
    bip = amos::biry(z, id, kode, &ierr);
    do_sferr("airy:", &bip, nz, ierr);
}

template <typename T>
void airye(std::complex<T> z, std::complex<T> &ai, std::complex<T> &aip, std::complex<T> &bi, std::complex<T> &bip) {
    int id = 0;
    int kode = 2; /* Exponential scaling */
    int nz, ierr;

    ai.real(NAN);
    ai.imag(NAN);
    bi.real(NAN);
    bi.imag(NAN);
    aip.real(NAN);
    aip.imag(NAN);
    bip.real(NAN);
    bip.imag(NAN);

    ai = amos::airy(z, id, kode, &nz, &ierr);
    do_sferr("airye:", &ai, nz, ierr);

    nz = 0;
    bi = amos::biry(z, id, kode, &ierr);
    do_sferr("airye:", &bi, nz, ierr);

    id = 1;
    aip = amos::airy(z, id, kode, &nz, &ierr);
    do_sferr("airye:", &aip, nz, ierr);

    nz = 0;
    bip = amos::biry(z, id, kode, &ierr);
    do_sferr("airye:", &bip, nz, ierr);
}

template <typename T>
void airye(T z, T &ai, T &aip, T &bi, T &bip) {
    int id = 0;
    int kode = 2; /* Exponential scaling */
    int nz, ierr;
    std::complex<T> cai, caip, cbi, cbip;

    cai.real(NAN);
    cai.imag(NAN);
    cbi.real(NAN);
    cbi.imag(NAN);
    caip.real(NAN);
    caip.imag(NAN);
    cbip.real(NAN);
    cbip.real(NAN);

    if (z < 0) {
        ai = NAN;
    } else {
        cai = amos::airy(z, id, kode, &nz, &ierr);
        do_sferr("airye:", &cai, nz, ierr);
        ai = cai.real();
    }

    nz = 0;
    cbi = amos::biry(z, id, kode, &ierr);
    do_sferr("airye:", &cbi, nz, ierr);
    bi = std::real(cbi);

    id = 1;
    if (z < 0) {
        aip = NAN;
    } else {
        caip = amos::airy(z, id, kode, &nz, &ierr);
        do_sferr("airye:", &caip, nz, ierr);
        aip = std::real(caip);
    }

    nz = 0;
    cbip = amos::biry(z, id, kode, &ierr);
    do_sferr("airye:", &cbip, nz, ierr);
    bip = cbip.real();
}

template <typename T>
void airy(T x, T &ai, T &aip, T &bi, T &bip) {
    /* For small arguments, use Cephes as it's slightly faster.
     * For large arguments, use AMOS as it's more accurate.
     */
    if (x < -10 || x > 10) {
        std::complex<T> zai, zaip, zbi, zbip;
        airy(std::complex(x), zai, zaip, zbi, zbip);
        ai = std::real(zai);
        aip = std::real(zaip);
        bi = std::real(zbi);
        bip = std::real(zbip);
    } else {
        cephes_airy(x, &ai, &aip, &bi, &bip);
    }
}

} // namespace special
