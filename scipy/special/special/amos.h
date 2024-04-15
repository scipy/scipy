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

inline sf_error_t ierr_to_sferr(int nz, int ierr) {
    /* Return sf_error equivalents for ierr values */

    if (nz != 0) {
        return SF_ERROR_UNDERFLOW;
    }

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

    return SF_ERROR_OK;
}

template <typename T>
void set_error_and_nan(const char *name, sf_error_t code, T &value) {
    if (code != SF_ERROR_OK) {
        set_error(name, code, nullptr);

        if (code == SF_ERROR_DOMAIN || code == SF_ERROR_OVERFLOW || code == SF_ERROR_NO_RESULT) {
            value = std::numeric_limits<T>::quiet_NaN();
        }
    }
}

template <typename T>
void set_error_and_nan(const char *name, sf_error_t code, std::complex<T> &value) {
    if (code != SF_ERROR_OK) {
        set_error(name, code, nullptr);

        if (code == SF_ERROR_DOMAIN || code == SF_ERROR_OVERFLOW || code == SF_ERROR_NO_RESULT) {
            value.real(std::numeric_limits<T>::quiet_NaN());
            value.imag(std::numeric_limits<T>::quiet_NaN());
        }
    }
}

template <typename T>
void airy(std::complex<T> z, std::complex<T> &ai, std::complex<T> &aip, std::complex<T> &bi, std::complex<T> &bip) {
    int id = 0;
    int ierr = 0;
    int kode = 1;
    int nz;

    ai = amos::airy(z, id, kode, &nz, &ierr);
    set_error_and_nan("airy:", ierr_to_sferr(nz, ierr), ai);

    nz = 0;
    bi = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airy:", ierr_to_sferr(nz, ierr), bi);

    id = 1;
    aip = amos::airy(z, id, kode, &nz, &ierr);
    set_error_and_nan("airy:", ierr_to_sferr(nz, ierr), aip);

    nz = 0;
    bip = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airy:", ierr_to_sferr(nz, ierr), bip);
}

template <typename T>
void airye(std::complex<T> z, std::complex<T> &ai, std::complex<T> &aip, std::complex<T> &bi, std::complex<T> &bip) {
    int id = 0;
    int kode = 2; /* Exponential scaling */
    int nz, ierr;

    ai = amos::airy(z, id, kode, &nz, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), ai);

    nz = 0;
    bi = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), bi);

    id = 1;
    aip = amos::airy(z, id, kode, &nz, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), aip);

    nz = 0;
    bip = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), bip);
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
        set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), cai);
        ai = std::real(cai);
    }

    nz = 0;
    cbi = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), cbi);
    bi = std::real(cbi);

    id = 1;
    if (z < 0) {
        aip = NAN;
    } else {
        caip = amos::airy(z, id, kode, &nz, &ierr);
        set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), caip);
        aip = std::real(caip);
    }

    nz = 0;
    cbip = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), cbip);
    bip = std::real(cbip);
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
