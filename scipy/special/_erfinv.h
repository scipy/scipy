#ifndef ERF_INV_H
#define ERF_INV_H

/*
 * mconf configures NANS, INFINITYs etc. for cephes and includes some standard
 * headers. Although erfinv and erfcinv are not defined in cephes, erf and erfc
 * are. We want to keep the behaviour consistent for the inverse functions and
 * so need to include mconf.
 */
#include "cephes/mconf.h"

/*
 * Inverse of the error function.
 *
 * Computes the inverse of the error function on the restricted domain
 * -1 < y < 1. This restriction ensures the existence of a unique result
 * such that erf(erfinv(y)) = y.
 */
double erfinv(double y) {
    if (cephes_isnan(y)) {
        sf_error("erfinv", SF_ERROR_DOMAIN, NULL);
        return NPY_NAN;
    }
    else if (y <= -1) {
        sf_error("erfinv", SF_ERROR_DOMAIN, NULL);
        return -NPY_INFINITY;
    }
    else if (y >= 1) {
        sf_error("erfinv", SF_ERROR_DOMAIN, NULL);
        return NPY_INFINITY;
    }
    else {
        return ndtri(0.5 * (y+1)) * NPY_SQRT1_2;
    }
}

/*
 * Inverse of the complementary error function.
 *
 * Computes the inverse of the complimentary error function on the restricted
 * domain 0 < y < 2. This restriction ensures the existence of a unique result
 * such that erfc(erfcinv(y)) = y.
 */
double erfcinv(double y) {
    if (cephes_isnan(y)) {
        sf_error("erfcinv", SF_ERROR_DOMAIN, NULL);
        return NPY_NAN;
    }
    else if (y <= 0) {
        sf_error("erfcinv", SF_ERROR_DOMAIN, NULL);
        return NPY_INFINITY;
    }
    else if (y >= 2) {
        sf_error("erfcinv", SF_ERROR_DOMAIN, NULL);
        return -NPY_INFINITY;
    }
    else {
        return -ndtri(0.5 * y) * NPY_SQRT1_2;
    }
}

#endif // ERF_INV_H
