/*
 * mconf configures NANS, INFINITYs etc. for cephes and includes some standard
 * headers. Although erfinv and erfcinv are not defined in cephes, erf and erfc
 * are. We want to keep the behaviour consistent for the inverse functions and
 * so need to include mconf.
 */
#include "mconf.h"

/*
 * Inverse of the error function.
 *
 * Computes the inverse of the error function on the restricted domain
 * -1 < y < 1. This restriction ensures the existence of a unique result
 * such that erf(erfinv(y)) = y.
 */
double erfinv(double y) {
    const double domain_lb = -1;
    const double domain_ub = 1;

    const double thresh = 1e-7;

    /* 
     * For small arguments, use the Taylor expansion
     * erf(y) = 2/\sqrt{\pi} (y - y^3 / 3 + O(y^5)),    y\to 0 
     * where we only retain the linear term.
     * Otherwise, y + 1 loses precision for |y| << 1.
     */
    if ((-thresh < y) && (y < thresh)){
        return y / M_2_SQRTPI;
    } 
    if ((domain_lb < y) && (y < domain_ub)) {
        return ndtri(0.5 * (y+1)) * NPY_SQRT1_2;
    }
    else if (y == domain_lb) {
        return -NPY_INFINITY;
    }
    else if (y == domain_ub) {
        return NPY_INFINITY;
    }
    else if (cephes_isnan(y)) {
        sf_error("erfinv", SF_ERROR_DOMAIN, NULL);
        return y;
    }
    else {
        sf_error("erfinv", SF_ERROR_DOMAIN, NULL);
        return NPY_NAN;
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
    const double domain_lb = 0;
    const double domain_ub = 2;

    if ((domain_lb < y) && (y < domain_ub)) {
        return -ndtri(0.5 * y) * NPY_SQRT1_2;
    }
    else if (y == domain_lb) {
        return NPY_INFINITY;
    }
    else if (y == domain_ub) {
        return -NPY_INFINITY;
    }
    else if (cephes_isnan(y)) {
        sf_error("erfcinv", SF_ERROR_DOMAIN, NULL);
        return y;
    }
    else {
        sf_error("erfcinv", SF_ERROR_DOMAIN, NULL);
        return NPY_NAN;
    }
}
