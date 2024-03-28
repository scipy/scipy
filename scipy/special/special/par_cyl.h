#pragma once

#include "specfun/specfun.h"

namespace special {

/*
 * If x > 0 return w1f and w1d. Otherwise set x = abs(x) and return
 * w2f and -w2d.
 */
template <typename T>
void pbwa(T a, T x, T *wf, T *wd) {
    int flag = 0;
    T w1f = 0.0, w1d = 0.0, w2f = 0.0, w2d = 0.0;

    if (x < -5 || x > 5 || a < -5 || a > 5) {
        /*
         * The Zhang and Jin implementation only uses Taylor series;
         * return NaN outside of the range which they are accurate.
         */
        *wf = std::numeric_limits<T>::quiet_NaN();
        *wd = std::numeric_limits<T>::quiet_NaN();
        set_error("pbwa", SF_ERROR_LOSS, NULL);
    } else {
        if (x < 0) {
            x = -x;
            flag = 1;
        }
        specfun::pbwa(a, x, &w1f, &w1d, &w2f, &w2d);
        if (flag) {
            *wf = w2f;
            *wd = -w2d;
        } else {
            *wf = w1f;
            *wd = w1d;
        }
    }
}

template <typename T>
void pbdv(T v, T x, T *pdf, T *pdd) {
    T *dv;
    T *dp;
    int num;

    if (isnan(v) || isnan(x)) {
        *pdf = std::numeric_limits<T>::quiet_NaN();
        *pdd = std::numeric_limits<T>::quiet_NaN();
    } else {
        /* NB. Indexing of DV/DP in specfun.f:PBDV starts from 0, hence +2 */
        num = std::abs((int) v) + 2;
        dv = (T *) malloc(sizeof(T) * 2 * num);
        if (dv == NULL) {
            set_error("pbdv", SF_ERROR_OTHER, "memory allocation error");
            *pdf = std::numeric_limits<T>::quiet_NaN();
            *pdd = std::numeric_limits<T>::quiet_NaN();
        } else {
            dp = dv + num;
            specfun::pbdv(x, v, dv, dp, pdf, pdd);
            free(dv);
        }
    }
}

template <typename T>
void pbvv(T v, T x, T *pvf, T *pvd) {
    T *vv;
    T *vp;
    int num;

    if (isnan(v) || isnan(x)) {
        *pvf = std::numeric_limits<T>::quiet_NaN();
        *pvd = std::numeric_limits<T>::quiet_NaN();
    } else {
        /* NB. Indexing of DV/DP in specfun.f:PBVV starts from 0, hence +2 */
        num = std::abs((int) v) + 2;
        vv = (T *) malloc(sizeof(T) * 2 * num);
        if (vv == NULL) {
            set_error("pbvv", SF_ERROR_OTHER, "memory allocation error");
            *pvf = std::numeric_limits<T>::quiet_NaN();
            *pvd = std::numeric_limits<T>::quiet_NaN();
        } else {
            vp = vv + num;
            specfun::pbvv(x, v, vv, vp, pvf, pvd);
            free(vv);
        }
    }
}

} // namespace special
