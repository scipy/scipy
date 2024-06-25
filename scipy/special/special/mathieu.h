#pragma once

#include "specfun/specfun.h"

namespace special {

template <typename T>
T sem_cva(T m, T q);

template <typename T>
void sem(T m, T q, T x, T &csf, T &csd);

/* Mathieu functions */
/* Characteristic values */
template <typename T>
T cem_cva(T m, T q) {
    int int_m, kd = 1;

    if ((m < 0) || (m != floor(m))) {
        set_error("cem_cva", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<T>::quiet_NaN();
    }
    int_m = (int) m;
    if (q < 0) {
        /* https://dlmf.nist.gov/28.2#E26 */
        if (int_m % 2 == 0) {
            return cem_cva(m, -q);
        } else {
            return sem_cva(m, -q);
        }
    }

    if (int_m % 2) {
        kd = 2;
    }
    return specfun::cva2(kd, int_m, q);
}

template <typename T>
T sem_cva(T m, T q) {
    int int_m, kd = 4;

    if ((m <= 0) || (m != floor(m))) {
        set_error("cem_cva", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<T>::quiet_NaN();
    }
    int_m = (int) m;
    if (q < 0) {
        /* https://dlmf.nist.gov/28.2#E26 */
        if (int_m % 2 == 0) {
            return sem_cva(m, -q);
        } else {
            return cem_cva(m, -q);
        }
    }
    if (int_m % 2) {
        kd = 3;
    }
    return specfun::cva2(kd, int_m, q);
}

/* Mathieu functions */
template <typename T>
void cem(T m, T q, T x, T &csf, T &csd) {
    int int_m, kf = 1, sgn;
    T f = 0.0, d = 0.0;
    if ((m < 0) || (m != floor(m))) {
        csf = std::numeric_limits<T>::quiet_NaN();
        csd = std::numeric_limits<T>::quiet_NaN();
        set_error("cem", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        if (q < 0) {
            /* https://dlmf.nist.gov/28.2#E34 */
            if (int_m % 2 == 0) {
                sgn = ((int_m / 2) % 2 == 0) ? 1 : -1;
                cem(m, -q, 90 - x, f, d);
                csf = sgn * f;
                csd = -sgn * d;

            } else {
                sgn = ((int_m / 2) % 2 == 0) ? 1 : -1;
                sem(m, -q, 90 - x, f, d);
                csf = sgn * f;
                csd = -sgn * d;
            }
        } else {
            specfun::mtu0(kf, int_m, q, x, &csf, &csd);
        }
    }
}

template <typename T>
void sem(T m, T q, T x, T &csf, T &csd) {
    int int_m, kf = 2, sgn;
    T f = 0.0, d = 0.0;
    if ((m < 0) || (m != floor(m))) {
        csf = std::numeric_limits<T>::quiet_NaN();
        csd = std::numeric_limits<T>::quiet_NaN();
        set_error("sem", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        if (int_m == 0) {
            csf = 0;
            csd = 0;
        } else if (q < 0) {
            /* https://dlmf.nist.gov/28.2#E34 */
            if (int_m % 2 == 0) {
                sgn = ((int_m / 2) % 2 == 0) ? -1 : 1;
                sem(m, -q, 90 - x, f, d);
                csf = sgn * f;
                csd = -sgn * d;
            } else {
                sgn = ((int_m / 2) % 2 == 0) ? 1 : -1;
                cem(m, -q, 90 - x, f, d);
                csf = sgn * f;
                csd = -sgn * d;
            }
        } else {
            specfun::mtu0(kf, int_m, q, x, &csf, &csd);
        }
    }
}

template <typename T>
void mcm1(T m, T q, T x, T &f1r, T &d1r) {
    int int_m, kf = 1, kc = 1;
    T f2r = 0.0, d2r = 0.0;

    if ((m < 0) || (m != floor(m)) || (q < 0)) {
        f1r = std::numeric_limits<T>::quiet_NaN();
        d1r = std::numeric_limits<T>::quiet_NaN();
        set_error("mcm1", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        specfun::mtu12(kf, kc, int_m, q, x, &f1r, &d1r, &f2r, &d2r);
    }
}

template <typename T>
void msm1(T m, T q, T x, T &f1r, T &d1r) {
    int int_m, kf = 2, kc = 1;
    T f2r = 0.0, d2r = 0.0;

    if ((m < 1) || (m != floor(m)) || (q < 0)) {
        f1r = std::numeric_limits<T>::quiet_NaN();
        d1r = std::numeric_limits<T>::quiet_NaN();
        set_error("msm1", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        specfun::mtu12(kf, kc, int_m, q, x, &f1r, &d1r, &f2r, &d2r);
    }
}

template <typename T>
void mcm2(T m, T q, T x, T &f2r, T &d2r) {
    int int_m, kf = 1, kc = 2;
    T f1r = 0.0, d1r = 0.0;

    if ((m < 0) || (m != floor(m)) || (q < 0)) {
        f2r = std::numeric_limits<T>::quiet_NaN();
        d2r = std::numeric_limits<T>::quiet_NaN();
        set_error("mcm2", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        specfun::mtu12(kf, kc, int_m, q, x, &f1r, &d1r, &f2r, &d2r);
    }
}

template <typename T>
void msm2(T m, T q, T x, T &f2r, T &d2r) {
    int int_m, kf = 2, kc = 2;
    T f1r = 0.0, d1r = 0.0;

    if ((m < 1) || (m != floor(m)) || (q < 0)) {
        f2r = std::numeric_limits<T>::quiet_NaN();
        d2r = std::numeric_limits<T>::quiet_NaN();
        set_error("msm2", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        specfun::mtu12(kf, kc, int_m, q, x, &f1r, &d1r, &f2r, &d2r);
    }
}

} // namespace special
