#pragma once

#include "specfun.h"

namespace xsf {

template <typename T>
T prolate_segv(T m, T n, T c) {
    int kd = 1;
    int int_m, int_n;
    T cv = 0.0, *eg;

    if ((m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (T *) malloc(sizeof(T) * (n - m + 2));
    if (eg == NULL) {
        set_error("pro_cv", SF_ERROR_MEMORY, "memory allocation error");
        return std::numeric_limits<T>::quiet_NaN();
    }
    specfun::Status status = specfun::segv(int_m, int_n, c, kd, &cv, eg);
    free(eg);
    if (status == specfun::Status::NoMemory) {
        set_error("pro_cv", SF_ERROR_MEMORY, "memory allocation error");
        return std::numeric_limits<T>::quiet_NaN();
    }
    return cv;
}

template <typename T>
T oblate_segv(T m, T n, T c) {
    int kd = -1;
    int int_m, int_n;
    T cv = 0.0, *eg;

    if ((m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (T *) malloc(sizeof(T) * (n - m + 2));
    if (eg == NULL) {
        set_error("obl_cv", SF_ERROR_MEMORY, "memory allocation error");
        return std::numeric_limits<T>::quiet_NaN();
    }
    specfun::Status status = specfun::segv(int_m, int_n, c, kd, &cv, eg);
    free(eg);
    if (status == specfun::Status::NoMemory) {
        set_error("obl_cv", SF_ERROR_MEMORY, "memory allocation error");
        return std::numeric_limits<T>::quiet_NaN();
    }
    return cv;
}

template <typename T>
void prolate_aswfa_nocv(T m, T n, T c, T x, T &s1f, T &s1d) {
    int kd = 1;
    int int_m, int_n;
    T cv = 0.0, *eg;

    if ((x >= 1) || (x <= -1) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("pro_ang1", SF_ERROR_DOMAIN, NULL);
        s1d = std::numeric_limits<T>::quiet_NaN();
        s1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (T *) malloc(sizeof(T) * (n - m + 2));
    if (eg == NULL) {
        set_error("pro_ang1", SF_ERROR_MEMORY, "memory allocation error");
        s1d = std::numeric_limits<T>::quiet_NaN();
        s1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    specfun::Status status = specfun::segv(int_m, int_n, c, kd, &cv, eg);
    free(eg);
    if (status == specfun::Status::NoMemory) {
        set_error("pro_ang1", SF_ERROR_MEMORY, "memory allocation error");
        s1d = std::numeric_limits<T>::quiet_NaN();
        s1f = std::numeric_limits<T>::quiet_NaN();
        return;        
    }
    if (specfun::aswfa(x, int_m, int_n, c, kd, cv, &s1f, &s1d)
            == specfun::Status::NoMemory) {
        set_error("prol_ang1", SF_ERROR_MEMORY, "memory allocation error");
        s1d = std::numeric_limits<T>::quiet_NaN();
        s1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
}

template <typename T>
void oblate_aswfa_nocv(T m, T n, T c, T x, T &s1f, T &s1d) {
    int kd = -1;
    int int_m, int_n;
    T cv = 0.0, *eg;

    if ((x >= 1) || (x <= -1) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("obl_ang1", SF_ERROR_DOMAIN, NULL);
        s1d = std::numeric_limits<T>::quiet_NaN();
        s1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (T *) malloc(sizeof(T) * (n - m + 2));
    if (eg == NULL) {
        set_error("obl_ang1", SF_ERROR_MEMORY, "memory allocation error");
        s1d = std::numeric_limits<T>::quiet_NaN();
        s1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    specfun::Status status = specfun::segv(int_m, int_n, c, kd, &cv, eg);
    free(eg);
    if (status == specfun::Status::NoMemory) {
        set_error("obl_ang1", SF_ERROR_MEMORY, "memory allocation error");
        s1d = std::numeric_limits<T>::quiet_NaN();
        s1f = std::numeric_limits<T>::quiet_NaN();
        return;        
    }
    if (specfun::aswfa(x, int_m, int_n, c, kd, cv, &s1f, &s1d)
            == specfun::Status::NoMemory) {
        set_error("obl_ang1", SF_ERROR_MEMORY, "memory allocation error");
        s1d = std::numeric_limits<T>::quiet_NaN();
        s1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
}

template <typename T>
void prolate_aswfa(T m, T n, T c, T cv, T x, T &s1f, T &s1d) {
    if ((x >= 1) || (x <= -1) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("pro_ang1_cv", SF_ERROR_DOMAIN, NULL);
        s1f = std::numeric_limits<T>::quiet_NaN();
        s1d = std::numeric_limits<T>::quiet_NaN();
    } else {
        specfun::Status status = specfun::aswfa(x, static_cast<int>(m),
                                                static_cast<int>(n), c, 1, cv,
                                                &s1f, &s1d);
        if (status == specfun::Status::NoMemory) {
            set_error("pro_ang1_cv", SF_ERROR_MEMORY, "memory allocation error");
            s1d = std::numeric_limits<T>::quiet_NaN();
            s1f = std::numeric_limits<T>::quiet_NaN();
            return;
        }
    }
}

template <typename T>
void oblate_aswfa(T m, T n, T c, T cv, T x, T &s1f, T &s1d) {
    if ((x >= 1) || (x <= -1) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("obl_ang1_cv", SF_ERROR_DOMAIN, NULL);
        s1f = std::numeric_limits<T>::quiet_NaN();
        s1d = std::numeric_limits<T>::quiet_NaN();
    } else {
        specfun::Status status = specfun::aswfa(x, static_cast<int>(m),
                                                static_cast<int>(n), c, -1, cv,
                                                &s1f, &s1d);
        if (status == specfun::Status::NoMemory) {
            set_error("obl_ang1_cv", SF_ERROR_MEMORY, "memory allocation error");
            s1d = std::numeric_limits<T>::quiet_NaN();
            s1f = std::numeric_limits<T>::quiet_NaN();
            return;
        }
    }
}

template <typename T>
void prolate_radial1_nocv(T m, T n, T c, T x, T &r1f, T &r1d) {
    int kf = 1, kd = 1;
    T r2f = 0.0, r2d = 0.0, cv = 0.0, *eg;
    int int_m, int_n;

    if ((x <= 1.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("pro_rad1", SF_ERROR_DOMAIN, NULL);
        r1d = std::numeric_limits<T>::quiet_NaN();
        r1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (T *) malloc(sizeof(T) * (n - m + 2));
    if (eg == NULL) {
        set_error("pro_rad1", SF_ERROR_MEMORY, "memory allocation error");
        r1d = std::numeric_limits<T>::quiet_NaN();
        r1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    specfun::Status status = specfun::segv(int_m, int_n, c, kd, &cv, eg);
    free(eg);
    if (status == specfun::Status::NoMemory) {
        set_error("pro_rad1", SF_ERROR_MEMORY, "memory allocation error");
        r1d = std::numeric_limits<T>::quiet_NaN();
        r1f = std::numeric_limits<T>::quiet_NaN();
        return;        
    }
    if (specfun::rswfp(int_m, int_n, c, x, cv, kf, &r1f, &r1d, &r2f, &r2d)
            == specfun::Status::NoMemory) {
        set_error("pro_rad1", SF_ERROR_MEMORY, "memory allocation error");
        r1d = std::numeric_limits<T>::quiet_NaN();
        r1f = std::numeric_limits<T>::quiet_NaN();
    }
}

template <typename T>
void prolate_radial2_nocv(T m, T n, T c, T x, T &r2f, T &r2d) {
    int kf = 2, kd = 1;
    T r1f = 0.0, r1d = 0.0, cv = 0.0, *eg;
    int int_m, int_n;

    if ((x <= 1.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("pro_rad2", SF_ERROR_DOMAIN, NULL);
        r2d = std::numeric_limits<T>::quiet_NaN();
        r2f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (T *) malloc(sizeof(T) * (n - m + 2));
    if (eg == NULL) {
        set_error("pro_rad2", SF_ERROR_MEMORY, "memory allocation error");
        r2d = std::numeric_limits<T>::quiet_NaN();
        r2f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    specfun::Status status = specfun::segv(int_m, int_n, c, kd, &cv, eg);
    free(eg);
    if (status == specfun::Status::NoMemory) {
        set_error("pro_rad2", SF_ERROR_MEMORY, "memory allocation error");
        r2d = std::numeric_limits<T>::quiet_NaN();
        r2f = std::numeric_limits<T>::quiet_NaN();
        return;        
    }
    if (specfun::rswfp(int_m, int_n, c, x, cv, kf, &r1f, &r1d, &r2f, &r2d)
            == specfun::Status::NoMemory) {
        set_error("pro_rad2", SF_ERROR_MEMORY, "memory allocation error");
        r2d = std::numeric_limits<T>::quiet_NaN();
        r2f = std::numeric_limits<T>::quiet_NaN();
    }
}

template <typename T>
void prolate_radial1(T m, T n, T c, T cv, T x, T &r1f, T &r1d) {
    int kf = 1;
    T r2f = 0.0, r2d = 0.0;
    int int_m, int_n;

    if ((x <= 1.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("pro_rad1_cv", SF_ERROR_DOMAIN, NULL);
        r1f = std::numeric_limits<T>::quiet_NaN();
        r1d = std::numeric_limits<T>::quiet_NaN();
    } else {
        int_m = (int) m;
        int_n = (int) n;
        if (specfun::rswfp(int_m, int_n, c, x, cv, kf, &r1f, &r1d, &r2f, &r2d)
                == specfun::Status::NoMemory) {
            set_error("pro_rad1_cv", SF_ERROR_MEMORY, NULL);
            r1f = std::numeric_limits<T>::quiet_NaN();
            r1d = std::numeric_limits<T>::quiet_NaN();
        }
    }
}

template <typename T>
void prolate_radial2(T m, T n, T c, T cv, T x, T &r2f, T &r2d) {
    int kf = 2;
    T r1f = 0.0, r1d = 0.0;
    int int_m, int_n;

    if ((x <= 1.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("pro_rad2_cv", SF_ERROR_DOMAIN, NULL);
        r2f = std::numeric_limits<double>::quiet_NaN();
        r2d = std::numeric_limits<double>::quiet_NaN();
    } else {
        int_m = (int) m;
        int_n = (int) n;
        if (specfun::rswfp(int_m, int_n, c, x, cv, kf, &r1f, &r1d, &r2f, &r2d)
                == specfun::Status::NoMemory) {
            set_error("pro_rad2_cv", SF_ERROR_MEMORY, NULL);
            r2f = std::numeric_limits<double>::quiet_NaN();
            r2d = std::numeric_limits<double>::quiet_NaN();
        }
    }
}

template <typename T>
void oblate_radial1_nocv(T m, T n, T c, T x, T &r1f, T &r1d) {
    int kf = 1, kd = -1;
    T r2f = 0.0, r2d = 0.0, cv = 0.0, *eg;
    int int_m, int_n;

    if ((x < 0.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("obl_rad1", SF_ERROR_DOMAIN, NULL);
        r1d = std::numeric_limits<T>::quiet_NaN();
        r1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (T *) malloc(sizeof(T) * (n - m + 2));
    if (eg == NULL) {
        set_error("obl_rad1", SF_ERROR_MEMORY, "memory allocation error");
        r1d = std::numeric_limits<T>::quiet_NaN();
        r1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    specfun::Status status = specfun::segv(int_m, int_n, c, kd, &cv, eg);
    free(eg);
    if (status == specfun::Status::NoMemory) {
        set_error("obl_rad1", SF_ERROR_MEMORY, "memory allocation error");
        r1d = std::numeric_limits<T>::quiet_NaN();
        r1f = std::numeric_limits<T>::quiet_NaN();
        return;        
    }
    if (specfun::rswfo(int_m, int_n, c, x, cv, kf, &r1f, &r1d, &r2f, &r2d)
            == specfun::Status::NoMemory) {
        set_error("obl_rad1", SF_ERROR_MEMORY, "memory allocation error");
        r1d = std::numeric_limits<T>::quiet_NaN();
        r1f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
}

template <typename T>
void oblate_radial2_nocv(T m, T n, T c, T x, T &r2f, T &r2d) {
    int kf = 2, kd = -1;
    T r1f = 0.0, r1d = 0.0, cv = 0.0, *eg;
    int int_m, int_n;

    if ((x < 0.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("obl_rad2", SF_ERROR_DOMAIN, NULL);
        r2d = std::numeric_limits<T>::quiet_NaN();
        r2f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (T *) malloc(sizeof(T) * (n - m + 2));
    if (eg == NULL) {
        set_error("obl_rad2", SF_ERROR_MEMORY, "memory allocation error");
        r2d = std::numeric_limits<T>::quiet_NaN();
        r2f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
    specfun::Status status = specfun::segv(int_m, int_n, c, kd, &cv, eg);
    free(eg);
    if (status == specfun::Status::NoMemory) {
        set_error("obl_rad2", SF_ERROR_MEMORY, "memory allocation error");
        r2d = std::numeric_limits<T>::quiet_NaN();
        r2f = std::numeric_limits<T>::quiet_NaN();
        return;        
    }
    if (specfun::rswfo(int_m, int_n, c, x, cv, kf, &r1f, &r1d, &r2f, &r2d)
            == specfun::Status::NoMemory){
        set_error("obl_rad2", SF_ERROR_MEMORY, "memory allocation error");
        r2d = std::numeric_limits<T>::quiet_NaN();
        r2f = std::numeric_limits<T>::quiet_NaN();
        return;
    }
}

template <typename T>
void oblate_radial1(T m, T n, T c, T cv, T x, T &r1f, T &r1d) {
    int kf = 1;
    T r2f = 0.0, r2d = 0.0;

    if ((x < 0.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("obl_rad1_cv", SF_ERROR_DOMAIN, NULL);
        r1f = std::numeric_limits<T>::quiet_NaN();
        r1d = std::numeric_limits<T>::quiet_NaN();
    } else {
        if (specfun::rswfo(static_cast<int>(m), static_cast<int>(n),
                           c, x, cv, kf, &r1f, &r1d, &r2f, &r2d)
                == specfun::Status::NoMemory) {
            set_error("obl_rad1_cv", SF_ERROR_MEMORY, "memory allocation error");
            r1d = std::numeric_limits<T>::quiet_NaN();
            r1f = std::numeric_limits<T>::quiet_NaN();
            return;     
        }
    }
}

template <typename T>
void oblate_radial2(T m, T n, T c, T cv, T x, T &r2f, T &r2d) {
    int kf = 2;
    T r1f = 0.0, r1d = 0.0;

    if ((x < 0.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("obl_rad2_cv", SF_ERROR_DOMAIN, NULL);
        r2f = std::numeric_limits<T>::quiet_NaN();
        r2d = std::numeric_limits<T>::quiet_NaN();
    } else {
        if (specfun::rswfo(static_cast<int>(m), static_cast<int>(n),
                           c, x, cv, kf, &r1f, &r1d, &r2f, &r2d)
                == specfun::Status::NoMemory) {
            set_error("obl_rad2_cv", SF_ERROR_MEMORY, "memory allocation error");
            r2d = std::numeric_limits<T>::quiet_NaN();
            r2f = std::numeric_limits<T>::quiet_NaN();
            return; 
        }
    }
}

} // namespace xsf
