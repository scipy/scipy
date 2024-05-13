#pragma once

#include <iostream>

#include "config.h"
#include "error.h"
#include "gamma.h"

namespace special {

/**
 * Compute the Legendre polynomials of degree n and argument z.
 *
 * @param[in] n degree of the polynomial.
 * @param[in] z the argument, a floating-point value
 * @param[in] callback a function to be called like callback(j, z, p, args...) for j = 0, j = 1, ..., j = n
 * @param[in] args arguments to forward to the callback
 */
template <typename T, typename Callable, typename... Args>
T legendre_p(int n, T z, Callable callback, Args &&...args) {
    T p = 1;
    callback(0, z, p, std::forward<Args>(args)...);

    if (n > 0) {
        T p_prev = p;
        p = z;
        callback(1, z, p, std::forward<Args>(args)...);

        for (int j = 2; j <= n; ++j) {
            T p_prev_prev = p_prev;
            p_prev = p;
            p = (T(2 * j - 1) * z * p_prev - T(j - 1) * p_prev_prev) / T(j);
            callback(j, z, p, std::forward<Args>(args)...);
        }
    }

    return p;
}

template <typename T, size_t N>
struct legendre_p_diff_callback {
    T p_prev[N + 1];
    T p_prev_prev[N + 1];

    template <typename Callable, typename... Args>
    void operator()(int j, T z, T p, Callable callback, Args &&...args) {
        T res[N + 1];

        res[0] = p;

        for (size_t r = 1; r <= N; ++r) {
            if (int(r) > j) {
                res[r] = 0;
            } else {
                res[r] = (T(2 * j - 1) * (T(r) * p_prev[r - 1] + z * p_prev[r]) - T(j - 1) * p_prev_prev[r]) / T(j);
            }
        }

        for (size_t i = 0; i < N + 1; ++i) {
            p_prev_prev[i] = p_prev[i];
            p_prev[i] = res[i];
        }

        callback(j, z, res, std::forward<Args>(args)...);
    }
};

template <typename T>
T legendre_p(int n, T z) {
    return legendre_p(n, z, [](int j, T z, T p) {});
}

template <typename T>
void legendre_p(int n, T z, T &res, T &res_jac) {
    legendre_p(n, z, legendre_p_diff_callback<T, 1>(), [&res, &res_jac](int j, T z, const T(&p)[2]) {
        res = p[0];
        res_jac = p[1];
    });
}

template <typename T>
void legendre_p(int n, T z, T &res, T &res_jac, T &res_hess) {
    legendre_p(n, z, legendre_p_diff_callback<T, 2>(), [&res, &res_jac, &res_hess](int j, T z, const T(&p)[3]) {
        res = p[0];
        res_jac = p[1];
        res_hess = p[2];
    });
}

template <typename T, typename OutputVec>
void legendre_p_all(T z, OutputVec res) {
    int n = res.extent(0) - 1;

    legendre_p(n, z, [res](int j, T z, T p) { res(j) = p; });
}

template <typename T, typename OutputVec1, typename OutputVec2>
void legendre_p_all(T z, OutputVec1 res, OutputVec2 res_jac) {
    int n = res.extent(0) - 1;

    legendre_p(n, z, legendre_p_diff_callback<T, 1>(), [res, res_jac](int j, T z, const T(&p)[2]) {
        res(j) = p[0];
        res_jac(j) = p[1];
    });
}

template <typename T, typename OutputVec1, typename OutputVec2, typename OutputVec3>
void legendre_p_all(T z, OutputVec1 res, OutputVec2 res_jac, OutputVec3 res_hess) {
    int n = res.extent(0) - 1;

    legendre_p(n, z, legendre_p_diff_callback<T, 2>(), [res, res_jac, res_hess](int j, T z, const T(&p)[3]) {
        res(j) = p[0];
        res_jac(j) = p[1];
        res_hess(j) = p[2];
    });
}

// Translated into C++ by SciPy developers in 2024.
// Original comments appear below.
//
// =====================================================
// Purpose: Compute the associated Legendre functions
//          Pmn(x) and their derivatives Pmn'(x) for
//          real argument
// Input :  x  --- Argument of Pmn(x)
//          m  --- Order of Pmn(x),  m = 0,1,2,...,n
//          n  --- Degree of Pmn(x), n = 0,1,2,...,N
//          mm --- Physical dimension of PM and PD
// Output:  PM(m,n) --- Pmn(x)
//          PD(m,n) --- Pmn'(x)
// =====================================================

template <typename T>
T assoc_legendre_p_diag(int m, int type, T x) {
    int m_abs = std::abs(m);

    T res = 1;

    T type_sign;
    if (type == 3) {
        type_sign = -1;
    } else {
        type_sign = 1;
    }

    // need to take care with complex arithmetic, due to signed zeros and sqrt(0j) != sqrt(-0j)
    bool m_odd = m_abs % 2;
    if (m_odd) {
        if (type == 3) {
            res *= std::sqrt(x * x - T(1));
            if (std::real(x) < 0) {
                res = -res;
            }
        } else {
            res *= -std::sqrt(T(1) - x * x);
        }
    }

    // unroll the loop to avoid the sqrt
    for (int j = 1 + m_odd; j <= m_abs; j += 2) {
        res *= T(2 * j - 1) * T(2 * j + 1) * type_sign * (T(1) - x * x);
    }

    if (m < 0) {
        res *= 1 / std::tgamma(2 * m_abs + 1);
        res *= std::pow(-1, m_abs);
        if (m_odd && type == 3) {
            res *= -1;
        }
    }

    return res;
}

template <typename T, typename Callable, typename... Args>
T assoc_legendre_p(int n, int m, int type, T x, Callable callback, Args &&...args) {
    int m_abs = std::abs(m);
    if (m_abs > n) {
        for (int j = 0; j <= n; ++j) {
            callback(j, m, type, x, 0, 0, std::forward<Args>(args)...);
        }

        return 0;
    }

    bool m_odd = m_abs % 2;
    if (m_odd) {
        callback(0, m, type, x, 0, 0, std::forward<Args>(args)...);
    }
    for (int j = 1 + m_odd; j <= m_abs; j += 2) {
        callback(j - 1, m, type, x, 0, 0, std::forward<Args>(args)...);
        callback(j, m, type, x, 0, 0, std::forward<Args>(args)...);
    }

    T p = assoc_legendre_p_diag(m, type, x);
    callback(m_abs, m, type, x, p, 0, std::forward<Args>(args)...);

    if (m_abs != n) {
        T p_prev = p;
        p = T(2 * (m_abs + 1) - 1) * x * p_prev / T(m_abs + 1 - m);
        callback(m_abs + 1, m, type, x, p, p_prev, std::forward<Args>(args)...);

        for (int j = m_abs + 2; j <= n; ++j) {
            T p_prev_prev = p_prev;
            p_prev = p;
            p = (T(2 * j - 1) * x * p_prev - T(m + j - 1) * p_prev_prev) / T(j - m);
            callback(j, m, type, x, p, p_prev, std::forward<Args>(args)...);
        }
    }

    return p;
}

template <typename T>
T assoc_legendre_p_jac_diag(int m, int type, T x) {
    if (m == 0) {
        return 0;
    }

    T type_sign;
    if (type == 3) {
        type_sign = -1;
    } else {
        type_sign = 1;
    }

    if (m == 1) {
        if (type == 3) {
            T res = x / std::sqrt(x * x - T(1));
            if (std::real(x) < 0) {
                res *= -1;
            }

            return res;
        }

        return x / std::sqrt(T(1) - x * x);
    }

    if (m == -1) {
        if (type == 3) {
            T res = x / (T(2) * std::sqrt(x * x - T(1)));
            if (std::real(x) < 0) {
                res *= -1;
            }

            return res;
        }

        return -x / (T(2) * std::sqrt(T(1) - x * x));
    }

    if (m < 0) {
        return x * type_sign * assoc_legendre_p_diag(m + 2, type, x) / T(4 * (m + 1));
    }

    return -T(4 * (m - 2) * m + 3) * T(m) * type_sign * x * assoc_legendre_p_diag(m - 2, type, x);
}

template <typename T>
struct remove_complex {
    using type = T;
};

template <typename T>
struct remove_complex<std::complex<T>> {
    using type = T;
};

template <typename T>
using remove_complex_t = typename remove_complex<T>::type;

template <typename T>
T assoc_legendre_p_jac_next(int n, int m, int type, T z, T p, T p_prev, T p_jac_prev, T p_jac_prev_prev) {
    int m_abs = std::abs(m);
    if (m_abs > n) {
        return 0;
    }

    T type_sign;
    if (type == 3) {
        type_sign = -1;
    } else {
        type_sign = 1;
    }

    if (std::abs(std::real(z)) == 1 && std::imag(z) == 0) {
        if (m == 0) {
            return T(n) * T(n + 1) * std::pow(z, T(n + 1)) / T(2);
        }

        if (m == 1) {
            return std::pow(z, T(n)) * std::numeric_limits<remove_complex_t<T>>::infinity();
        }

        if (m == 2) {
            return -type_sign * T(n + 2) * T(n + 1) * T(n) * T(n - 1) * std::pow(z, T(n + 1)) / T(4);
        }

        if (m == -2) {
            return -type_sign * std::pow(z, T(n + 1)) / T(4);
        }

        if (m == -1) {
            return -std::pow(z, T(n)) * std::numeric_limits<remove_complex_t<T>>::infinity();
        }

        return 0;
    }

    if (m_abs == n) {
        return assoc_legendre_p_jac_diag(m, type, z);
    }

    if (m_abs < n) {
        return (T(2 * n - 1) * (p_prev + z * p_jac_prev) - T(n + m - 1) * p_jac_prev_prev) / T(n - m);
    }

    return 0;
}

template <typename T>
T assoc_legendre_p_hess_diag(int m, int type, T z) {
    if (m == 0) {
        return 0;
    }

    if (m == 1) {
        return 1 / (std::sqrt(1 - z * z) * (1 - z * z));
    }

    if (m == -1) {
        return -1 / (2 * std::sqrt(1 - z * z) * (1 - z * z));
    }

    if (m == 2) {
        return -6;
    }

    if (m == -2) {
        return -1 / T(4);
    }

    if (m == 3) {
        return 45 * (1 - 2 * z * z) / std::sqrt(1 - z * z);
    }

    if (m == -3) {
        return (2 * z * z - 1) / (16 * std::sqrt(1 - z * z));
    }

    if (m < 0) {
        return (T(m + 1) * z * z + 1) * assoc_legendre_p_diag(m + 4, type, z) / T(16 * T(m + 1) * T(m + 2) * T(m + 3));
    }

    return T(2 * m - 1) * T(2 * m - 3) * T(2 * m - 5) * T(2 * m - 7) * m * ((m - 1) * z * z - 1) *
           assoc_legendre_p_diag(m - 4, type, z);
}

template <typename T>
T assoc_legendre_p_hess_next(
    int n, int m, int type, T z, T p, T p_prev, T p_jac, T p_jac_prev, T p_hess_prev, T p_hess_prev_prev
) {
    // need to complete these
    if (std::abs(z) == 1) {
        if (m == 0) {
            return T(n + 2) * T(n + 1) * T(n) * T(n - 1) / T(8);
        }

        if (m == 1) {
            return std::numeric_limits<T>::infinity();
        }

        if (m == 2) {
            return -T((n + 1) * n - 3) * T(n + 2) * T(n + 1) * T(n) * T(n - 1) / T(12);
        }

        if (m == 3) {
            return std::numeric_limits<T>::infinity();
        }

        if (m == 4) {
            return T(n + 4) * T(n + 3) * T(n + 2) * T(n + 1) * T(n) * T(n - 1) * T(n - 2) * T(n - 3) / T(48);
        }

        if (m == -4) {
            return 0;
        }

        if (m == -3) {
            return -std::numeric_limits<T>::infinity();
        }

        if (m == -2) {
            return -T(1) / T(4);
        }

        if (m == -1) {
            return -std::numeric_limits<T>::infinity();
        }

        return 0;
    }

    int m_abs = std::abs(m);
    if (m_abs == n) {
        return assoc_legendre_p_hess_diag(m, type, z);
    }

    if (m_abs < n) {
        return (T(2 * n - 1) * (z * p_hess_prev + T(2) * p_jac_prev) - T(n + m - 1) * p_hess_prev_prev) / T(n - m);
    }

    return 0;
}

template <typename T, size_t N>
struct assoc_legendre_p_diff_callback;

template <typename T>
struct assoc_legendre_p_diff_callback<T, 1> {
    T p_jac_prev;
    T p_jac_prev_prev;

    assoc_legendre_p_diff_callback() : p_jac_prev(0), p_jac_prev_prev(0) {}

    template <typename Callable, typename... Args>
    void operator()(int j, int i, int type, T z, T p, T p_prev, Callable callback, Args &&...args) {
        T p_jac = assoc_legendre_p_jac_next(j, i, type, z, p, p_prev, p_jac_prev, p_jac_prev_prev);
        callback(j, i, type, z, {p, p_jac}, {p_prev, p_jac_prev}, std::forward<Args>(args)...);

        p_jac_prev_prev = p_jac_prev;
        p_jac_prev = p_jac;
    }
};

template <typename T>
struct assoc_legendre_p_diff_callback<T, 2> {
    assoc_legendre_p_diff_callback<T, 1> callback_jac;
    T p_hess_prev;
    T p_hess_prev_prev;

    assoc_legendre_p_diff_callback() : p_hess_prev(0), p_hess_prev_prev(0) {}

    template <typename Callable, typename... Args>
    void operator()(int j, int i, int type, T z, T p, T p_prev, Callable callback, Args &&...args) {
        T p_jac;
        T p_jac_prev;
        callback_jac(
            j, i, type, z, p, p_prev,
            [&p_jac, &p_jac_prev](int j, int i, int type, T z, const T(&p)[2], const T(&p_prev)[2]) {
                p_jac = p[1];
                p_jac_prev = p_prev[1];
            }
        );

        T p_hess =
            assoc_legendre_p_hess_next(j, i, type, z, p, p_prev, p_jac, p_jac_prev, p_hess_prev, p_hess_prev_prev);
        callback(j, i, type, z, {p, p_jac, p_hess}, {p_prev, p_jac_prev, p_hess_prev}, std::forward<Args>(args)...);

        p_hess_prev_prev = p_hess_prev;
        p_hess_prev = p_hess;
    }
};

/*
    if (std::abs(x) == 1) {
        for (long i = 1; i <= n; i++) {
            p(0, i) = std::pow(x, i);
        }
    }
*/

template <typename T>
T assoc_legendre_p(int n, int m, int type, T x) {
    return assoc_legendre_p(n, m, type, x, [](int j, int i, int type, T x, T value, T value_prev) {});
}

template <typename T>
void assoc_legendre_p(int n, int m, int type, T z, T &res, T &res_jac) {
    assoc_legendre_p(
        n, m, type, z, assoc_legendre_p_diff_callback<T, 1>(),
        [&res, &res_jac](int j, int i, int type, T z, const T(&p)[2], const T(&p_prev)[2]) {
            res = p[0];
            res_jac = p[1];
        }
    );
}

template <typename T>
void assoc_legendre_p(int n, int m, int type, T z, T &res, T &res_jac, T &res_hess) {
    assoc_legendre_p(
        n, m, type, z, assoc_legendre_p_diff_callback<T, 2>(),
        [&res, &res_jac, &res_hess](int j, int i, int type, T z, const T(&p)[3], const T(&p_prev)[3]) {
            res = p[0];
            res_jac = p[1];
            res_hess = p[2];
        }
    );
}

template <typename T, typename OutputMat>
void assoc_legendre_p_all(int type, T x, OutputMat res) {
    int m = (res.extent(0) - 1) / 2;
    int n = res.extent(1) - 1;

    for (int i = -m; i <= m; ++i) {
        assoc_legendre_p(n, i, type, x, [res](int j, int i, int type, T x, T p, T p_prev) {
            if (i < 0) {
                int i_abs = std::abs(i);
                res(res.extent(0) - i_abs, j) = p;
            } else {
                res(i, j) = p;
            }
        });
    }
}

template <typename T, typename OutputMat1, typename OutputMat2>
void assoc_legendre_p_all(int type, T x, OutputMat1 res, OutputMat2 res_jac) {
    int m = (res.extent(0) - 1) / 2;
    int n = res.extent(1) - 1;

    for (int i = -m; i <= m; ++i) {
        assoc_legendre_p(
            n, i, type, x, assoc_legendre_p_diff_callback<T, 1>(),
            [res, res_jac](int j, int i, int type, T z, const T(&p)[2], const T(&p_prev)[2]) {
                if (i < 0) {
                    int i_abs = std::abs(i);
                    res(res.extent(0) - i_abs, j) = p[0];
                    res_jac(res.extent(0) - i_abs, j) = p[1];
                } else {
                    res(i, j) = p[0];
                    res_jac(i, j) = p[1];
                }
            }
        );
    }
}

template <typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void assoc_legendre_p_all(int type, T x, OutputMat1 res, OutputMat2 res_jac, OutputMat3 res_hess) {
    int m = (res.extent(0) - 1) / 2;
    int n = res.extent(1) - 1;

    for (int i = -m; i <= m; ++i) {
        assoc_legendre_p(
            n, i, type, x, assoc_legendre_p_diff_callback<T, 2>(),
            [res, res_jac, res_hess](int j, int i, int type, T z, const T(&p)[3], const T(&p_prev)[3]) {
                if (i < 0) {
                    int i_abs = std::abs(i);
                    res(res.extent(0) - i_abs, j) = p[0];
                    res_jac(res.extent(0) - i_abs, j) = p[1];
                    res_hess(res.extent(0) - i_abs, j) = p[2];
                } else {
                    res(i, j) = p[0];
                    res_jac(i, j) = p[1];
                    res_hess(i, j) = p[2];
                }
            }
        );
    }
}

template <typename T, typename Callable>
T sph_legendre_p(unsigned int n, unsigned int m, T phi, Callable callback) {
    if (m > n) {
        for (unsigned int j = 0; j < n; ++j) {
            callback(j, m, phi, 0);
        }

        return 0;
    }

    T x = std::cos(phi);

    long ls = (std::abs(x) > 1 ? -1 : 1);
    T xq = std::sqrt(ls * (1 - x * x));
    if (x < -1) {
        xq = -xq; // Ensure connection to the complex-valued function for |x| > 1
    }

    T p = 1 / (T(2) * std::sqrt(M_PI));
    for (unsigned int j = 0; j < m; ++j) {
        callback(j, m, phi, 0);

        T fac = std::sqrt(T(2 * j + 3) / T(2 * (j + 1)));
        p *= -ls * fac * xq;
    }

    callback(m, m, phi, p);

    if (m != n) {
        T p_prev = p;
        p = std::sqrt(T(2 * (m + 1) + 1)) * x * p_prev;
        callback(m + 1, m, phi, p);

        for (unsigned int j = m + 2; j <= n; ++j) {
            T fac_prev = std::sqrt(T(2 * j + 1) * T(4 * (j - 1) * (j - 1) - 1) / (T(2 * j - 3) * T(j * j - m * m)));
            T fac_prev_prev =
                std::sqrt(T(2 * j + 1) * T((j - 1) * (j - 1) - m * m) / (T(2 * j - 3) * T(j * j - m * m)));

            T p_prev_prev = p_prev;
            p_prev = p;
            p = fac_prev * x * p_prev - fac_prev_prev * p_prev_prev;
            callback(j, m, phi, p);
        }
    }

    return p;
}

template <typename T>
T sph_legendre_p(unsigned int n, unsigned int m, T phi) {
    return sph_legendre_p(n, m, phi, [](unsigned int j, unsigned int i, T phi, T value) {});
}

template <typename T, typename OutMat>
void sph_legendre_p_all(T phi, OutMat p) {
    unsigned int m = p.extent(0) - 1;
    unsigned int n = p.extent(1) - 1;

    for (unsigned int i = 0; i <= m; ++i) {
        sph_legendre(n, i, phi, [p](unsigned int j, unsigned int i, T phi, T value) { p(i, j) = value; });
    }
}

// ====================================================
// Purpose: Compute Legendre functions Qn(x) & Qn'(x)
// Input :  x  --- Argument of Qn(x)
//          n  --- Degree of Qn(x)  ( n = 0,1,2,…)
// Output:  QN(n) --- Qn(x)
//          QD(n) --- Qn'(x)
// ====================================================

template <typename T, typename OutputVec1, typename OutputVec2>
void lqn(T x, OutputVec1 qn, OutputVec2 qd) {
    int n = qn.size() - 1;

    T x2, q0, q1, qf, qc1, qc2, qr, qf0, qf1, qf2;
    const T eps = 1.0e-14;

    if (fabs(x) == 1.0) {
        for (int k = 0; k <= n; k++) {
            qn[k] = 1.0e300;
            qd[k] = 1.0e300;
        }
        return;
    }

    if (x <= 1.021) {
        x2 = fabs((1.0 + x) / (1.0 - x));
        q0 = 0.5 * log(x2);
        q1 = x * q0 - 1.0;
        qn[0] = q0;
        qn[1] = q1;
        qd[0] = 1.0 / (1.0 - x * x);
        qd[1] = qn[0] + x * qd[0];

        for (int k = 2; k <= n; k++) {
            qf = ((2.0 * k - 1.0) * x * q1 - (k - 1.0) * q0) / k;
            qn[k] = qf;
            qd[k] = (qn[k - 1] - x * qf) * k / (1.0 - x * x);
            q0 = q1;
            q1 = qf;
        }
    } else {
        qc1 = 0.0;
        qc2 = 1.0 / x;

        for (int j = 1; j <= n; j++) {
            qc2 *= j / ((2.0 * j + 1.0) * x);
            if (j == n - 1)
                qc1 = qc2;
        }

        for (int l = 0; l <= 1; l++) {
            int nl = n + l;
            qf = 1.0;
            qr = 1.0;

            for (int k = 1; k <= 500; k++) {
                qr = qr * (0.5 * nl + k - 1.0) * (0.5 * (nl - 1) + k) / ((nl + k - 0.5) * k * x * x);
                qf += qr;
                if (fabs(qr / qf) < eps)
                    break;
            }

            if (l == 0) {
                qn[n - 1] = qf * qc1;
            } else {
                qn[n] = qf * qc2;
            }
        }

        qf2 = qn[n];
        qf1 = qn[n - 1];

        for (int k = n; k >= 2; k--) {
            qf0 = ((2 * k - 1.0) * x * qf1 - k * qf2) / (k - 1.0);
            qn[k - 2] = qf0;
            qf2 = qf1;
            qf1 = qf0;
        }

        qd[0] = 1.0 / (1.0 - x * x);

        for (int k = 1; k <= n; k++) {
            qd[k] = k * (qn[k - 1] - x * qn[k]) / (1.0 - x * x);
        }
    }
}

// ==================================================
// Purpose: Compute the Legendre functions Qn(z) and
//          their derivatives Qn'(z) for a complex
//          argument
// Input :  x --- Real part of z
//          y --- Imaginary part of z
//          n --- Degree of Qn(z), n = 0,1,2,...
// Output:  CQN(n) --- Qn(z)
//          CQD(n) --- Qn'(z)
// ==================================================

template <typename T, typename OutputVec1, typename OutputVec2>
void lqn(std::complex<T> z, OutputVec1 cqn, OutputVec2 cqd) {
    int n = cqn.size() - 1;

    std::complex<T> cq0, cq1, cqf0 = 0.0, cqf1, cqf2;

    if (std::real(z) == 1) {
        for (int k = 0; k <= n; ++k) {
            cqn(k) = 1e300;
            cqd(k) = 1e300;
        }
        return;
    }
    int ls = ((std::abs(z) > 1.0) ? -1 : 1);

    cq0 = std::log(static_cast<T>(ls) * (static_cast<T>(1) + z) / (static_cast<T>(1) - z)) / static_cast<T>(2);
    cq1 = z * cq0 - static_cast<T>(1);

    cqn(0) = cq0;
    cqn(1) = cq1;

    if (std::abs(z) < 1.0001) {
        cqf0 = cq0;
        cqf1 = cq1;
        for (int k = 2; k <= n; k++) {
            cqf2 = (static_cast<T>(2 * k - 1) * z * cqf1 - static_cast<T>(k - 1) * cqf0) / static_cast<T>(k);
            cqn(k) = cqf2;
            cqf0 = cqf1;
            cqf1 = cqf2;
        }
    } else {
        int km;
        if (std::abs(z) > 1.1) {
            km = 40 + n;
        } else {
            km = (int) ((40 + n) * floor(-1.0 - 1.8 * log(std::abs(z - static_cast<T>(1)))));
        }

        cqf2 = 0.0;
        cqf1 = 1.0;
        for (int k = km; k >= 0; k--) {
            cqf0 = (static_cast<T>(2 * k + 3) * z * cqf1 - static_cast<T>(k + 2) * cqf2) / static_cast<T>(k + 1);
            if (k <= n) {
                cqn[k] = cqf0;
            }
            cqf2 = cqf1;
            cqf1 = cqf0;
        }
        for (int k = 0; k <= n; ++k) {
            cqn[k] *= cq0 / cqf0;
        }
    }
    cqd(0) = (cqn(1) - z * cqn(0)) / (z * z - static_cast<T>(1));

    for (int k = 1; k <= n; ++k) {
        cqd(k) = (static_cast<T>(k) * z * cqn(k) - static_cast<T>(k) * cqn(k - 1)) / (z * z - static_cast<T>(1));
    }
}

// ==========================================================
// Purpose: Compute the associated Legendre functions of the
//          second kind, Qmn(x) and Qmn'(x)
// Input :  x  --- Argument of Qmn(x)
//          m  --- Order of Qmn(x)  ( m = 0,1,2,… )
//          n  --- Degree of Qmn(x) ( n = 0,1,2,… )
//          mm --- Physical dimension of QM and QD
// Output:  QM(m,n) --- Qmn(x)
//          QD(m,n) --- Qmn'(x)
// ==========================================================

template <typename T, typename OutputMat1, typename OutputMat2>
void lqmn(T x, OutputMat1 qm, OutputMat2 qd) {
    int m = qm.extent(0) - 1;
    int n = qm.extent(1) - 1;

    double q0, q1, q10, qf, qf0, qf1, qf2, xs, xq;
    int i, j, k, km, ls;

    if (fabs(x) == 1.0) {
        for (i = 0; i < (m + 1); i++) {
            for (j = 0; j < (n + 1); j++) {
                qm(i, j) = 1e300;
                qd(i, j) = 1e300;
            }
        }
        return;
    }
    ls = 1;
    if (fabs(x) > 1.0) {
        ls = -1;
    }
    xs = ls * (1.0 - x * x);
    xq = sqrt(xs);
    q0 = 0.5 * log(fabs((x + 1.0) / (x - 1.0)));
    if (fabs(x) < 1.0001) {
        qm(0, 0) = q0;
        qm(0, 1) = x * q0 - 1.0;
        qm(1, 0) = -1.0 / xq;
        qm(1, 1) = -ls * xq * (q0 + x / (1. - x * x));
        for (i = 0; i <= 1; i++) {
            for (j = 2; j <= n; j++) {
                qm(i, j) = ((2.0 * j - 1.) * x * qm(i, j - 1) - (j + i - 1) * qm(i, j - 2)) / (j - i);
            }
        }
        /* 15 */
        for (i = 2; i <= m; i++) {
            for (j = 0; j <= n; j++) {
                qm(i, j) = -2.0 * (i - 1.0) * x / xq * qm(i - 1, j) - ls * (j + i - 1.0) * (j - i + 2.0) * qm(i - 2, j);
            }
        }
    } else {
        if (fabs(x) > 1.1) {
            km = 40 + m + n;
        } else {
            km = (40 + m + n) * ((int) (-1. - 1.8 * log(x - 1.)));
        }
        qf2 = 0.0;
        qf1 = 1.0;
        qf0 = 0.0;
        for (k = km; k >= 0; k--) {
            qf0 = ((2.0 * k + 3.0) * x * qf1 - (k + 2.0) * qf2) / (k + 1.0);
            if (k <= n) {
                qm(0, k) = qf0;
            }
            qf2 = qf1;
            qf1 = qf0;
        }

        for (k = 0; k <= n; k++) {
            qm(0, k) *= q0 / qf0;
        }

        qf2 = 0.0;
        qf1 = 1.0;
        for (k = km; k >= 0; k--) {
            qf0 = ((2.0 * k + 3.0) * x * qf1 - (k + 1.0) * qf2) / (k + 2.0);
            if (k <= n) {
                qm(1, k) = qf0;
            }
            qf2 = qf1;
            qf1 = qf0;
        }

        q10 = -1.0 / xq;
        for (k = 0; k <= n; k++) {
            qm(1, k) *= q10 / qf0;
        }

        for (j = 0; j <= n; j++) {
            q0 = qm(0, j);
            q1 = qm(1, j);
            for (i = 0; i <= (m - 2); i++) {
                qf = -2. * (i + 1.) * x / xq * q1 + (j - i) * (j + i + 1.) * q0;
                qm(i + 2, j) = qf;
                q0 = q1;
                q1 = qf;
            }
        }
    }

    qd(0, 0) = ls / xs;
    for (j = 1; j <= n; j++) {
        qd(0, j) = ls * j * (qm(0, j - 1) - x * qm(0, j)) / xs;
    }

    for (i = 1; i <= m; i++) {
        for (j = 0; j <= n; j++) {
            qd(i, j) = ls * i * x / xs * qm(i, j) + (i + j) * (j - i + 1.) / xq * qm(i - 1, j);
        }
    }
}

// =======================================================
// Purpose: Compute the associated Legendre functions of
//          the second kind, Qmn(z) and Qmn'(z), for a
//          complex argument
// Input :  x  --- Real part of z
//          y  --- Imaginary part of z
//          m  --- Order of Qmn(z)  ( m = 0,1,2,… )
//          n  --- Degree of Qmn(z) ( n = 0,1,2,… )
//          mm --- Physical dimension of CQM and CQD
// Output:  CQM(m,n) --- Qmn(z)
//          CQD(m,n) --- Qmn'(z)
// =======================================================

template <typename T, typename OutputMat1, typename OutputMat2>
void lqmn(std::complex<T> z, OutputMat1 cqm, OutputMat2 cqd) {
    int m = cqm.extent(0) - 1;
    int n = cqm.extent(1) - 1;

    int i, j, k, km, ls;
    std::complex<T> cq0, cq1, cq10, cqf0 = 0, cqf, cqf1, cqf2, zq, zs;

    if ((std::abs(std::real(z)) == 1) && (std::imag(z) == 0)) {
        for (i = 0; i < (m + 1); i++) {
            for (j = 0; j < (n + 1); j++) {
                cqm(i, j) = 1e300;
                cqd(i, j) = 1e300;
            }
        }

        return;
    }

    T xc = std::abs(z);
    ls = 0;
    if ((std::imag(z) == 0) || (xc < 1)) {
        ls = 1;
    }
    if (xc > 1) {
        ls = -1;
    }
    zs = static_cast<T>(ls) * (static_cast<T>(1) - z * z);
    zq = std::sqrt(zs);

    cq0 = std::log(static_cast<T>(ls) * (static_cast<T>(1) + z) / (static_cast<T>(1) - z)) / static_cast<T>(2);
    if (xc < 1.0001) {
        cqm(0, 0) = cq0;
        cqm(1, 0) = -static_cast<T>(1) / zq;
        cqm(0, 1) = z * cq0 - static_cast<T>(1);
        cqm(1, 1) = -zq * (cq0 + z / (static_cast<T>(1) - z * z));

        for (i = 0; i <= 1; i++) {
            for (j = 2; j <= n; j++) {
                cqm(i, j) =
                    (static_cast<T>(2 * j - 1) * z * cqm(i, j - 1) - static_cast<T>(j + i - 1) * cqm(i, j - 2)) /
                    static_cast<T>(j - i);
            }
        }

        for (i = 2; i <= m; i++) {
            for (j = 0; j <= n; j++) {
                cqm(i, j) = -2 * static_cast<T>(i - 1) * z / zq * cqm(i - 1, j) -
                            static_cast<T>(ls * (j + i - 1) * (j - i + 2)) * cqm(i - 2, j);
            }
        }
    } else {
        if (xc > 1.1) {
            km = 40 + m + n;
        } else {
            km = (40 + m + n) * ((int) (-1.0 - 1.8 * log(xc - 1.)));
        }
        cqf2 = 0.0;
        cqf1 = 1.0;
        for (k = km; k >= 0; k--) {
            cqf0 = (static_cast<T>(2 * k + 3) * z * cqf1 - static_cast<T>(k + 2) * cqf2) / static_cast<T>(k + 1);
            if (k <= n) {
                cqm(0, k) = cqf0;
            }
            cqf2 = cqf1;
            cqf1 = cqf0;
        }

        for (k = 0; k <= n; k++) {
            cqm(0, k) *= cq0 / cqf0;
        }

        cqf2 = 0.0;
        cqf1 = 1.0;
        for (k = km; k >= 0; k--) {
            cqf0 = (static_cast<T>(2 * k + 3) * z * cqf1 - static_cast<T>(k + 1) * cqf2) / static_cast<T>(k + 2);
            if (k <= n) {
                cqm(1, k) = cqf0;
            }
            cqf2 = cqf1;
            cqf1 = cqf0;
        }

        cq10 = -static_cast<T>(1) / zq;
        for (k = 0; k <= n; k++) {
            cqm(1, k) *= cq10 / cqf0;
        }

        for (j = 0; j <= n; j++) {
            cq0 = cqm(0, j);
            cq1 = cqm(1, j);
            for (i = 0; i <= (m - 2); i++) {
                cqf = -static_cast<T>(2 * (i + 1)) * z / zq * cq1 + static_cast<T>((j - i) * (j + i + 1)) * cq0;
                cqm(i + 2, j) = cqf;
                cq0 = cq1;
                cq1 = cqf;
            }
        }

        cqd(0, 0) = static_cast<T>(ls) / zs;
        for (j = 1; j <= n; j++) {
            cqd(0, j) = ls * static_cast<T>(j) * (cqm(0, j - 1) - z * cqm(0, j)) / zs;
        }

        for (i = 1; i <= m; i++) {
            for (j = 0; j <= n; j++) {
                cqd(i, j) = static_cast<T>(ls * i) * z / zs * cqm(i, j) +
                            static_cast<T>((i + j) * (j - i + 1)) / zq * cqm(i - 1, j);
            }
        }
    }
}

} // namespace special
