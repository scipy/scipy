#pragma once

#include <iostream>

#include "error.h"
#include "gamma.h"
#include "recur.h"

namespace special {

template <typename T>
struct legendre_p_recurrence {
    T z;

    void operator()(int n, T (&res)[2]) const {
        res[0] = -T(n - 1) / T(n);
        res[1] = T(2 * n - 1) * z / T(n);
    }

    void operator()(int n, T (&res)[2], T (&res_jac)[2]) const {
        (*this)(n, res);

        res_jac[0] = 0;
        res_jac[1] = T(2 * n - 1) / T(n);
    }

    void operator()(int n, T (&res)[2], T (&res_jac)[2], T (&res_hess)[2]) const {
        (*this)(n, res, res_jac);

        res_hess[0] = 0;
        res_hess[1] = 0;
    }

    void init(T *res) {
        res[0] = 1;
        res[1] = z;
    }

    void init(T *res, T *res_jac) {
        init(res);

        res_jac[0] = 0;
        res_jac[1] = 1;
    }

    void init(T *res, T *res_jac, T *res_hess) {
        init(res, res_jac);

        res_hess[0] = 0;
        res_hess[1] = 0;
    }
};

/**
 * Compute the Legendre polynomial of degree n.
 *
 * @param n degree of the polynomial
 * @param z argument of the polynomial, either real or complex
 * @param callback a function to be called as callback(j, z, p, args...) for 0 <= j <= n
 * @param args arguments to forward to the callback
 *
 * @return value of the polynomial
 */
template <typename T, typename Callable>
void legendre_p_recur(int n, T z, T (&res)[3], Callable callback) {
    legendre_p_recurrence<T> r{z};
    r.init(res);

    forward_recur(r, res, 0, n + 1, callback);
}

template <typename T, typename Callable>
void legendre_p_recur(int n, T z, T (&res)[3], T (&res_jac)[3], Callable callback) {
    legendre_p_recurrence<T> r{z};
    r.init(res, res_jac);

    forward_recur(r, res, res_jac, 0, n + 1, callback);
}

template <typename T, typename Callable>
void legendre_p_recur(int n, T z, T (&res)[3], T (&res_jac)[3], T (&res_hess)[3], Callable callback) {
    legendre_p_recurrence<T> r{z};
    r.init(res, res_jac, res_hess);

    forward_recur(r, res, res_jac, res_hess, 0, n + 1, callback);
}

/**
 * Compute the Legendre polynomial of degree n.
 *
 * @param n degree of the polynomial
 * @param z argument of the polynomial, either real or complex
 *
 * @return value of the polynomial
 */
template <typename T>
T legendre_p(int n, T z) {
    T p[3] = {};
    legendre_p_recur(n, z, p, [](int j, auto r, const T(&p)[3]) {});

    return p[2];
}

/**
 * Compute the Legendre polynomial of degree n.
 *
 * @param n degree of the polynomial
 * @param z argument of the polynomial, either real or complex
 * @param res value of the polynomial
 * @param res_jac value of the polynomial first derivative
 */
template <typename T>
void legendre_p(int n, T z, T &res, T &res_jac) {
    T p[3] = {};
    T p_jac[3] = {};
    legendre_p_recur(n, z, p, p_jac, [](int j, auto r, const T(&p)[3], const T(&p_jac)[3]) {});

    res = p[2];
    res_jac = p_jac[2];
}

/**
 * Compute the Legendre polynomial of degree n.
 *
 * @param n degree of the polynomial
 * @param z argument of the polynomial, either real or complex
 * @param res value of the polynomial
 * @param res_jac value of the polynomial first derivative
 * @param res_hess value of the polynomial second derivative
 */
template <typename T>
void legendre_p(int n, T z, T &res, T &res_jac, T &res_hess) {
    T p[3] = {};
    T p_jac[3] = {};
    T p_hess[3] = {};
    legendre_p_recur(
        n, z, p, p_jac, p_hess, [](int j, auto r, const T(&p)[3], const T(&p_jac)[3], const T(&p_hess)[3]) {}
    );

    res = p[2];
    res_jac = p_jac[2];
    res_hess = p_hess[2];
}

/**
 * Compute all Legendre polynomials of degree j, where 0 <= j <= n.
 *
 * @param z argument of the polynomials, either real or complex
 * @param res a view into a multidimensional array with element type T and size n + 1 to store the value of each
 *            polynomial
 */
template <typename T, typename OutputVec>
void legendre_p_all(T z, OutputVec res) {
    int n = res.extent(0) - 1;

    T p[3] = {};
    legendre_p_recur(n, z, p, [res](int j, auto r, const T(&p)[3]) { res(j) = p[2]; });
}

/**
 * Compute all Legendre polynomials of degree j, where 0 <= j <= n.
 *
 * @param z argument of the polynomials, either real or complex
 * @param res a view into a multidimensional array with element type T and size n + 1 to store the value of each
 *            polynomial
 * @param res_jac same as res, but for the polynomial first derivatives
 */
template <typename T, typename OutputVec1, typename OutputVec2>
void legendre_p_all(T z, OutputVec1 res, OutputVec2 res_jac) {
    int n = res.extent(0) - 1;

    T p[3] = {};
    T p_jac[3] = {};
    legendre_p_recur(n, z, p, p_jac, [res, res_jac](int j, auto r, const T(&p)[3], const T(&p_jac)[3]) {
        res(j) = p[2];
        res_jac(j) = p_jac[2];
    });
}

/**
 * Compute all Legendre polynomials of degree j, where 0 <= j <= n.
 *
 * @param z argument of the polynomials, either real or complex
 * @param res a view into a multidimensional array with element type T and size n + 1 to store the value of each
 *            polynomial
 * @param res_jac same as res, but for the polynomial first derivatives
 * @param res_hess same as res, but for the polynomial second derivatives
 */
template <typename T, typename OutputVec1, typename OutputVec2, typename OutputVec3>
void legendre_p_all(T z, OutputVec1 res, OutputVec2 res_jac, OutputVec3 res_hess) {
    int n = res.extent(0) - 1;

    T p[3] = {};
    T p_jac[3] = {};
    T p_hess[3] = {};
    legendre_p_recur(
        n, z, p, p_jac, p_hess,
        [res, res_jac, res_hess](int j, auto r, const T(&p)[3], const T(&p_jac)[3], const T(&p_hess)[3]) {
            res(j) = p[2];
            res_jac(j) = p_jac[2];
            res_hess(j) = p_hess[2];
        }
    );
}

template <typename T>
struct diag_recurrence {
    int m;
    int type;
    T z;
    T type_sign;

    diag_recurrence(int m, int type, T z) : m(m), type(type), z(z) {
        if (type == 3) {
            type_sign = -1;
        } else {
            type_sign = 1;
        }
    }

    void operator()(int n, T (&res)[2]) {
        T fac = type_sign * T((2 * n - 1) * (2 * n - 3));

        // other square roots can be avoided if each iteration increments by 2

        res[0] = fac * (T(1) - z * z);
        res[1] = 0;
    }

    void operator()(int n, T (&res)[2], T (&res_jac)[2]) {
        (*this)(n, res);

        T fac = type_sign * T((2 * n - 1) * (2 * n - 3));

        res_jac[0] = -T(2) * fac * z;
        res_jac[1] = 0;
    }

    void operator()(int n, T (&res)[2], T (&res_jac)[2], T (&res_hess)[2]) {
        (*this)(n, res, res_jac);

        T fac = type_sign * T((2 * n - 1) * (2 * n - 3));

        res_hess[0] = -T(2) * fac;
        res_hess[1] = 0;
    }

    void init2(T (&res)[3]) {
        res[0] = 1;
        if (type == 3) {
            res[1] = std::sqrt(z * z - T(1)); // do not modify, see function comment
            if (std::real(z) < 0) {
                res[1] = -res[1];
            }
        } else {
            res[1] = -std::sqrt(T(1) - z * z); // do not modify, see function comment
        }
    }

    void init(T (&res)[3]) {
        init2(res);

        int m_abs = std::abs(m);
        bool m_odd = m_abs % 2;

        T fac = 1;
        if (m < 0) {
            fac *= std::pow(-1, m_abs) / std::tgamma(2 * m_abs + 1);
            if (m_odd && type == 3) {
                fac *= -1;
            }
        }

        for (size_t j = 0; j < 3; ++j) {
            res[j] *= fac;
        }
    }

    void init2(T (&res)[3], T (&res_jac)[3]) {
        init2(res);

        T type_sign;
        if (type == 3) {
            type_sign = -1;
        } else {
            type_sign = 1;
        }

        res_jac[0] = 0;
        res_jac[1] = -type_sign * z / res[1];
    }

    void init(T (&res)[3], T (&res_jac)[3]) {
        init2(res, res_jac);

        int m_abs = std::abs(m);
        bool m_odd = m_abs % 2;

        T fac = 1;
        if (m < 0) {
            fac *= std::pow(-1, m_abs) / std::tgamma(2 * m_abs + 1);
            if (m_odd && type == 3) {
                fac *= -1;
            }
        }

        for (size_t j = 0; j < 3; ++j) {
            res[j] *= fac;
            res_jac[j] *= fac;
        }
    }

    void init2(T (&res)[3], T (&res_jac)[3], T (&res_hess)[3]) {
        init2(res, res_jac);

        res_hess[0] = 0;
        res_hess[1] = -type_sign * (T(1) - z * res_jac[1] / res[1]) / res[1];
    }

    void init(T (&res)[3], T (&res_jac)[3], T (&res_hess)[3]) {
        init2(res, res_jac, res_hess);

        int m_abs = std::abs(m);
        bool m_odd = m_abs % 2;

        T fac = 1;
        if (m < 0) {
            fac *= std::pow(-1, m_abs) / std::tgamma(2 * m_abs + 1);
            if (m_odd && type == 3) {
                fac *= -1;
            }
        }

        for (size_t j = 0; j < 3; ++j) {
            res[j] *= fac;
            res_jac[j] *= fac;
            res_hess[j] *= fac;
        }
    }
};

template <typename T>
struct assoc_legendre_p_recurrence {
    int m;
    int type;
    T z;
    T type_sign;

    assoc_legendre_p_recurrence(int m, int type, T z) : m(m), type(type), z(z) {
        if (type == 3) {
            type_sign = -1;
        } else {
            type_sign = 1;
        }
    }

    bool is_pm1() const { return std::abs(std::real(z)) == 1 && std::imag(z) == 0; }

    void operator()(int n, T (&res)[2]) const {
        res[0] = -T(n + m - 1) / T(n - m);
        res[1] = T(2 * n - 1) * z / T(n - m);
    }

    void operator()(int n, T (&res)[2], T (&res_jac)[2]) const {
        (*this)(n, res);

        res_jac[0] = 0;
        res_jac[1] = T(2 * n - 1) / T(n - m);
    }

    void operator()(int n, T (&res)[2], T (&res_jac)[2], T (&res_hess)[2]) const {
        (*this)(n, res, res_jac);

        res_hess[0] = 0;
        res_hess[1] = 0;
    }

    void init(T (&res)[3]) const {
        int m_abs = std::abs(m);

        diag_recurrence<T> r{m, type, z};

        T vals[3];
        r.init(vals);

        forward_recur(r, vals, 0, m_abs + 1, [](auto, auto, auto) {});

        res[0] = vals[2];
        res[1] = T(2 * (m_abs + 1) - 1) * z * res[0] / T(m_abs + 1 - m);
    }

    void init(T (&res)[3], T (&res_jac)[3]) const {
        int m_abs = std::abs(m);

        diag_recurrence<T> r{m, type, z};

        T vals[3];
        T vals_jac[3];
        r.init(vals, vals_jac);

        forward_recur(r, vals, vals_jac, 0, m_abs + 1, [](auto, auto, auto, auto) {});

        res[0] = vals[2];
        res[1] = T(2 * (m_abs + 1) - 1) * z * res[0] / T(m_abs + 1 - m);

        res_jac[0] = vals_jac[2];
        res_jac[1] = T(2 * (m_abs + 1) - 1) * (res[0] + z * res_jac[0]) / T(m_abs + 1 - m);
    }

    void init(T (&res)[3], T (&res_jac)[3], T (&res_hess)[3]) const {
        int m_abs = std::abs(m);

        diag_recurrence<T> r{m, type, z};

        T vals[3];
        T vals_jac[3];
        T vals_hess[3];
        r.init(vals, vals_jac, vals_hess);

        forward_recur(r, vals, vals_jac, vals_hess, 0, m_abs + 1, [](auto, auto, auto, auto, auto) {});

        res[0] = vals[2];
        res[1] = T(2 * (m_abs + 1) - 1) * z * res[0] / T(m_abs + 1 - m);

        res_jac[0] = vals_jac[2];
        res_jac[1] = T(2 * (m_abs + 1) - 1) * (res[0] + z * res_jac[0]) / T(m_abs + 1 - m);

        res_hess[0] = vals_hess[2];
        res_hess[1] = T(2 * (m_abs + 1) - 1) * (res_jac[0] + res_jac[0] + z * res_hess[0]) / T(m_abs + 1 - m);
    }

    void at_pm1(int n, T &res) {
        if (m == 0) {
            res = 1;
        } else {
            res = 0;
        }
    }

    void at_pm1(int n, T &res, T &res_jac) {
        at_pm1(n, res);

        if (std::abs(m) > n) {
            res_jac = 0;
        } else if (m == 0) {
            res_jac = T(n) * T(n + 1) * std::pow(z, T(n + 1)) / T(2);
        } else if (m == 1) {
            res_jac = std::pow(z, T(n)) * std::numeric_limits<remove_complex_t<T>>::infinity();
        } else if (m == 2) {
            res_jac = -type_sign * T(n + 2) * T(n + 1) * T(n) * T(n - 1) * std::pow(z, T(n + 1)) / T(4);
        } else if (m == -2) {
            res_jac = -type_sign * std::pow(z, T(n + 1)) / T(4);
        } else if (m == -1) {
            res_jac = -std::pow(z, T(n)) * std::numeric_limits<remove_complex_t<T>>::infinity();
        } else {
            res_jac = 0;
        }
    }

    void at_pm1(int n, T &res, T &res_jac, T &res_hess) {
        at_pm1(n, res, res_jac);

        // need to complete these
        if (std::abs(m) > n) {
            res_hess = 0;
        } else if (m == 0) {
            res_hess = T(n + 2) * T(n + 1) * T(n) * T(n - 1) / T(8);
        } else if (m == 1) {
            res_hess = std::numeric_limits<remove_complex_t<T>>::infinity();
        } else if (m == 2) {
            res_hess = -T((n + 1) * n - 3) * T(n + 2) * T(n + 1) * T(n) * T(n - 1) / T(12);
        } else if (m == 3) {
            res_hess = std::numeric_limits<remove_complex_t<T>>::infinity();
        } else if (m == 4) {
            res_hess = T(n + 4) * T(n + 3) * T(n + 2) * T(n + 1) * T(n) * T(n - 1) * T(n - 2) * T(n - 3) / T(48);
        } else if (m == -4) {
            res_hess = 0;
        } else if (m == -3) {
            res_hess = -std::numeric_limits<remove_complex_t<T>>::infinity();
        } else if (m == -2) {
            res_hess = -T(1) / T(4);
        } else if (m == -1) {
            res_hess = -std::numeric_limits<remove_complex_t<T>>::infinity();
        } else {
            res_hess = 0;
        }
    }
};

/**
 * Compute the associated Legendre polynomial of degree n and order n.
 *
 * We need to be careful with complex arithmetic, in particular the square roots
 * should not be modified. This is because the sign bit of a real or imaginary part,
 * even if it is equal to zero, can affect the branch cut.
 */

/**
 * Compute the associated Legendre polynomial of degree n and order m.
 *
 * @param n degree of the polynomial
 * @param m order of the polynomial
 * @param type specifies the branch cut of the polynomial, either 1, 2, or 3
 * @param z argument of the polynomial, either real or complex
 * @param callback a function to be called as callback(j, m, type, z, p, p_prev, args...) for 0 <= j <= n
 * @param args arguments to forward to the callback
 *
 * @return value of the polynomial
 */
template <typename T, typename Callable>
void assoc_legendre_p_recur(int n, int m, int type, T z, T (&res)[3], Callable callback) {
    assoc_legendre_p_recurrence<T> r{m, type, z};

    int m_abs = std::abs(m);
    if (m_abs > n) {
        for (int j = 0; j <= n; ++j) {
            callback(j, r, res);
        }
    } else {
        for (int j = 0; j < m_abs; ++j) {
            callback(j, r, res);
        }

        if (r.is_pm1()) {
            for (int j = m_abs; j <= n; ++j) {
                forward_recur_shift(res);

                r.at_pm1(j, res[2]);
                callback(j, r, res);
            }
        } else {
            r.init(res);
            forward_recur(r, res, m_abs, n + 1, callback);
        }
    }
}

template <typename T, typename Callable>
void assoc_legendre_p_recur(int n, int m, int type, T z, T (&res)[3], T (&res_jac)[3], Callable callback) {
    assoc_legendre_p_recurrence<T> r{m, type, z};

    int m_abs = std::abs(m);
    if (m_abs > n) {
        for (int j = 0; j <= n; ++j) {
            callback(j, r, res, res_jac);
        }
    } else {
        for (int j = 0; j < m_abs; ++j) {
            callback(j, r, res, res_jac);
        }

        if (r.is_pm1()) {
            for (int j = m_abs; j <= n; ++j) {
                forward_recur_shift(res);
                forward_recur_shift(res_jac);

                r.at_pm1(j, res[2], res_jac[2]);
                callback(j, r, res, res_jac);
            }
        } else {
            r.init(res, res_jac);
            forward_recur(r, res, res_jac, m_abs, n + 1, callback);
        }
    }
}

template <typename T, typename Callable>
void assoc_legendre_p_recur(
    int n, int m, int type, T z, T (&res)[3], T (&res_jac)[3], T (&res_hess)[3], Callable callback
) {
    assoc_legendre_p_recurrence<T> r{m, type, z};

    int m_abs = std::abs(m);
    if (m_abs > n) {
        for (int j = 0; j <= n; ++j) {
            callback(j, r, res, res_jac, res_hess);
        }
    } else {
        for (int j = 0; j < m_abs; ++j) {
            callback(j, r, res, res_jac, res_hess);
        }

        if (r.is_pm1()) {
            for (int j = m_abs; j <= n; ++j) {
                forward_recur_shift(res);
                forward_recur_shift(res_jac);
                forward_recur_shift(res_hess);

                r.at_pm1(j, res[2], res_jac[2], res_hess[2]);
                callback(j, r, res, res_jac, res_hess);
            }
        } else {
            r.init(res, res_jac, res_hess);
            forward_recur(r, res, res_jac, res_hess, m_abs, n + 1, callback);
        }
    }
}

/**
 * Compute the associated Legendre polynomial of degree n and order m.
 *
 * @param n degree of the polynomial
 * @param m order of the polynomial
 * @param type specifies the branch cut of the polynomial, either 1, 2, or 3
 * @param z argument of the polynomial, either real or complex
 *
 * @return value of the polynomial
 */
template <typename T>
T assoc_legendre_p(int n, int m, int type, T z) {
    T p[3] = {};
    assoc_legendre_p_recur(n, m, type, z, p, [](int j, const auto &r, const T(&p)[3]) {});

    return p[2];
}

template <typename T>
void assoc_legendre_p(int n, int m, int type, T z, T &res, T &res_jac) {
    T p[3] = {};
    T p_jac[3] = {};
    assoc_legendre_p_recur(n, m, type, z, p, p_jac, [](int j, const auto &r, const T(&p)[3], const T(&p_jac)[3]) {});

    res = p[2];
    res_jac = p_jac[2];
}

template <typename T>
void assoc_legendre_p(int n, int m, int type, T z, T &res, T &res_jac, T &res_hess) {
    T p[3] = {};
    T p_jac[3] = {};
    T p_hess[3] = {};
    assoc_legendre_p_recur(
        n, m, type, z, p, p_jac, p_hess,
        [](int j, const auto &r, const T(&p)[3], const T(&p_jac)[3], const T(&p_hess)[3]) {}
    );

    res = p[2];
    res_jac = p_jac[2];
    res_hess = p_hess[2];
}

/**
 * Compute all associated Legendre polynomials of degree j and order i, where 0 <= j <= n and -m <= i <= m.
 *
 * @param type specifies the branch cut of the polynomial, either 1, 2, or 3
 * @param z argument of the polynomial, either real or complex
 * @param res a view into the output with element type T and extents (2 * m + 1, n + 1)
 *
 * @return value of the polynomial
 */
template <typename T, typename OutputMat>
void assoc_legendre_p_all(int type, T z, OutputMat res) {
    int m = (res.extent(0) - 1) / 2;
    int n = res.extent(1) - 1;

    for (int i = -m; i <= m; ++i) {
        T p[3] = {};
        assoc_legendre_p_recur(n, i, type, z, p, [res](int j, auto r, const T(&p)[3]) {
            int i_offset;
            if (r.m < 0) {
                i_offset = res.extent(0);
            } else {
                i_offset = 0;
            }

            res(r.m + i_offset, j) = p[2];
        });
    }
}

template <typename T, typename OutputMat1, typename OutputMat2>
void assoc_legendre_p_all(int type, T z, OutputMat1 res, OutputMat2 res_jac) {
    int m = (res.extent(0) - 1) / 2;
    int n = res.extent(1) - 1;

    for (int i = -m; i <= m; ++i) {
        T p[3] = {};
        T p_jac[3] = {};
        assoc_legendre_p_recur(
            n, i, type, z, p, p_jac,
            [res, res_jac](int j, auto r, const T(&p)[3], const T(&p_jac)[3]) {
                int i_offset;
                if (r.m < 0) {
                    i_offset = res.extent(0);
                } else {
                    i_offset = 0;
                }

                res(r.m + i_offset, j) = p[2];
                res_jac(r.m + i_offset, j) = p_jac[2];
            }
        );
    }
}

template <typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void assoc_legendre_p_all(int type, T z, OutputMat1 res, OutputMat2 res_jac, OutputMat3 res_hess) {
    int m = (res.extent(0) - 1) / 2;
    int n = res.extent(1) - 1;

    for (int i = -m; i <= m; ++i) {
        T p[3] = {};
        T p_jac[3] = {};
        T p_hess[3] = {};
        assoc_legendre_p_recur(
            n, i, type, z, p, p_jac, p_hess,
            [res, res_jac, res_hess](int j, auto r, const T(&p)[3], const T(&p_jac)[3], const T(&p_hess)[3]) {
                int i_offset;
                if (r.m < 0) {
                    i_offset = res.extent(0);
                } else {
                    i_offset = 0;
                }

                res(r.m + i_offset, j) = p[2];
                res_jac(r.m + i_offset, j) = p_jac[2];
                res_hess(r.m + i_offset, j) = p_hess[2];
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
