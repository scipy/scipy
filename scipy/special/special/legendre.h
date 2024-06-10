#pragma once

#include "error.h"
#include "recur.h"
#include "tuple_utility.h"

namespace special {

template <typename T>
struct legendre_p_initializer_n {
    T z;

    void operator()(std::tuple<T (&)[2]> res) const { tuple_assign(res, {{1, z}}); }

    void operator()(std::tuple<T (&)[2], T (&)[2]> res) const { tuple_assign(res, {{1, z}, {0, 1}}); }

    void operator()(std::tuple<T (&)[2], T (&)[2], T (&)[2]> res) const { tuple_assign(res, {{1, z}, {0, 1}, {0, 0}}); }
};

template <typename T>
struct legendre_p_recurrence_n {
    T z;

    void operator()(int n, std::tuple<T (&)[2]> res) const {
        T fac0 = -T(n - 1) / T(n);
        T fac1 = T(2 * n - 1) / T(n);

        tuple_assign(res, {{fac0, fac1 * z}});
    }

    void operator()(int n, std::tuple<T (&)[2], T (&)[2]> res) const {
        T fac0 = -T(n - 1) / T(n);
        T fac1 = T(2 * n - 1) / T(n);

        tuple_assign(res, {{fac0, fac1 * z}, {0, fac1}});
    }

    void operator()(int n, std::tuple<T (&)[2], T (&)[2], T (&)[2]> res) const {
        T fac0 = -T(n - 1) / T(n);
        T fac1 = T(2 * n - 1) / T(n);

        tuple_assign(res, {{fac0, fac1 * z}, {0, fac1}, {0, 0}});
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
template <typename T, typename... OutputVals, typename Func>
void legendre_p_for_each_n(int n, T z, std::tuple<OutputVals (&)[2]...> res, Func f) {
    legendre_p_initializer_n<T> init_n{z};
    init_n(res);

    legendre_p_recurrence_n<T> re_n{z};
    forward_recur(0, n + 1, re_n, res, f);
}

/**
 * Compute the Legendre polynomial of degree n.
 *
 * @param n degree of the polynomial
 * @param z argument of the polynomial, either real or complex
 * @param res result (gradient)
 */
template <typename T, typename... OutputVals>
void legendre_p(int n, T z, std::tuple<OutputVals &...> res) {
    std::tuple<OutputVals[2]...> res_n;
    legendre_p_for_each_n(n, z, tuple_ref_each(res_n), [](int n, std::tuple<OutputVals(&)[2]...> res_n) {});

    res = tuple_access_each(res_n, 1);
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
    T res;
    legendre_p(n, z, std::tie(res));

    return res;
}

/**
 * Compute all Legendre polynomials of degree j, where 0 <= j <= n.
 *
 * @param z argument of the polynomials, either real or complex
 * @param res a view into a multidimensional array with element type T and size n + 1 to store the value of each
 *            polynomial
 */
template <typename T, typename... OutputVecs>
void legendre_p_all(T z, std::tuple<OutputVecs &...> res) {
    auto &res0 = std::get<0>(res);
    int n = res0.extent(0) - 1;

    std::tuple<typename OutputVecs::value_type[2]...> res_n;
    legendre_p_for_each_n(n, z, tuple_ref_each(res_n), [&res](int n, auto res_n) {
        tuple_call_each(res, n) = tuple_access_each(res_n, 1);
    });
}

struct assoc_legendre_unnorm_policy {};

struct assoc_legendre_norm_policy {};

constexpr assoc_legendre_unnorm_policy assoc_legendre_unnorm;

constexpr assoc_legendre_norm_policy assoc_legendre_norm;

template <typename T, typename NormPolicy>
struct assoc_legendre_p_initializer_m_m_abs;

template <typename T>
struct assoc_legendre_p_initializer_m_m_abs<T, assoc_legendre_unnorm_policy> {
    bool m_signbit;
    int type;
    T z;
    T type_sign;
    T w;

    assoc_legendre_p_initializer_m_m_abs(bool m_signbit, int type, T z) : m_signbit(m_signbit), type(type), z(z) {
        if (type == 3) {
            type_sign = -1;

            w = std::sqrt(z * z - T(1)); // do not modify, see function comment
            if (std::real(z) < 0) {
                w = -w;
            }
        } else {
            type_sign = 1;

            w = -std::sqrt(T(1) - z * z); // do not modify, see function comment
            if (m_signbit) {
                w = -w;
            }
        }
    }

    void operator()(std::tuple<T (&)[2]> res) const {
        tuple_assign(res, {{1, w}});

        if (m_signbit) {
            std::get<0>(res)[1] /= 2;
        }
    }

    void operator()(std::tuple<T (&)[2], T (&)[2]> res) const {
        tuple_assign(res, {{1, w}, {0, -type_sign * z / w}});

        if (m_signbit) {
            std::get<0>(res)[1] /= 2;
            std::get<1>(res)[1] /= 2;
        }
    }

    void operator()(std::tuple<T (&)[2], T (&)[2], T (&)[2]> res) const {
        tuple_assign(res, {{1, w}, {0, -type_sign * z / w}, {0, T(1) / ((z * z - T(1)) * w)}});

        if (m_signbit) {
            std::get<0>(res)[1] /= 2;
            std::get<1>(res)[1] /= 2;
            std::get<2>(res)[1] /= 2;
        }
    }
};

template <typename T>
struct assoc_legendre_p_initializer_m_m_abs<T, assoc_legendre_norm_policy> {
    bool m_signbit;
    int type;
    T z;
    T type_sign;
    T w;

    assoc_legendre_p_initializer_m_m_abs(bool m_signbit, int type, T z) : m_signbit(m_signbit), type(type), z(z) {
        if (type == 3) {
            type_sign = -1;

            w = std::sqrt(z * z - T(1)); // do not modify, see function comment
            if (std::real(z) < 0) {
                w = -w;
            }
        } else {
            type_sign = 1;

            w = -std::sqrt(T(1) - z * z); // do not modify, see function comment
            if (m_signbit) {
                w = -w;
            }
        }
    }

    void operator()(std::tuple<T (&)[2]> res) const {
        tuple_assign(res, {{T(1) / std::sqrt(T(2)), std::sqrt(T(3)) * w / T(2)}});
    }

    void operator()(std::tuple<T (&)[2], T (&)[2]> res) const {
        tuple_assign(
            res,
            {{T(1) / std::sqrt(T(2)), std::sqrt(T(3)) * w / T(2)}, {0, -type_sign * std::sqrt(T(3)) * z / (T(2) * w)}}
        );
    }

    void operator()(std::tuple<T (&)[2], T (&)[2], T (&)[2]> res) const {
        tuple_assign(
            res, {{T(1) / std::sqrt(T(2)), std::sqrt(T(3)) * w / T(2)},
                  {0, -type_sign * std::sqrt(T(3)) * z / (T(2) * w)},
                  {0, std::sqrt(T(3)) * T(1) / ((z * z - T(1)) * w) / (T(2))}}
        );
    }
};

template <typename T, typename NormPolicy>
struct assoc_legendre_p_recurrence_m_m_abs;

template <typename T>
struct assoc_legendre_p_recurrence_m_m_abs<T, assoc_legendre_unnorm_policy> {
    bool m_signbit;
    int type;
    T z;
    T type_sign;

    assoc_legendre_p_recurrence_m_m_abs(bool m_signbit, int type, T z) : m_signbit(m_signbit), type(type), z(z) {
        if (type == 3) {
            type_sign = -1;
        } else {
            type_sign = 1;
        }
    }

    // other square roots can be avoided if each iteration increments by 2

    void operator()(int n, std::tuple<T (&)[2]> res) const {
        T fac;
        if (m_signbit) {
            fac = type_sign / T((2 * n) * (2 * n - 2));
        } else {
            fac = type_sign * T((2 * n - 1) * (2 * n - 3));
        }

        tuple_assign(res, {{fac * (T(1) - z * z), 0}});
    }

    void operator()(int n, std::tuple<T (&)[2], T (&)[2]> res) const {
        T fac;
        if (m_signbit) {
            fac = type_sign / T((2 * n) * (2 * n - 2));
        } else {
            fac = type_sign * T((2 * n - 1) * (2 * n - 3));
        }

        tuple_assign(res, {{fac * (T(1) - z * z), 0}, {-T(2) * fac * z, 0}});
    }

    void operator()(int n, std::tuple<T (&)[2], T (&)[2], T (&)[2]> res) const {
        T fac;
        if (m_signbit) {
            fac = type_sign / T((2 * n) * (2 * n - 2));
        } else {
            fac = type_sign * T((2 * n - 1) * (2 * n - 3));
        }

        tuple_assign(res, {{fac * (T(1) - z * z), 0}, {-T(2) * fac * z, 0}, {-T(2) * fac, 0}});
    }
};

template <typename T>
struct assoc_legendre_p_recurrence_m_m_abs<T, assoc_legendre_norm_policy> {
    bool m_signbit;
    int type;
    T z;
    T type_sign;

    assoc_legendre_p_recurrence_m_m_abs(bool m_signbit, int type, T z) : m_signbit(m_signbit), type(type), z(z) {
        if (type == 3) {
            type_sign = -1;
        } else {
            type_sign = 1;
        }
    }

    void operator()(int n, std::tuple<T (&)[2]> res) const {
        T fac = type_sign * std::sqrt(T((2 * n + 1) * (2 * n - 1)) / T(4 * n * (n - 1)));

        tuple_assign(res, {{fac * (T(1) - z * z), 0}});
    }

    void operator()(int n, std::tuple<T (&)[2], T (&)[2]> res) const {
        T fac = type_sign * std::sqrt(T((2 * n + 1) * (2 * n - 1)) / T(4 * n * (n - 1)));

        tuple_assign(res, {{fac * (T(1) - z * z), 0}, {-T(2) * fac * z, 0}});
    }

    void operator()(int n, std::tuple<T (&)[2], T (&)[2], T (&)[2]> res) const {
        T fac = type_sign * std::sqrt(T((2 * n + 1) * (2 * n - 1)) / T(4 * n * (n - 1)));

        tuple_assign(res, {{fac * (T(1) - z * z), 0}, {-T(2) * fac * z, 0}, {-T(2) * fac, 0}});
    }
};

template <typename NormPolicy, typename T, typename... OutputVals, typename Func>
void assoc_legendre_p_for_each_m_m_abs(
    NormPolicy norm, int m, int type, T z, std::tuple<OutputVals (&)[2]...> res, Func f
) {
    int m_abs = std::abs(m);
    bool m_signbit = std::signbit(m);

    assoc_legendre_p_initializer_m_m_abs<T, NormPolicy> init_m_m_abs{m_signbit, type, z};
    init_m_m_abs(res);

    assoc_legendre_p_recurrence_m_m_abs<T, NormPolicy> re_m_m_abs{m_signbit, type, z};
    forward_recur(0, m_abs + 1, re_m_m_abs, res, f);
}

template <typename T, typename NormPolicy>
struct assoc_legendre_p_recurrence_n;

template <typename T>
struct assoc_legendre_p_recurrence_n<T, assoc_legendre_unnorm_policy> {
    int m;
    int type;
    T z;

    void operator()(int n, std::tuple<T (&)[2]> res) const {
        T fac0 = -T(n + m - 1) / T(n - m);
        T fac1 = T(2 * n - 1) / T(n - m);

        tuple_assign(res, {{fac0, fac1 * z}});
    }

    void operator()(int n, std::tuple<T (&)[2], T (&)[2]> res) const {
        T fac0 = -T(n + m - 1) / T(n - m);
        T fac1 = T(2 * n - 1) / T(n - m);

        tuple_assign(res, {{fac0, fac1 * z}, {0, fac1}});
    }

    void operator()(int n, std::tuple<T (&)[2], T (&)[2], T (&)[2]> res) const {
        T fac0 = -T(n + m - 1) / T(n - m);
        T fac1 = T(2 * n - 1) / T(n - m);

        tuple_assign(res, {{fac0, fac1 * z}, {0, fac1}, {0, 0}});
    }
};

template <typename T>
struct assoc_legendre_p_recurrence_n<T, assoc_legendre_norm_policy> {
    int m;
    int type;
    T z;

    void operator()(int n, std::tuple<T (&)[2]> res) const {
        T fac0 = -std::sqrt(T((2 * n + 1) * ((n - 1) * (n - 1) - m * m)) / T((2 * n - 3) * (n * n - m * m)));
        T fac1 = std::sqrt(T((2 * n + 1) * (4 * (n - 1) * (n - 1) - 1)) / T((2 * n - 3) * (n * n - m * m)));

        tuple_assign(res, {{fac0, fac1 * z}});
    }

    void operator()(int n, std::tuple<T (&)[2], T (&)[2]> res) const {
        T fac0 = -std::sqrt(T((2 * n + 1) * ((n - 1) * (n - 1) - m * m)) / T((2 * n - 3) * (n * n - m * m)));
        T fac1 = std::sqrt(T((2 * n + 1) * (4 * (n - 1) * (n - 1) - 1)) / T((2 * n - 3) * (n * n - m * m)));

        tuple_assign(res, {{fac0, fac1 * z}, {0, fac1}});
    }

    void operator()(int n, std::tuple<T (&)[2], T (&)[2], T (&)[2]> res) const {
        T fac0 = -std::sqrt(T((2 * n + 1) * ((n - 1) * (n - 1) - m * m)) / T((2 * n - 3) * (n * n - m * m)));
        T fac1 = std::sqrt(T((2 * n + 1) * (4 * (n - 1) * (n - 1) - 1)) / T((2 * n - 3) * (n * n - m * m)));

        tuple_assign(res, {{fac0, fac1 * z}, {0, fac1}, {0, 0}});
    }
};

/**
 * Compute the associated Legendre polynomial of degree n and order n.
 *
 * We need to be careful with complex arithmetic, in particular the square roots
 * should not be modified. This is because the sign bit of a real or imaginary part,
 * even if it is equal to zero, can affect the branch cut.
 */

template <typename T, typename NormPolicy>
struct assoc_legendre_p_initializer_n;

template <typename T>
struct assoc_legendre_p_initializer_n<T, assoc_legendre_unnorm_policy> {
    int m;
    int type;
    T z;
    bool diag;

    void operator()(std::tuple<T (&)[2]> res) const {
        auto &[res0] = res;

        if (diag) {
            assoc_legendre_p_for_each_m_m_abs(assoc_legendre_unnorm, m, type, z, res, [](int n, std::tuple<T(&)[2]> p) {
            });

            res0[0] = res0[1];
        }

        int m_abs = std::abs(m);

        res0[1] = T(2 * (m_abs + 1) - 1) * z * res0[0] / T(m_abs + 1 - m);
    }

    void operator()(std::tuple<T (&)[2], T (&)[2]> res) const {
        auto &[res0, res1] = res;

        if (diag) {
            assoc_legendre_p_for_each_m_m_abs(
                assoc_legendre_unnorm, m, type, z, res, [](int n, std::tuple<T(&)[2], T(&)[2]> p) {}
            );

            res0[0] = res0[1];
            res1[0] = res1[1];
        }

        int m_abs = std::abs(m);

        res0[1] = T(2 * (m_abs + 1) - 1) * z * res0[0] / T(m_abs + 1 - m);
        res1[1] = T(2 * (m_abs + 1) - 1) * (res0[0] + z * res1[0]) / T(m_abs + 1 - m);
    }

    void operator()(std::tuple<T (&)[2], T (&)[2], T (&)[2]> res) const {
        auto &[res0, res1, res2] = res;

        if (diag) {
            assoc_legendre_p_for_each_m_m_abs(
                assoc_legendre_unnorm, m, type, z, res, [](int n, std::tuple<T(&)[2], T(&)[2], T(&)[2]> p) {}
            );

            res0[0] = res0[1];
            res1[0] = res1[1];
            res2[0] = res2[1];
        }

        int m_abs = std::abs(m);

        res0[1] = T(2 * (m_abs + 1) - 1) * z * res0[0] / T(m_abs + 1 - m);
        res1[1] = T(2 * (m_abs + 1) - 1) * (res0[0] + z * res1[0]) / T(m_abs + 1 - m);
        res2[1] = T(2 * (m_abs + 1) - 1) * (T(2) * res1[0] + z * res2[0]) / T(m_abs + 1 - m);
    }
};

template <typename T>
struct assoc_legendre_p_initializer_n<T, assoc_legendre_norm_policy> {
    int m;
    int type;
    T z;
    bool diag;

    void operator()(std::tuple<T (&)[2]> res) const {
        auto &[res0] = res;
        if (diag) {
            assoc_legendre_p_for_each_m_m_abs(assoc_legendre_unnorm, m, type, z, res, [](int n, auto) {});
            tuple_access_each(res, 0) = tuple_access_each(res, 1);
        }

        int m_abs = std::abs(m);
        res0[1] = std::sqrt(T(2 * m_abs + 3)) * z * res0[0];
    }

    void operator()(std::tuple<T (&)[2], T (&)[2]> res) const {
        auto &[res0, res1] = res;
        if (diag) {
            assoc_legendre_p_for_each_m_m_abs(assoc_legendre_unnorm, m, type, z, res, [](int n, auto) {});
            tuple_access_each(res, 0) = tuple_access_each(res, 1);
        }

        int m_abs = std::abs(m);
        res0[1] = std::sqrt(T(2 * m_abs + 3)) * z * res0[0];
        res1[1] = std::sqrt(T(2 * m_abs + 3)) * (res0[0] + z * res1[0]);
    }

    void operator()(std::tuple<T (&)[2], T (&)[2], T (&)[2]> res) const {
        auto &[res0, res1, res2] = res;
        if (diag) {
            assoc_legendre_p_for_each_m_m_abs(assoc_legendre_unnorm, m, type, z, res, [](int n, auto) {});
            tuple_access_each(res, 0) = tuple_access_each(res, 1);
        }

        int m_abs = std::abs(m);

        tuple_assign(
            tuple_access_each(res, 1),
            {std::sqrt(T(2 * m_abs + 3)) * z * res0[0], std::sqrt(T(2 * m_abs + 3)) * (res0[0] + z * res1[0]),
             std::sqrt(T(2 * m_abs + 3)) * (T(2) * res1[0] + z * res2[0])}
        );
    }
};

template <typename T>
T assoc_legendre_p_pm1(assoc_legendre_unnorm_policy norm, int n, int m, int type, T z) {
    if (m == 0) {
        return 1;
    }

    return 0;
}

template <typename T>
T assoc_legendre_p_pm1(assoc_legendre_norm_policy norm, int n, int m, int type, T z) {
    if (m == 0) {
        return 1;
    }

    return 0;
}

template <typename NormPolicy, typename T>
void assoc_legendre_p_pm1(NormPolicy norm, int n, int m, int type, T z, std::tuple<T &> res) {
    std::get<0>(res) = assoc_legendre_p_pm1(norm, n, m, type, z);
}

template <typename NormPolicy, typename T>
void assoc_legendre_p_pm1(NormPolicy norm, int n, int m, int type, T z, std::tuple<T &, T &> res_) {
    auto &[res, res_jac] = res_;
    res = assoc_legendre_p_pm1(norm, n, m, type, z);

    T type_sign;
    if (type == 3) {
        type_sign = -1;
    } else {
        type_sign = 1;
    }

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

template <typename NormPolicy, typename T>
void assoc_legendre_p_pm1(NormPolicy norm, int n, int m, int type, T z, std::tuple<T &, T &, T &> res_) {
    auto &[res, res_jac, res_hess] = res_;
    assoc_legendre_p_pm1(norm, n, m, type, z, std::tie(res, res_jac));

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
template <typename NormPolicy, typename T, typename... OutputVals, typename Func>
void assoc_legendre_p_for_each_n(
    NormPolicy norm, int n, int m, int type, T z, bool recur_m_m_abs, std::tuple<OutputVals (&)[2]...> res, Func f
) {
    std::tuple<OutputVals...> res_m_m_abs;
    if (!recur_m_m_abs) {
        res_m_m_abs = tuple_access_each(res, 0);
    }

    tuple_fill_each(res, 0);

    int m_abs = std::abs(m);
    if (m_abs > n) {
        for (int j = 0; j <= n; ++j) {
            f(j, res);
        }
    } else {
        for (int j = 0; j < m_abs; ++j) {
            f(j, res);
        }

        if (std::abs(std::real(z)) == 1 && std::imag(z) == 0) {
            for (int j = m_abs; j <= n; ++j) {
                forward_recur_shift_left(res);
                assoc_legendre_p_pm1(norm, j, m, type, z, tuple_access_each(res, 1));

                f(j, res);
            }
        } else {
            if (!recur_m_m_abs) {
                tuple_access_each(res, 0) = res_m_m_abs;
            }

            assoc_legendre_p_initializer_n<T, NormPolicy> init_n{m, type, z, recur_m_m_abs};
            init_n(res);

            assoc_legendre_p_recurrence_n<T, NormPolicy> re_n{m, type, z};
            forward_recur(m_abs, n + 1, re_n, res, f);
        }
    }
}

template <typename NormPolicy, typename T, typename... OutputVals, typename Func>
void assoc_legendre_p_for_each_n_m(
    NormPolicy norm, int n, int m, int type, T z, std::tuple<OutputVals (&)[2]...> res, Func f
) {
    std::tuple<OutputVals[2]...> res_m_m_abs;
    assoc_legendre_p_for_each_m_m_abs(
        norm, m, type, z, tuple_ref_each(res_m_m_abs),
        [norm, n, type, z, &res, f](int i, std::tuple<OutputVals(&)[2]...> res_m_m_abs) {
            tuple_access_each(res, 0) = tuple_access_each(res_m_m_abs, 1);

            assoc_legendre_p_for_each_n(
                norm, n, i, type, z, false, res,
                [f, i](int j, std::tuple<OutputVals(&)[2]...> res_n) { f(j, i, res_n); }
            );
        }
    );

    assoc_legendre_p_for_each_m_m_abs(
        norm, -m, type, z, tuple_ref_each(res_m_m_abs),
        [norm, n, type, z, &res, f](int i_abs, std::tuple<OutputVals(&)[2]...> res_m_m_abs) {
            tuple_access_each(res, 0) = tuple_access_each(res_m_m_abs, 1);

            int i = -i_abs;
            assoc_legendre_p_for_each_n(
                norm, n, i, type, z, false, res,
                [f, i](int j, std::tuple<OutputVals(&)[2]...> res_n) { f(j, i, res_n); }
            );
        }
    );
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
template <typename NormPolicy, typename T, typename... OutputVals>
void multi_assoc_legendre_p(NormPolicy norm, int n, int m, int type, T z, std::tuple<OutputVals &...> res) {
    std::tuple<OutputVals[2]...> res_n;
    assoc_legendre_p_for_each_n(
        norm, n, m, type, z, true, tuple_ref_each(res_n), [](int n, std::tuple<OutputVals(&)[2]...> res_n) {}
    );

    res = tuple_access_each(res_n, 1);
}

template <typename NormPolicy, typename T>
T multi_assoc_legendre_p(NormPolicy norm, int n, int m, int type, T z) {
    T res;
    multi_assoc_legendre_p(norm, n, m, type, z, std::tie(res));

    return res;
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
template <typename NormPolicy, typename T, typename... OutputMats>
void multi_assoc_legendre_p_all(NormPolicy norm, int type, T z, std::tuple<OutputMats &...> res) {
    auto &res0 = std::get<0>(res);
    int n = res0.extent(0) - 1;
    int m = (res0.extent(1) - 1) / 2;

    std::tuple<typename OutputMats::value_type[2]...> p;
    assoc_legendre_p_for_each_n_m(norm, n, m, type, z, tuple_ref_each(p), [&res](int n, int m, auto res_n_m) {
        if (m >= 0) {
            tuple_call_each(res, n, m) = tuple_access_each(res_n_m, 1);
        } else {
            auto &res0 = std::get<0>(res);
            tuple_call_each(res, n, m + res0.extent(1)) = tuple_access_each(res_n_m, 1);
        }
    });
}

template <typename NormPolicy, typename T, typename... OutputVals>
void assoc_legendre_p(NormPolicy norm, int n, int m, T z, std::tuple<OutputVals &...> res) {
    int type;
    if (std::abs(z) <= 1) {
        type = 2;
    } else {
        type = 3;
    }

    multi_assoc_legendre_p(norm, n, m, type, z, res);
}

template <typename NormPolicy, typename T>
T assoc_legendre_p(NormPolicy norm, int n, int m, T z) {
    T res;
    assoc_legendre_p(norm, n, m, z, std::tie(res));

    return res;
}

template <typename NormPolicy, typename T, typename... OutputMats>
void assoc_legendre_p_all(NormPolicy norm, T z, std::tuple<OutputMats &...> res) {
    int type;
    if (std::abs(z) <= 1) {
        type = 2;
    } else {
        type = 3;
    }

    multi_assoc_legendre_p_all(norm, type, z, res);
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
