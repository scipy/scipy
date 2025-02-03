#pragma once

#include "dual.h"
#include "error.h"
#include "recur.h"

namespace xsf {

template <typename T>
struct legendre_p_initializer_n {
    T z;

    void operator()(T (&res)[2]) const {
        res[0] = T(1);
        res[1] = z;
    }
};

template <typename T>
struct legendre_p_recurrence_n {
    T z;

    void operator()(int n, T (&res)[2]) const {
        using value_type = remove_dual_t<T>;
        value_type fac0 = -value_type(n - 1) / value_type(n);
        value_type fac1 = value_type(2 * n - 1) / value_type(n);

        res[0] = fac0;
        res[1] = fac1 * z;
    }
};

/**
 * Compute the Legendre polynomial of degree n.
 *
 * @param n degree of the polynomial
 * @param z argument of the polynomial, either real or complex
 * @param f a function to be called as callback(j, res) for 0 <= j <= n
 * @param res value and derivatives of the polynomial
 */
template <typename T, typename Func>
void legendre_p_for_each_n(int n, T z, T (&res)[2], Func f) {
    legendre_p_initializer_n<T> init_n{z};
    init_n(res);

    legendre_p_recurrence_n<T> re_n{z};
    forward_recur(0, n + 1, re_n, res, f);
}

/*
 * Compute the Legendre polynomial of degree n.
 *
 * @param n degree of the polynomial
 * @param z argument of the polynomial, either real or complex
 *
 * @return value of the polynomial
 */
template <typename T>
T legendre_p(int n, T z) {
    T res_n[2];
    legendre_p_for_each_n(n, z, res_n, [](int n, const T(&res_n)[2]) {});

    return res_n[1];
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

    T res_n[2];
    legendre_p_for_each_n(n, z, res_n, [&res](int n, const T(&res_n)[2]) { res(n) = res_n[1]; });
}

struct assoc_legendre_unnorm_policy {};

struct assoc_legendre_norm_policy {};

constexpr assoc_legendre_unnorm_policy assoc_legendre_unnorm;

constexpr assoc_legendre_norm_policy assoc_legendre_norm;

template <typename T, typename NormPolicy>
struct assoc_legendre_p_initializer_m_abs_m;

template <typename T>
struct assoc_legendre_p_initializer_m_abs_m<T, assoc_legendre_unnorm_policy> {
    bool m_signbit;
    T z;
    int branch_cut;
    T w;

    assoc_legendre_p_initializer_m_abs_m(bool m_signbit, T z, int branch_cut)
        : m_signbit(m_signbit), z(z), branch_cut(branch_cut) {
        if (branch_cut == 3) {
            w = sqrt(z - T(1)) * sqrt(z + T(1)); // form for this branch cut
        } else {
            w = -sqrt(T(1) - z * z); // form for this branch cut
            if (m_signbit) {
                w = -w;
            }
        }
    }

    void operator()(T (&res)[2]) const {
        res[0] = T(1);
        res[1] = w;

        if (m_signbit) {
            res[1] /= 2;
        }
    }
};

template <typename T>
struct assoc_legendre_p_initializer_m_abs_m<T, assoc_legendre_norm_policy> {
    bool m_signbit;
    T z;
    int branch_cut;
    T w;

    assoc_legendre_p_initializer_m_abs_m(bool m_signbit, T z, int branch_cut)
        : m_signbit(m_signbit), z(z), branch_cut(branch_cut) {
        if (branch_cut == 3) {

            w = sqrt(z - T(1)) * sqrt(z + T(1)); // form for this branch cut
        } else {
            w = -sqrt(T(1) - z * z); // form for this branch cut
            if (m_signbit) {
                w = -w;
            }
        }
    }

    void operator()(T (&res)[2]) const {
        res[0] = T(1) / sqrt(T(2));
        res[1] = sqrt(T(3)) * w / T(2);
    }
};

template <typename T, typename NormPolicy>
struct assoc_legendre_p_recurrence_m_abs_m;

template <typename T>
struct assoc_legendre_p_recurrence_m_abs_m<T, assoc_legendre_unnorm_policy> {
    T z;
    int branch_cut;
    T branch_cut_sign;

    assoc_legendre_p_recurrence_m_abs_m(T z, int branch_cut) : z(z), branch_cut(branch_cut) {
        if (branch_cut == 3) {
            branch_cut_sign = T(-1);
        } else {
            branch_cut_sign = T(1);
        }
    }

    // other square roots can be avoided if each iteration increments by 2
    void operator()(int m, T (&res)[2]) const {
        int m_abs = abs(m);

        T fac;
        if (m < 0) {
            fac = branch_cut_sign / T((2 * m_abs) * (2 * m_abs - 2));
        } else {
            fac = branch_cut_sign * T((2 * m_abs - 1) * (2 * m_abs - 3));
        }

        res[0] = fac * (T(1) - z * z);
        res[1] = T(0);
    }
};

template <typename T>
struct assoc_legendre_p_recurrence_m_abs_m<T, assoc_legendre_norm_policy> {
    T z;
    int branch_cut;
    T branch_cut_sign;

    assoc_legendre_p_recurrence_m_abs_m(T z, int branch_cut) : z(z), branch_cut(branch_cut) {
        if (branch_cut == 3) {
            branch_cut_sign = T(-1);
        } else {
            branch_cut_sign = T(1);
        }
    }

    void operator()(int m, T (&res)[2]) const {
        int m_abs = abs(m);

        T fac = branch_cut_sign * sqrt(T((2 * m_abs + 1) * (2 * m_abs - 1)) / T(4 * m_abs * (m_abs - 1)));

        res[0] = fac * (T(1) - z * z);
        res[1] = T(0);
    }
};

template <typename NormPolicy, typename T, typename Func>
void assoc_legendre_p_for_each_m_abs_m(NormPolicy norm, int m, T z, int branch_cut, T (&res)[2], Func f) {
    bool m_signbit;
    if (m < 0) {
        m_signbit = true;
    } else {
        m_signbit = false;
    }

    assoc_legendre_p_initializer_m_abs_m<T, NormPolicy> init_m_abs_m{m_signbit, z, branch_cut};
    init_m_abs_m(res);

    assoc_legendre_p_recurrence_m_abs_m<T, NormPolicy> re_m_abs_m{z, branch_cut};
    if (m >= 0) {
        forward_recur(0, m + 1, re_m_abs_m, res, f);
    } else {
        backward_recur(0, m - 1, re_m_abs_m, res, f);
    }
}

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
    T z;
    int type;

    void operator()(const T &res_m_abs_m, T (&res)[2]) const {
        int m_abs = abs(m);
        T fac = T(2 * (m_abs + 1) - 1) / T(m_abs + 1 - m);

        res[0] = res_m_abs_m;
        res[1] = fac * z * res_m_abs_m;
    }
};

template <typename T>
struct assoc_legendre_p_initializer_n<T, assoc_legendre_norm_policy> {
    int m;
    T z;
    int type;

    void operator()(const T &res_m_abs_m, T (&res)[2]) const {
        T fac = sqrt(T(2 * abs(m) + 3));

        res[0] = res_m_abs_m;
        res[1] = fac * z * res_m_abs_m;
    }
};

template <typename T, typename NormPolicy>
struct assoc_legendre_p_recurrence_n;

template <typename T>
struct assoc_legendre_p_recurrence_n<T, assoc_legendre_unnorm_policy> {
    int m;
    T z;
    int type;

    void operator()(int n, T (&res)[2]) const {
        using value_type = remove_dual_t<T>;
        value_type fac0 = -value_type(n + m - 1) / value_type(n - m);
        value_type fac1 = value_type(2 * n - 1) / value_type(n - m);

        res[0] = fac0;
        res[1] = fac1 * z;
    }
};

template <typename T>
struct assoc_legendre_p_recurrence_n<T, assoc_legendre_norm_policy> {
    int m;
    T z;
    int type;

    void operator()(int n, T (&res)[2]) const {
        using value_type = remove_dual_t<T>;
        value_type fac0 =
            -sqrt(value_type((2 * n + 1) * ((n - 1) * (n - 1) - m * m)) / value_type((2 * n - 3) * (n * n - m * m)));
        value_type fac1 =
            sqrt(value_type((2 * n + 1) * (4 * (n - 1) * (n - 1) - 1)) / value_type((2 * n - 3) * (n * n - m * m)));

        res[0] = fac0;
        res[1] = fac1 * z;
    }
};

template <typename NormPolicy, typename T>
void assoc_legendre_p_pm1(NormPolicy norm, int n, int m, T z, int branch_cut, T &res) {
    if (m == 0) {
        res = T(1);
    } else {
        res = T(0);
    }
}

template <typename NormPolicy, typename T, size_t Order>
void assoc_legendre_p_pm1(NormPolicy norm, int n, int m, dual<T, Order> z, int branch_cut, dual<T, Order> &res) {
    if (m == 0) {
        res[0] = T(1);
    } else {
        res[0] = T(0);
    }

    T branch_cut_sign;
    if (branch_cut == 3) {
        branch_cut_sign = -1;
    } else {
        branch_cut_sign = 1;
    }

    if (Order >= 1) {
        if (abs(m) > n) {
            res[1] = 0;
        } else if (m == 0) {
            res[1] = T(n) * T(n + 1) * std::pow(z[0], T(n + 1)) / T(2);
        } else if (m == 1) {
            res[1] = std::pow(z[0], T(n)) * std::numeric_limits<remove_complex_t<T>>::infinity();
        } else if (m == 2) {
            res[1] = -branch_cut_sign * T(n + 2) * T(n + 1) * T(n) * T(n - 1) * std::pow(z[0], T(n + 1)) / T(4);
        } else if (m == -2) {
            res[1] = -branch_cut_sign * std::pow(z[0], T(n + 1)) / T(4);
        } else if (m == -1) {
            res[1] = -std::pow(z[0], T(n)) * std::numeric_limits<remove_complex_t<T>>::infinity();
        } else {
            res[1] = 0;
        }

        if (Order >= 2) {
            if (abs(m) > n) {
                res[2] = 0;
            } else if (m == 0) {
                res[2] = T(n + 2) * T(n + 1) * T(n) * T(n - 1) / T(8);
            } else if (m == 1) {
                res[2] = std::numeric_limits<remove_complex_t<T>>::infinity();
            } else if (m == 2) {
                res[2] = -T((n + 1) * n - 3) * T(n + 2) * T(n + 1) * T(n) * T(n - 1) / T(12);
            } else if (m == 3) {
                res[2] = std::numeric_limits<remove_complex_t<T>>::infinity();
            } else if (m == 4) {
                res[2] = T(n + 4) * T(n + 3) * T(n + 2) * T(n + 1) * T(n) * T(n - 1) * T(n - 2) * T(n - 3) / T(48);
            } else if (m == -4) {
                res[2] = 0;
            } else if (m == -3) {
                res[2] = -std::numeric_limits<remove_complex_t<T>>::infinity();
            } else if (m == -2) {
                res[2] = -T(1) / T(4);
            } else if (m == -1) {
                res[2] = -std::numeric_limits<remove_complex_t<T>>::infinity();
            } else {
                res[2] = 0;
            }
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
 * @param callback a function to be called as callback(j, m, type, z, p, p_prev, args...) for 0 <= j <= n
 * @param args arguments to forward to the callback
 *
 * @return value of the polynomial
 */
template <typename NormPolicy, typename T, typename Func>
void assoc_legendre_p_for_each_n(
    NormPolicy norm, int n, int m, T z, int branch_cut, const T &res_m_abs_m, T (&res)[2], Func f
) {
    res[0] = 0;
    res[1] = 0;

    int m_abs = abs(m);
    if (m_abs > n) {
        for (int j = 0; j <= n; ++j) {
            f(j, res);
        }
    } else {
        for (int j = 0; j < m_abs; ++j) {
            f(j, res);
        }

        if (abs(real(z)) == 1 && imag(z) == 0) {
            for (int j = m_abs; j <= n; ++j) {
                forward_recur_shift_left(res);
                assoc_legendre_p_pm1(norm, j, m, z, branch_cut, res[1]);

                f(j, res);
            }
        } else {
            assoc_legendre_p_initializer_n<T, NormPolicy> init_n{m, z, branch_cut};
            init_n(res_m_abs_m, res);

            assoc_legendre_p_recurrence_n<T, NormPolicy> re_n{m, z, branch_cut};
            forward_recur(m_abs, n + 1, re_n, res, f);
        }
    }
}

template <typename NormPolicy, typename T, typename Func>
void assoc_legendre_p_for_each_n(NormPolicy norm, int n, int m, T z, int branch_cut, T (&res)[2], Func f) {
    assoc_legendre_p_for_each_m_abs_m(norm, m, z, branch_cut, res, [](int m, const T(&res)[2]) {});

    T res_m_abs_m = res[1];
    assoc_legendre_p_for_each_n(norm, n, m, z, branch_cut, res_m_abs_m, res, f);
}

template <typename NormPolicy, typename T, typename Func>
void assoc_legendre_p_for_each_n_m(NormPolicy norm, int n, int m, T z, int branch_cut, T (&res)[2], Func f) {
    T res_m_abs_m[2];
    assoc_legendre_p_for_each_m_abs_m(
        norm, m, z, branch_cut, res_m_abs_m,
        [norm, n, z, branch_cut, &res, f](int m, const T(&res_m_abs_m)[2]) {
            res[0] = res_m_abs_m[1];

            assoc_legendre_p_for_each_n(
                norm, n, m, z, branch_cut, res_m_abs_m[1], res, [f, m](int n, const T(&res_n)[2]) { f(n, m, res_n); }
            );
        }
    );
    assoc_legendre_p_for_each_m_abs_m(
        norm, -m, z, branch_cut, res_m_abs_m,
        [norm, n, z, branch_cut, &res, f](int m, const T(&res_m_abs_m)[2]) {
            res[0] = res_m_abs_m[1];

            assoc_legendre_p_for_each_n(
                norm, n, m, z, branch_cut, res_m_abs_m[1], res, [f, m](int n, const T(&res_n)[2]) { f(n, m, res_n); }
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
template <typename NormPolicy, typename T>
T assoc_legendre_p(NormPolicy norm, int n, int m, T z, int branch_cut) {
    T res_n[2];
    assoc_legendre_p_for_each_n(norm, n, m, z, branch_cut, res_n, [](int n, const T(&res_n)[2]) {});

    return res_n[1];
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
template <typename NormPolicy, typename T, typename OutputMat>
void assoc_legendre_p_all(NormPolicy norm, T z, int branch_cut, OutputMat res) {
    int n = res.extent(0) - 1;
    int m = (res.extent(1) - 1) / 2;

    T p[2];
    assoc_legendre_p_for_each_n_m(norm, n, m, z, branch_cut, p, [&res](int n, int m, const T(&res_n_m)[2]) {
        if (m >= 0) {
            res(n, m) = res_n_m[1];
        } else {
            res(n, m + res.extent(1)) = res_n_m[1];
        }
    });
}

template <typename T>
struct sph_legendre_p_initializer_m_abs_m {
    bool m_signbit;
    T theta;
    T theta_sin;

    sph_legendre_p_initializer_m_abs_m(bool m_signbit, T theta)
        : m_signbit(m_signbit), theta(theta), theta_sin(sin(theta)) {}

    void operator()(T (&res)[2]) const {
        T fac0 = T(1) / (T(2) * sqrt(T(M_PI)));
        T fac1 = -sqrt(T(3)) / (T(2) * sqrt(T(2) * T(M_PI)));
        if (m_signbit) {
            fac1 = -fac1;
        }

        res[0] = fac0;
        res[1] = fac1 * abs(theta_sin);
    }
};

template <typename T>
struct sph_legendre_p_recurrence_m_abs_m {
    T theta;
    T theta_sin;

    sph_legendre_p_recurrence_m_abs_m(T theta) : theta(theta), theta_sin(sin(theta)) {}

    void operator()(int m, T (&res)[2]) const {
        int m_abs = abs(m);

        T fac = sqrt(T((2 * m_abs + 1) * (2 * m_abs - 1)) / T(4 * m_abs * (m_abs - 1)));

        res[0] = fac * theta_sin * theta_sin;
        res[1] = 0;
    }
};

template <typename T, typename Func>
void sph_legendre_p_for_each_m_abs_m(int m, T theta, T (&res)[2], Func f) {
    bool m_signbit;
    if (m < 0) {
        m_signbit = true;
    } else {
        m_signbit = false;
    }

    sph_legendre_p_initializer_m_abs_m<T> init_m_abs_m{m_signbit, theta};
    init_m_abs_m(res);

    sph_legendre_p_recurrence_m_abs_m<T> re_m_abs_m{theta};
    if (m >= 0) {
        forward_recur(0, m + 1, re_m_abs_m, res, f);
    } else {
        backward_recur(0, m - 1, re_m_abs_m, res, f);
    }
}

template <typename T>
struct sph_legendre_p_initializer_n {
    int m;
    T theta;
    T theta_cos;

    sph_legendre_p_initializer_n(int m, T theta) : m(m), theta(theta), theta_cos(cos(theta)) {}

    void operator()(const T &res_m_abs_m, T (&res)[2]) const {
        T fac = sqrt(T(2 * abs(m) + 3));

        res[0] = res_m_abs_m;
        res[1] = fac * theta_cos * res_m_abs_m;
    }
};

template <typename T>
struct sph_legendre_p_recurrence_n {
    int m;
    T theta;
    T theta_cos;

    sph_legendre_p_recurrence_n(int m, T theta) : m(m), theta(theta), theta_cos(cos(theta)) {}

    void operator()(int n, T (&res)[2]) const {
        using value_type = remove_dual_t<T>;
        value_type fac0 =
            -sqrt(value_type((2 * n + 1) * ((n - 1) * (n - 1) - m * m)) / value_type((2 * n - 3) * (n * n - m * m)));
        value_type fac1 =
            sqrt(value_type((2 * n + 1) * (4 * (n - 1) * (n - 1) - 1)) / value_type((2 * n - 3) * (n * n - m * m)));

        res[0] = fac0;
        res[1] = fac1 * theta_cos;
    }
};

/**
 * Compute the spherical Legendre polynomial of degree n and order m.
 *
 * @param n degree of the polynomial
 * @param m order of the polynomial
 * @param theta z = cos(theta) argument of the polynomial, either real or complex
 * @param callback a function to be called as callback(j, m, type, z, p, p_prev, args...) for 0 <= j <= n
 * @param args arguments to forward to the callback
 *
 * @return value of the polynomial
 */
template <typename T, typename Func>
void sph_legendre_p_for_each_n(int n, int m, T theta, const T &res_m_abs_m, T (&res)[2], Func f) {
    res[0] = 0;
    res[1] = 0;

    int m_abs = abs(m);
    if (m_abs > n) {
        for (int j = 0; j <= n; ++j) {
            f(j, res);
        }
    } else {
        for (int j = 0; j < m_abs; ++j) {
            f(j, res);
        }

        sph_legendre_p_initializer_n<T> init_n{m, theta};
        init_n(res_m_abs_m, res);

        sph_legendre_p_recurrence_n<T> re_n{m, theta};
        forward_recur(m_abs, n + 1, re_n, res, f);
    }
}

template <typename T, typename Func>
void sph_legendre_p_for_each_n(int n, int m, T theta, T (&res)[2], Func f) {
    sph_legendre_p_for_each_m_abs_m(m, theta, res, [](int m, auto) {});

    T res_m_abs_m = res[1];
    sph_legendre_p_for_each_n(n, m, theta, res_m_abs_m, res, f);
}

template <typename T, typename Func>
void sph_legendre_p_for_each_n_m(int n, int m, T theta, T (&res)[2], Func f) {
    T res_m_abs_m[2];
    sph_legendre_p_for_each_m_abs_m(m, theta, res_m_abs_m, [n, theta, &res, f](int m, const T(&res_m_abs_m)[2]) {
        res[0] = res_m_abs_m[1];

        sph_legendre_p_for_each_n(n, m, theta, res_m_abs_m[1], res, [f, m](int n, const T(&res_n)[2]) {
            f(n, m, res_n);
        });
    });
    sph_legendre_p_for_each_m_abs_m(-m, theta, res_m_abs_m, [n, theta, &res, f](int m, const T(&res_m_abs_m)[2]) {
        res[0] = res_m_abs_m[1];

        sph_legendre_p_for_each_n(n, m, theta, res_m_abs_m[1], res, [f, m](int n, const T(&res_n)[2]) {
            f(n, m, res_n);
        });
    });
}

template <typename T>
T sph_legendre_p(int n, int m, T theta) {
    T res_n[2];
    sph_legendre_p_for_each_n(n, m, theta, res_n, [](int n, const T(&res_n)[2]) {});

    return res_n[1];
}

template <typename T, typename OutputMat>
void sph_legendre_p_all(T theta, OutputMat res) {
    int n_max = res.extent(0) - 1;
    int m_max = (res.extent(1) - 1) / 2;

    T res_n_m[2];
    sph_legendre_p_for_each_n_m(n_max, m_max, theta, res_n_m, [m_max, &res](int n, int m, const T(&res_n_m)[2]) {
        if (m >= 0) {
            res(n, m) = res_n_m[1];
        } else {
            res(n, m + 2 * m_max + 1) = res_n_m[1];
        }
    });
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

    if (real(z) == 1) {
        for (int k = 0; k <= n; ++k) {
            cqn(k) = 1e300;
            cqd(k) = 1e300;
        }
        return;
    }
    int ls = ((abs(z) > 1.0) ? -1 : 1);

    cq0 = std::log(static_cast<T>(ls) * (static_cast<T>(1) + z) / (static_cast<T>(1) - z)) / static_cast<T>(2);
    cq1 = z * cq0 - static_cast<T>(1);

    cqn(0) = cq0;
    cqn(1) = cq1;

    if (abs(z) < 1.0001) {
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
        if (abs(z) > 1.1) {
            km = 40 + n;
        } else {
            km = (int) ((40 + n) * floor(-1.0 - 1.8 * log(abs(z - static_cast<T>(1)))));
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

    if ((abs(real(z)) == 1) && (imag(z) == 0)) {
        for (i = 0; i < (m + 1); i++) {
            for (j = 0; j < (n + 1); j++) {
                cqm(i, j) = 1e300;
                cqd(i, j) = 1e300;
            }
        }

        return;
    }

    T xc = abs(z);
    ls = 0;
    if ((imag(z) == 0) || (xc < 1)) {
        ls = 1;
    }
    if (xc > 1) {
        ls = -1;
    }
    zs = static_cast<T>(ls) * (static_cast<T>(1) - z * z);
    zq = sqrt(zs);

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

} // namespace xsf
