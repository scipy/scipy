#pragma once

#include <iostream>

#include "config.h"
#include "error.h"
#include "gamma.h"

namespace special {

// Translated into C++ by SciPy developers in 2024.
//
// ===============================================
// Purpose: Compute Legendre polynomials Pn(z)
//          and their derivatives Pn'(z)
// Input :  x --- Argument of Pn(z)
//          n --- Degree of Pn(z) ( n = 0,1,...)
// Output:  PN(n) --- Pn(z)
//          PD(n) --- Pn'(z)
// ===============================================

template <typename T, typename Callable>
T legendre_p(unsigned int n, T z, Callable callback) {
    T p = 1;
    callback(0, z, p);

    if (n > 0) {
        T p_prev = p;
        p = z;
        callback(1, z, p);

        for (unsigned int j = 2; j <= n; ++j) {
            T p_prev_prev = p_prev;
            p_prev = p;
            p = (T(2 * j - 1) * z * p_prev - T(j - 1) * p_prev_prev) / T(j);
            callback(j, z, p);
        }
    }

    return p;
}

template <typename T>
T legendre_p(unsigned int n, T z) {
    return legendre_p(n, z, [](unsigned int j, T z, T value) {});
}

template <typename T, typename OutputVec>
void legendre_p_all(T z, OutputVec p) {
    unsigned int n = p.extent(0) - 1;

    legendre_p(n, z, [p](unsigned int j, T z, T value) { p(j) = value; });
}

template <typename T>
T legendre_p_jac_next(unsigned int n, T z, T p, T p_prev) {
    if (n == 0) {
        return 0;
    }

    if (std::abs(std::real(z)) == 1 && std::imag(z) == 0) {
        return T(n) * T(n + 1) * std::pow(std::real(z), T(n + 1)) / T(2);
    }

    return T(n) * (p_prev - z * p) / (T(1) - z * z);
}

template <typename T, typename Callable>
void legendre_p_jac(unsigned int n, T z, Callable callback) {
    T value_jac;
    T value_prev = std::numeric_limits<T>::quiet_NaN();
    legendre_p(n, z, [&value_jac, &value_prev, &callback](unsigned int j, T z, T value) {
        value_jac = legendre_p_jac_next(j, z, value, value_prev);
        value_prev = value;
        callback(j, z, value, value_jac);
    });
}

template <typename T>
void legendre_p_jac(unsigned int n, T z, T &p, T &p_jac) {
    legendre_p_jac(n, z, [&p, &p_jac](unsigned int j, T z, T value, T value_jac) {
        p = value;
        p_jac = value_jac;
    });
}

template <typename T, typename OutputVec1, typename OutputVec2>
void legendre_p_jac_all(T z, OutputVec1 p, OutputVec2 p_jac) {
    unsigned int n = p.extent(0) - 1;

    legendre_p_jac(n, z, [p, p_jac](unsigned int j, T z, T value, T value_jac) {
        p(j) = value;
        p_jac(j) = value_jac;
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

template <typename T, typename Callable>
T assoc_legendre_p(unsigned int n, unsigned int m, T x, Callable callback) {
    unsigned int m_abs = m;

    if (m_abs > n) {
        for (unsigned int j = 0; j < n; ++j) {
            callback(j, m, x, 0);
        }

        return 0;
    }

    unsigned int ls = (std::abs(x) > 1 ? -1 : 1);
    T xq = std::sqrt(ls * (1 - x * x));
    if (x < -1) {
        xq = -xq; // Ensure connection to the complex-valued function for |x| > 1
    }

    T p = 1;
    for (unsigned int j = 0; j < m; ++j) {
        callback(j, m, x, 0);

        p *= -ls * T(2 * j + 1) * xq;
    }

    callback(m, m, x, p);

    if (m != n) {
        T p_prev = p;
        p = T(2 * m + 1) * x * p_prev;
        callback(m + 1, m, x, p);

        for (unsigned int j = m + 2; j <= n; ++j) {
            T p_prev_prev = p_prev;
            p_prev = p;
            p = (T(2 * j - 1) * x * p_prev - T(m + j - 1) * p_prev_prev) / T(j - m);
            callback(j, m, x, p);
        }
    }

    return p;
}

/*
    if (std::abs(x) == 1) {
        for (long i = 1; i <= n; i++) {
            p(0, i) = std::pow(x, i);
        }
    }
*/

template <typename T>
T assoc_legendre_p(unsigned int n, unsigned int m, T x) {
    return assoc_legendre_p(n, m, x, [](unsigned int j, unsigned int i, T x, T value) {});
}

template <typename T, typename OutputMat>
void assoc_legendre_p_all(T x, OutputMat p) {
    unsigned int m = p.extent(0) - 1;
    unsigned int n = p.extent(1) - 1;

    for (unsigned int i = 0; i <= m; ++i) {
        assoc_legendre_p(n, i, x, [p](unsigned int j, unsigned int i, T x, T value) { p(i, j) = value; });
    }
}

template <typename T, typename Callable>
T assoc_legendre_p(unsigned int n, int m, T x, Callable callback) {
    unsigned int m_abs = std::abs(m);
    bool m_signbit = std::signbit(m);

    return assoc_legendre_p(n, m_abs, x, [m_signbit, &callback](unsigned int j, unsigned int i_abs, T x, T value_unm) {
        int i = i_abs;
        T value = value_unm;
        if (m_signbit) {
            i = -i;
            if (i_abs <= j) {
                value *= gamma_ratio(j - i_abs + 1, j + i_abs + 1);
                if (std::abs(x) < 1) {
                    value *= std::pow(-1, i_abs);
                }
            }
        }

        callback(j, i, x, value);
    });
}

template <typename T, typename OutputMat>
void assoc_legendre_p_all(bool m_signbit, T x, OutputMat p) {
    unsigned int m_abs = p.extent(0) - 1;
    unsigned int n = p.extent(1) - 1;

    int m = m_abs;
    if (m_signbit) {
        m = -m;
    }

    assoc_legendre_p(n, m, x, [p](unsigned int j, int i, T x, T value) {
        unsigned int i_abs = std::abs(i);
        p(i_abs, j) = value;
    });
}

template <typename T>
T assoc_legendre_p_jac_next(unsigned int n, int m, T x, T p_curr, T p_prev) {
    if (std::abs(x) == 1) {
        if (m == 0) {
            return std::pow(x, n + 1) * n * (n + 1) / 2;
        }

        if (m == 1) {
            return std::numeric_limits<T>::infinity();
        }

        if (m == 2) {
            return -(n + 2) * (n + 1) * n * (n - 1) * std::pow(x, n + 1) / 4;
        }

        return 0;
    }

    if (m > n) {
        return 0;
    }

    if (m == n) {
        return n * x * p_curr / (x * x - 1);
    }

    int ls = (std::abs(x) > 1 ? -1 : 1);
    T xq = std::sqrt(ls * (1 - x * x));
    // Ensure connection to the complex-valued function for |x| > 1
    if (x < -1) {
        xq = -xq;
    }
    T xs = ls * (1 - x * x);

    return (n * x * p_curr - (n + m) * p_prev) / (x * x - 1);
}

template <typename T>
T assoc_legendre_p_jac_next(unsigned int n, unsigned int m, bool m_signbit, T x, T p, T p_prev) {
    if (m_signbit) {
        return assoc_legendre_p_jac_next(m, n, x, p, p_prev);
    }

    return assoc_legendre_p_jac_next(m, n, x, p, p_prev);
}

template <typename T, typename Callable>
T assoc_legendre_p_jac(unsigned int n, unsigned int m, T x, Callable callback) {
    T value_jac;
    T value_prev = std::numeric_limits<T>::quiet_NaN();
    assoc_legendre_p(n, m, x, [&value_jac, &value_prev, &callback](unsigned int j, unsigned int i, T x, T value) {
        value_jac = assoc_legendre_p_jac_next(j, i, x, value, value_prev);
        value_prev = value;
        callback(j, i, x, value_jac);
    });

    return value_jac;
}

template <typename T, typename Callable>
T assoc_legendre_p_jac(unsigned int n, int m, T x, Callable callback) {
    T value_jac;
    T value_prev = std::numeric_limits<T>::quiet_NaN();
    assoc_legendre_p(n, m, x, [&value_jac, &value_prev, &callback](unsigned int j, int i, T x, T value) {
        value_jac = assoc_legendre_p_jac_next(j, i, x, value, value_prev);
        value_prev = value;
        callback(j, i, x, value, value_jac);
    });

    return value_jac;
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

// Translated into C++ by SciPy developers in 2024.
//
// =========================================================
// Purpose: Compute the associated Legendre functions Pmn(z)
//          and their derivatives Pmn'(z) for a complex
//          argument
// Input :  x     --- Real part of z
//          y     --- Imaginary part of z
//          m     --- Order of Pmn(z),  m = 0,1,2,...,n
//          n     --- Degree of Pmn(z), n = 0,1,2,...,N
//          mm    --- Physical dimension of CPM and CPD
//          ntype --- type of cut, either 2 or 3
// Output:  CPM(m,n) --- Pmn(z)
//          CPD(m,n) --- Pmn'(z)
// =========================================================

template <typename T, typename OutputMat1, typename OutputMat2>
void clpmn(std::complex<T> z, long ntype, OutputMat1 cpm, OutputMat2 cpd) {
    int m = cpm.extent(0) - 1;
    int n = cpm.extent(1) - 1;

    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            cpm(i, j) = 0;
            cpd(i, j) = 0;
        }
    }

    cpm(0, 0) = 1;
    if (n == 0) {
        return;
    }

    if ((std::abs(std::real(z)) == 1) && (std::imag(z) == 0)) {
        for (int i = 1; i <= n; i++) {
            cpm(0, i) = std::pow(std::real(z), i);
            cpd(0, i) = i * (i + 1) * std::pow(std::real(z), i + 1) / 2;
        }
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                if (i == 1) {
                    cpd(i, j) = std::numeric_limits<T>::infinity();
                } else if (i == 2) {
                    cpd(i, j) = -(j + 2) * (j + 1) * j * (j - 1) * std::pow(std::real(z), j + 1) / 4;
                }
            }
        }
        return;
    }

    std::complex<T> zq, zs;
    int ls;
    if (ntype == 2) {
        // sqrt(1 - z**2) with branch cut on |x|>1
        zs = (static_cast<T>(1) - z * z);
        zq = -std::sqrt(zs);
        ls = -1;
    } else {
        // sqrt(z**2 - 1) with branch cut between [-1, 1]
        zs = (z * z - static_cast<T>(1));
        zq = std::sqrt(zs);
        if (std::real(z) < 0) {
            zq = -zq;
        }
        ls = 1;
    }

    for (int i = 1; i <= m; i++) {
        // DLMF 14.7.15
        cpm(i, i) = static_cast<T>(2 * i - 1) * zq * cpm(i - 1, i - 1);
    }

    for (int i = 0; i <= (m > n - 1 ? n - 1 : m); i++) {
        // DLMF 14.10.7
        cpm(i, i + 1) = static_cast<T>(2 * i + 1) * z * cpm(i, i);
    }

    for (int i = 0; i <= m; i++) {
        for (int j = i + 2; j <= n; j++) {
            // DLMF 14.10.3
            cpm(i, j) = (static_cast<T>(2 * j - 1) * z * cpm(i, j - 1) - static_cast<T>(i + j - 1) * cpm(i, j - 2)) /
                        static_cast<T>(j - i);
        }
    }

    cpd(0, 0) = 0;
    for (int j = 1; j <= n; j++) {
        // DLMF 14.10.5
        cpd(0, j) = ls * static_cast<T>(j) * (z * cpm(0, j) - cpm(0, j - 1)) / zs;
    }

    for (int i = 1; i <= m; i++) {
        for (int j = i; j <= n; j++) {
            // derivative of DLMF 14.7.11 & DLMF 14.10.6 for type 3
            // derivative of DLMF 14.7.8 & DLMF 14.10.1 for type 2
            cpd(i, j) = static_cast<T>(ls) * (-static_cast<T>(i) * z * cpm(i, j) / zs +
                                              static_cast<T>((j + i) * (j - i + 1)) / zq * cpm(i - 1, j));
        }
    }
}

template <typename T, typename OutputMat1, typename OutputMat2>
void clpmn(std::complex<T> z, long ntype, bool m_signbit, OutputMat1 cpm, OutputMat2 cpd) {
    clpmn(z, ntype, cpm, cpd);

    int m = cpm.extent(0) - 1;
    int n = cpm.extent(1) - 1;

    if (m_signbit) {
        for (int j = 0; j < n + 1; ++j) {
            for (int i = 0; i < m + 1; ++i) {
                T fac = 0;
                if (i <= j) {
                    fac = std::tgamma(j - i + 1) / std::tgamma(j + i + 1);
                    if (ntype == 2) {
                        fac *= std::pow(-1, i);
                    }
                }

                cpm(i, j) *= fac;
                cpd(i, j) *= fac;
            }
        }
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
