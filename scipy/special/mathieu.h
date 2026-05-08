#pragma once

#include <algorithm>
#include <vector>

#include <xsf/error.h>
#include <xsf/mathieu.h>
#include <xsf/mathieu_very_new.h>

#include "tridiagonal.h"

namespace special {

template <xsf::mathieu::Parity FuncParity, typename T>
struct mathieu_cv {

    std::vector<double> D;
    std::vector<double> E;
    eigvalsh_tridiagonal solver;

    T operator()(T m, T q) {
        using namespace xsf::mathieu;
        auto constexpr Even = Parity::Even;
        auto constexpr Odd = Parity::Odd;

        if constexpr (FuncParity == Even) {
            if ((m < 0) || (m != std::floor(m))) {
                xsf::set_error("mathieu_a", SF_ERROR_DOMAIN, NULL);
                return std::numeric_limits<double>::quiet_NaN();
            }
        } else {
            if ((m <= 0) || (m != std::floor(m))) {
                xsf::set_error("mathieu_b", SF_ERROR_DOMAIN, NULL);
                return std::numeric_limits<double>::quiet_NaN();
            }
        }

        if (m > 500) {
            // Fallback to Zhang and Jin implementation for very high order.
            if constexpr (FuncParity == Even) {
                return xsf::cem_cva(m, q);
            }
            return xsf::sem_cva(m, q);
        }

        auto int_m = static_cast<CBLAS_INT>(m);
        auto N = int_m + 25;

        D.resize(N);
        E.resize(N - 1);

        if (int_m % 2) {
            make_matrix<FuncParity, Odd>(q, xsf::numpy::as_mdspan(D), xsf::numpy::as_mdspan(E));
        } else {
            make_matrix<FuncParity, Even>(q, xsf::numpy::as_mdspan(D), xsf::numpy::as_mdspan(E));
        }

        auto status = solver(D, E);
        if (status != SF_ERROR_OK) {
            if constexpr (FuncParity == Even) {
                xsf::set_error("mathieu_a", status, NULL);
            } else {
                xsf::set_error("mathieu_b", status, NULL);
            }
            return std::numeric_limits<double>::quiet_NaN();
        }

        auto idx = cv_index<FuncParity>(int_m);
        return static_cast<T>(D[idx]);
    }
};

template <xsf::mathieu::Parity FuncParity>
struct mathieu_coeffs {
    std::vector<double> D;
    std::vector<double> E;
    std::vector<double> Z;
    eigh_tridiagonal solver;

    sf_error_t operator()(int m, double q, std::vector<double> &X) {
        using namespace xsf::mathieu;
        auto constexpr Even = Parity::Even;
        auto constexpr Odd = Parity::Odd;

        auto N = X.size();

        if (N == 0) {
            return SF_ERROR_OK;
        }

        D.resize(N);
        E.resize(N - 1);
        Z.resize(N * N);

        if (m % 2) {
            make_matrix<FuncParity, Odd>(q, xsf::numpy::as_mdspan(D), xsf::numpy::as_mdspan(E));
        } else {
            make_matrix<FuncParity, Even>(q, xsf::numpy::as_mdspan(D), xsf::numpy::as_mdspan(E));
        }

        auto status = solver(D, E, Z);
        if (status != SF_ERROR_OK) {
            return status;
        }

        auto col = cv_index<FuncParity>(m);
        std::copy_n(Z.begin() + col * N, N, X.begin());
        if constexpr (FuncParity == Even) {
            if (m % 2 == 0) {
                X[0] /= M_SQRT2;
            }
        }

        return SF_ERROR_OK;
    }
};

template <xsf::mathieu::Parity FuncParity, typename T>
struct mathieu_xem {
    double last_q{std::numeric_limits<double>::quiet_NaN()};
    int last_m{-1};
    mathieu_coeffs<FuncParity> get_coefs;
    std::vector<double> coefs;
    void operator()(T m, T q, T x, T &out, T &out_diff) {
        using namespace xsf::mathieu;
        auto constexpr Even = Parity::Even;
        auto constexpr Odd = Parity::Odd;

        double q_d = static_cast<double>(q);
        /* SciPy's mathieu functions take angles in degrees.
         * Convert to radians for internal calculation.*/
        double v = static_cast<double>(x) * (M_PI / 180.0);
        double out_d, out_diff_d;

        if constexpr (FuncParity == Even) {
            if ((m < 0) || m != std::floor(m)) {
                out = std::numeric_limits<T>::quiet_NaN();
                out_diff = std::numeric_limits<T>::quiet_NaN();
                xsf::set_error("mathieu_cem", SF_ERROR_DOMAIN, NULL);
                last_m = -1; // invalidate cache upon error
                return;
            }
        } else {
            if ((m <= 0) || m != std::floor(m)) {
                out = std::numeric_limits<T>::quiet_NaN();
                out_diff = std::numeric_limits<T>::quiet_NaN();
                xsf::set_error("mathieu_sem", SF_ERROR_DOMAIN, NULL);
                last_m = -1; // invalidate cache upon error
                return;
            }
        }

        auto int_m = static_cast<int>(m);

        if (q_d != last_q || int_m != last_m) {
            auto N = get_partial_sum_N(int_m, q_d);
            coefs.resize(N);
            auto status = get_coefs(int_m, q_d, coefs);
            if (status != SF_ERROR_OK) {
                out = std::numeric_limits<T>::quiet_NaN();
                out_diff = std::numeric_limits<T>::quiet_NaN();
                last_m = -1; // invalidate cache upon error
                if constexpr (FuncParity == Even) {
                    xsf::set_error("mathieu_cem", status, NULL);
                } else {
                    xsf::set_error("mathieu_sem", status, NULL);
                }
                return;
            }
            last_q = q;
            last_m = int_m;
        }
        if (int_m % 2) {
            sum_fourier_series<FuncParity, Odd>(xsf::numpy::as_mdspan(coefs), v, out_d, out_diff_d);
        } else {
            sum_fourier_series<FuncParity, Even>(xsf::numpy::as_mdspan(coefs), v, out_d, out_diff_d);
        }
        out = static_cast<T>(out_d);
        /* Need to map to the derivative with respect to an angle in radians. */
        out_diff = static_cast<T>(out_diff_d * (M_PI / 180.0));
    }
};

} // namespace special
