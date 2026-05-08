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

    sf_error_t operator()(double q, int m, std::vector<double> &AA) {
        using namespace xsf::mathieu;
        auto constexpr Even = Parity::Even;
        auto constexpr Odd = Parity::Odd;

        auto N = AA.size();

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

        auto idx = cv_index<FuncParity>(m);
        std::copy_n(Z.begin() + idx * N, N, AA.begin());

        return SF_ERROR_OK;
    }
};

// scalar versions for cython special.
template <typename T>
T mathieu_a(T m, T q) {
    return mathieu_cv<xsf::mathieu::Parity::Even, T>{}(m, q);
}

} // namespace special
