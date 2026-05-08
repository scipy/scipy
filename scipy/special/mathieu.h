#pragma once

#include <algorithm>
#include <vector>

#include <xsf/error.h>
#include <xsf/mathieu.h>
#include <xsf/mathieu_very_new.h>

#include "tridiagonal.h"

namespace special {

template <xsf::mathieu::Parity FuncParity>
double mathieu_cv(double m, double q) {
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

    std::vector<double> D(N);
    std::vector<double> E(N - 1);

    if (int_m % 2) {
        make_matrix<FuncParity, Odd>(q, xsf::numpy::as_mdspan(D), xsf::numpy::as_mdspan(E));
    } else {
        make_matrix<FuncParity, Even>(q, xsf::numpy::as_mdspan(D), xsf::numpy::as_mdspan(E));
    }

    auto status = eigvalsh_tridiagonal(D, E);
    if (status != SF_ERROR_OK) {
        if constexpr (FuncParity == Even) {
            xsf::set_error("mathieu_a", status, NULL);
        } else {
            xsf::set_error("mathieu_b", status, NULL);
        }
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto idx = cv_index<FuncParity>(int_m);
    return D[idx];
}

template <xsf::mathieu::Parity FuncParity>
sf_error_t mathieu_coeffs(double q, int m, std::vector<double> &AA) {
    using namespace xsf::mathieu;
    auto constexpr Even = Parity::Even;
    auto constexpr Odd = Parity::Odd;

    auto N = AA.size();

    if (N == 0) {
        return SF_ERROR_OK;
    }

    std::vector<double> D(N);
    std::vector<double> E(N - 1);
    std::vector<double> Z(N * N);

    if (m % 2) {
        make_matrix<FuncParity, Odd>(q, xsf::numpy::as_mdspan(D), xsf::numpy::as_mdspan(E));
    } else {
        make_matrix<FuncParity, Even>(q, xsf::numpy::as_mdspan(D), xsf::numpy::as_mdspan(E));
    }

    auto status = eigh_tridiagonal(D, E, Z);
    if (status != SF_ERROR_OK) {
        return status;
    }

    auto idx = cv_index<FuncParity>(m);
    std::copy_n(Z.begin() + idx * N, N, AA.begin());

    return SF_ERROR_OK;
}

double mathieu_a(double m, double q) { return mathieu_cv<xsf::mathieu::Parity::Even>(m, q); }
float mathieu_a(float m, float q) { return mathieu_a(static_cast<double>(m), static_cast<double>(q)); }

double mathieu_b(double m, double q) { return mathieu_cv<xsf::mathieu::Parity::Odd>(m, q); }
float mathieu_b(float m, float q) { return mathieu_b(static_cast<double>(m), static_cast<double>(q)); }

} // namespace special
