#pragma once

#include <vector>

#include <xsf/error.h>
#include <xsf/mathieu_very_new.h>

#include "tridiagonal.h"

namespace special {

double mathieu_a(double m, double q) {
    if ((m < 0) || (m != std::floor(m))) {
      xsf::set_error("mathieu_a", SF_ERROR_DOMAIN, NULL);
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (m > 500) {
      // Don't support absurdly larger orders for now.
      xsf::set_error("mathieu_a", SF_ERROR_NO_RESULT, NULL);
      return std::numeric_limits<double>::quiet_NaN();
    }

    auto int_m = static_cast<CBLAS_INT>(m);
    auto N = int_m + 25;

    std::vector<double> D(N);
    std::vector<double> E(N-1);

    if (int_m % 2) {
      xsf::mathieu::make_matrix_eo(q, xsf::numpy::as_mdspan(D), xsf::numpy::as_mdspan(E));
    } else {
      xsf::mathieu::make_matrix_ee(q, xsf::numpy::as_mdspan(D), xsf::numpy::as_mdspan(E));
    }

    auto status = eigvalsh_tridiagonal(D, E);
    if (status != SF_ERROR_OK) {
      xsf::set_error("mathieu_a", status, NULL);
    }

    auto idx = int_m % 2 ? (int_m - 1) / 2 : int_m / 2;
    return D[idx];
}

float mathieu_a(float m, float q) {
    return mathieu_a(static_cast<double>(m), static_cast<double>(q)); 
}

}
