#pragma once

#include "config.h"
#include "mdspan.h"
#include "specfun/specfun.h"

namespace special {

template <typename T>
void lpn(T x, std::mdspan<T, std::dextents<int, 1>, std::layout_stride> pn,
         std::mdspan<T, std::dextents<int, 1>, std::layout_stride> pd);

// Translated into C++ by SciPy developers in 2024.
// Original comments appear below.
//
// ===============================================
// Purpose: Compute Legendre polynomials Pn(x)
//          and their derivatives Pn'(x)
// Input :  x --- Argument of Pn(x)
//          n --- Degree of Pn(x) ( n = 0,1,...)
// Output:  PN(n) --- Pn(x)
//          PD(n) --- Pn'(x)
// ===============================================

template <>
void lpn(double x, std::mdspan<double, std::dextents<int, 1>, std::layout_stride> pn,
         std::mdspan<double, std::dextents<int, 1>, std::layout_stride> pd) {
    int n = pn.size() - 1;

    int k;
    double p0, p1, pf;
    pn[0] = 1.0;
    pn[1] = x;
    pd[0] = 0.0;
    pd[1] = 1.0;
    p0 = 1.0;
    p1 = x;
    for (k = 2; k <= n; k++) {
        pf = (2.0 * k - 1.0) / k * x * p1 - (k - 1.0) / k * p0;
        pn[k] = pf;
        if (fabs(x) == 1.0) {
            pd[k] = 0.5 * pow(x, k + 1) * k * (k + 1);
        } else {
            pd[k] = k * (p1 - x * pf) / (1.0 - x * x);
        }
        p0 = p1;
        p1 = pf;
    }
}

// Translated into C++ by SciPy developers in 2024.
// Original comments appear below.
//
// ==================================================
// Purpose: Compute Legendre polynomials Pn(z) and
//          their derivatives Pn'(z) for a complex
//          argument
// Input :  x --- Real part of z
//          y --- Imaginary part of z
//          n --- Degree of Pn(z), n = 0,1,2,...
// Output:  CPN(n) --- Pn(z)
//          CPD(n) --- Pn'(z)
// ==================================================

template <>
void lpn(std::complex<double> z, std::mdspan<std::complex<double>, std::dextents<int, 1>, std::layout_stride> cpn,
         std::mdspan<std::complex<double>, std::dextents<int, 1>, std::layout_stride> cpd) {
    int n = cpn.size() - 1;

    int k;
    std::complex<double> cp0, cp1, cpf;

    cpn[0] = 1.0;
    cpn[1] = z;
    cpd[0] = 0.0;
    cpd[1] = 1.0;
    cp0 = 1.0;
    cp1 = z;
    for (k = 2; k <= n; k++) {
        cpf = (2.0 * k - 1.0) / k * z * cp1 - (k - 1.0) / k * cp0;
        cpn[k] = cpf;
        if (z == 1.0) {
            cpd[k] = 0.5 * pow(z.real(), k + 1) * k * (k + 1.0);
        } else {
            cpd[k] = static_cast<double>(k) * (cp1 - z * cpf) / (1.0 - z * z);
        }
        cp0 = cp1;
        cp1 = cpf;
    }
    return;
}

void lpmn(double x, std::mdspan<double, std::dextents<int, 2>, std::layout_stride> pm,
          std::mdspan<double, std::dextents<int, 2>, std::layout_stride> pd) {
    int m = pm.extent(0) - 1;
    int n = pm.extent(1) - 1;

    return specfun::lpmn(m, n, x, pm.data_handle(), pd.data_handle());
}

} // namespace special
