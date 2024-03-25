#pragma once

#include "config.h"
#include "mdspan.h"
#include "specfun/specfun.h"

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

template <typename T>
void lpn(T z, std::mdspan<T, std::dextents<int, 1>, std::layout_stride> pn,
         std::mdspan<T, std::dextents<int, 1>, std::layout_stride> pd) {
    int n = pn.extent(0) - 1;

    pn[0] = 1;
    pn[1] = z;
    pd[0] = 0;
    pd[1] = 1;

    T p0 = 1;
    T p1 = z;
    T pf;
    for (int k = 2; k <= n; k++) {
        pf = static_cast<T>(2 * k - 1) * z * p1 / static_cast<T>(k) - static_cast<T>(k - 1) * p0 / static_cast<T>(k);
        pn(k) = pf;
        if (std::abs(std::real(z)) == 1.0 && std::imag(z) == 0) {
            pd(k) = std::pow(std::real(z), k + 1) * k * (k + 1) / 2;
        } else {
            pd(k) = static_cast<T>(k) * (p1 - z * pf) / (static_cast<T>(1) - z * z);
        }
        p0 = p1;
        p1 = pf;
    }
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

void lpmn(double x, std::mdspan<double, std::dextents<int, 2>, std::layout_stride> pm,
          std::mdspan<double, std::dextents<int, 2>, std::layout_stride> pd) {
    int m = pm.extent(0) - 1;
    int n = pm.extent(1) - 1;

    specfun::lpmn(m, n, x, pm.data_handle(), pd.data_handle());
}

// Translated into C++ by SciPy developers in 2024.
// Original comments appear below.
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

void clpmn(std::complex<double> z, long ntype,
           std::mdspan<std::complex<double>, std::dextents<int, 2>, std::layout_stride> cpm,
           std::mdspan<std::complex<double>, std::dextents<int, 2>, std::layout_stride> cpd) {
    int m = cpm.extent(0) - 1;
    int n = cpm.extent(1) - 1;

    specfun::clpmn(z, m, n, ntype, cpm.data_handle(), cpd.data_handle());
}

} // namespace special
