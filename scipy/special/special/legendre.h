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

    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            cpm(i, j) = 0.0;
            cpd(i, j) = 0.0;
        }
    }

    cpm(0, 0) = 1.0;
    if (n == 0) {
        return;
    }

    if ((std::abs(std::real(z)) == 1.0) && (std::imag(z) == 0.0)) {
        for (int i = 1; i <= n; i++) {
            cpm(0, i) = pow(std::real(z), i);
            cpd(0, i) = 0.5 * i * (i + 1) * pow(std::real(z), i + 1);
        }
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                if (i == 1) {
                    cpd(i, j) = INFINITY;
                } else if (i == 2) {
                    cpd(i, j) = -0.25 * (j + 2) * (j + 1) * j * (j - 1) * pow(std::real(z), j + 1);
                }
            }
        }
        return;
    }

    std::complex<double> zq, zs;
    int ls;
    if (ntype == 2) {
        // sqrt(1 - z**2) with branch cut on |x|>1
        zs = (1.0 - z * z);
        zq = -std::sqrt(zs);
        ls = -1;
    } else {
        // sqrt(z**2 - 1) with branch cut between [-1, 1]
        zs = (z * z - 1.0);
        zq = std::sqrt(zs);
        if (std::real(z) < 0.) {
            zq = -zq;
        }
        ls = 1;
    }

    /*
        for (int i = 1; i <= m; i++) {
            // DLMF 14.7.15
            cpm[i*(n + 2)] = (2.*i - 1.)*zq*cpm[(i-1)*(n+2)];
        }
        for (int i = 0; i <= (m > n-1 ? n-1 : m); i++) {
            // DLMF 14.10.7
            cpm[i*(n + 2) + 1] = (2.*i + 1)*z*cpm[i*(n + 2)];
        }
    */

    for (int i = 1; i <= m; i++) {
        // DLMF 14.7.15
        cpm(i, i) = (2. * i - 1.) * zq * cpm(i - 1, i - 1);
    }

    for (int i = 0; i <= (m > n - 1 ? n - 1 : m); i++) {
        // DLMF 14.10.7
        cpm(i, i + 1) = (2. * i + 1) * z * cpm(i, i);
    }

    for (int i = 0; i <= m; i++) {
        for (int j = i + 2; j <= n; j++) {
            // DLMF 14.10.3
            cpm(i, j) = ((2. * j - 1) * z * cpm(i, j - 1) - static_cast<double>(i + j - 1) * cpm(i, j - 2)) /
                        static_cast<double>(j - i);
        }
    }

    cpd(0, 0) = 0.0;
    for (int j = 1; j <= n; j++) {
        // DLMF 14.10.5
        cpd(0, j) = ls * static_cast<double>(j) * (z * cpm(0, j) - cpm(0, j - 1)) / zs;
    }

    for (int i = 1; i <= m; i++) {
        for (int j = i; j <= n; j++) {
            // derivative of DLMF 14.7.11 & DLMF 14.10.6 for type 3
            // derivative of DLMF 14.7.8 & DLMF 14.10.1 for type 2
            cpd(i, j) = static_cast<double>(ls) *
                        (-static_cast<double>(i) * z * cpm(i, j) / zs + (j + i) * (j - i + 1.0) / zq * cpm(i - 1, j));
        }
    }
}

} // namespace special
