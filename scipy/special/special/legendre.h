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

    int i, j, ls;
    double xq, xs;

    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++) {
            pm(i, j) = 0;
            pd(i, j) = 0;
        }
    }

    pm(0, 0) = 1;

    if (n == 0) {
        return;
    }

    if (fabs(x) == 1.0) {
        for (i = 1; i <= n; i++) {
            pm(0, i) = pow(x, i);
            pd(0, i) = 0.5 * i * (i + 1.0) * pow(x, i + 1);
        }

        for (i = 1; i <= m; i++) {
            for (j = 1; j <= n; j++) {
                if (i == 1) {
                    pd(i, j) = INFINITY;
                } else if (i == 2) {
                    pd(i, j) = -0.25 * (j + 2) * (j + 1) * j * (j - 1) * pow(x, j + 1);
                }
            }
        }
        return;
    }

    ls = (fabs(x) > 1.0 ? -1 : 1);
    xq = sqrt(ls * (1.0 - x * x));
    // Ensure connection to the complex-valued function for |x| > 1
    if (x < -1.0) {
        xq = -xq;
    }
    xs = ls * (1.0 - x * x);
    /* 30 */
    for (i = 1; i <= m; ++i) {
        pm(i, 0) = -ls * (2.0 * i - 1.0) * xq * pm(i - 1, 0);
    }
    /* 35 */
    for (i = 0; i <= (m > (n - 1) ? n - 1 : m); i++) {
        pm(i, 1) = (2.0 * i + 1.0) * x * pm(i, 0);
    }
    /* 40 */
    for (i = 0; i <= m; i++) {
        for (j = i + 2; j <= n; j++) {
            pm(i, j) = ((2.0 * j - 1.0) * x * pm(i, j - 1) - (i + j - 1.0) * pm(i, j - 2)) / (j - i);
        }
    }

    pd(0, 0) = 0.0;
    /* 45 */
    for (j = 1; j <= n; j++) {
        pd(0, j) = ls * j * (pm(0, j - 1) - x * pm(0, j)) / xs;
    }
    /* 50 */
    for (i = 1; i <= m; i++) {
        for (j = i; j <= n; j++) {
            pd(i, j) = ls * i * x * pm(i, j) / xs + (j + i) * (j - i + 1.0) / xq * pm(i - 1, j);
        }
    }
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

    int i, j, ls;
    std::complex<double> zq, zs;
    double x = z.real();
    double y = z.imag();

    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++) {
            cpm(i, j) = 0.0;
            cpd(i, j) = 0.0;
        }
    }

    cpm(0, 0) = 1.0;
    if (n == 0) {
        return;
    }

    if ((fabs(x) == 1.0) && (y == 0.0)) {
        for (i = 1; i <= n; i++) {
            cpm(0, i) = pow(x, i);
            cpd(0, i) = 0.5 * i * (i + 1) * pow(x, i + 1);
        }
        for (i = 1; i <= m; i++) {
            for (j = 1; j <= n; j++) {
                if (i == 1) {
                    cpd(i, j) = INFINITY;
                } else if (i == 2) {
                    cpd(i, j) = -0.25 * (j + 2) * (j + 1) * j * (j - 1) * pow(x, j + 1);
                }
            }
        }
        return;
    }

    if (ntype == 2) {
        // sqrt(1 - z**2) with branch cut on |x|>1
        zs = (1.0 - z * z);
        zq = -std::sqrt(zs);
        ls = -1;
    } else {
        // sqrt(z**2 - 1) with branch cut between [-1, 1]
        zs = (z * z - 1.0);
        zq = std::sqrt(zs);
        if (x < 0.) {
            zq = -zq;
        }
        ls = 1;
    }

    for (i = 1; i <= m; i++) {
        // DLMF 14.7.15
        cpm(i, i) = (2. * i - 1.) * zq * cpm(i - 1, i - 1);
    }

    for (i = 0; i <= (m > n - 1 ? n - 1 : m); i++) {
        // DLMF 14.10.7
        cpm(i, i + 1) = (2. * i + 1) * z * cpm(i, i);
    }

    for (i = 0; i <= m; i++) {
        for (j = i + 2; j <= n; j++) {
            // DLMF 14.10.3
            cpm(i, j) = ((2. * j - 1) * z * cpm(i, j - 1) - static_cast<double>(i + j - 1) * cpm(i, j - 2)) /
                        static_cast<double>(j - i);
        }
    }

    cpd(0, 0) = 0.0;
    for (j = 1; j <= n; j++) {
        // DLMF 14.10.5
        cpd(0, j) = ls * static_cast<double>(j) * (z * cpm(0, j) - cpm(0, j - 1)) / zs;
    }

    for (i = 1; i <= m; i++) {
        for (j = i; j <= n; j++) {
            // derivative of DLMF 14.7.11 & DLMF 14.10.6 for type 3
            // derivative of DLMF 14.7.8 & DLMF 14.10.1 for type 2
            cpd(i, j) = static_cast<double>(ls) *
                        (-static_cast<double>(i) * z * cpm(i, j) / zs + (j + i) * (j - i + 1.0) / zq * cpm(i - 1, j));
        }
    }
}

} // namespace special
