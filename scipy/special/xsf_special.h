#pragma once

#include "Python.h"

#include "xsf/bessel.h"
#include "xsf/dual.h"
#include "xsf/legendre.h"
#include "xsf/specfun.h"
#include "xsf/sph_harm.h"

// This header exists to add behaviors to special functions from the xsf library,
// either because they involve some Python-specific features or because there are
// legacy behaviors that need to be preserved

namespace {

template <typename T, size_t N>
xsf::dual<T, N> legendre_p(long long int n, T z) {
    return xsf::legendre_p(n, xsf::dual_var<N>(z));
}

template <typename T, typename OutputVec>
void legendre_p_all(T z, OutputVec res) {
    static constexpr size_t N = OutputVec::value_type::max_order();

    xsf::legendre_p_all(xsf::dual_var<N>(z), res);
}

using xsf::assoc_legendre_norm;
using xsf::assoc_legendre_unnorm;

template <size_t N, typename NormPolicy, typename T>
xsf::dual<T, N> assoc_legendre_p(NormPolicy norm, long long int n, long long int m, T z, long long int branch_cut) {
    return xsf::assoc_legendre_p(norm, n, m, xsf::dual_var<N>(z), branch_cut);
}

template <typename NormPolicy, typename T, typename OutputMat>
void assoc_legendre_p_all(NormPolicy norm, T z, long long int branch_cut, OutputMat res) {
    static constexpr size_t N = OutputMat::value_type::max_order();

    xsf::assoc_legendre_p_all(norm, xsf::dual_var<N>(z), branch_cut, res);
}

template <typename T, size_t N>
xsf::dual<T, N> sph_legendre_p(long long int n, long long int m, T theta) {
    return xsf::sph_legendre_p(n, m, xsf::dual_var<N>(theta));
}

template <typename T, typename OutputMat>
void sph_legendre_p_all(T theta, OutputMat res) {
    static constexpr size_t N = OutputMat::value_type::max_order();

    xsf::sph_legendre_p_all(xsf::dual_var<N>(theta), res);
}

template <typename T, size_t N>
xsf::dual<std::complex<T>, N, N> sph_harm_y(long long int n, long long int m, T theta, T phi) {
    return xsf::sph_harm_y(n, m, xsf::dual_var<N, N>(theta, 0), xsf::dual_var<N, N>(phi, 1));
}

template <typename T>
std::complex<T> sph_harm(long long int m, long long int n, T theta, T phi) {
    if (n < 0) {
        xsf::set_error("sph_harm", SF_ERROR_ARG, "n should not be negative");
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::abs(m) > n) {
        xsf::set_error("sph_harm", SF_ERROR_ARG, "m should not be greater than n");
        return std::numeric_limits<T>::quiet_NaN();
    }

    return xsf::sph_harm_y(n, m, phi, theta);
}

template <typename T>
std::complex<T> sph_harm(T m, T n, T theta, T phi) {
    if (static_cast<long>(m) != m || static_cast<long>(n) != n) {
        PyGILState_STATE gstate = PyGILState_Ensure();
        PyErr_WarnEx(PyExc_RuntimeWarning, "floating point number truncated to an integer", 1);
        PyGILState_Release(gstate);
    }

    return sph_harm(static_cast<long>(m), static_cast<long>(n), theta, phi);
}

template <typename T, typename OutputMat>
void sph_harm_y_all(T theta, T phi, OutputMat res) {
    static constexpr size_t N = OutputMat::value_type::max_order();

    xsf::sph_harm_y_all(xsf::dual_var<N, N>(theta, 0), xsf::dual_var<N, N>(phi, 1), res);
}

} // namespace
