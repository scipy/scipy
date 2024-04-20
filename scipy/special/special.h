#pragma once

#include "Python.h"

#include "special/legendre.h"
#include "special/specfun.h"
#include "special/sph_harm.h"

namespace {

template <typename T, typename OutputVec1, typename OutputVec2>
void lpn(T z, OutputVec1 p, OutputVec2 p_jac) {
    unsigned int n = p.extent(0) - 1;

    special::legendre_p_jac(n, z, [p, p_jac](unsigned int j, T z, T value, T value_jac) {
        p(j) = value;
        p_jac(j) = value_jac;
    });
}

template <typename T, typename OutputMat1, typename OutputMat2>
void lpmn(T x, bool m_signbit, OutputMat1 p, OutputMat2 p_jac) {
    unsigned int m = p.extent(0) - 1;
    unsigned int n = p.extent(1) - 1;

    for (unsigned int i = 0; i <= m; ++i) {
        special::assoc_legendre_p_jac(
            n, i, m_signbit, x,
            [p, p_jac](unsigned int j, unsigned int i, bool i_signbit, T x, T value, T value_jac) {
                p(i, j) = value;
                p_jac(i, j) = value_jac;
            }
        );
    }
}

template <typename T>
std::complex<T> sph_harm(long m, long n, T theta, T phi) {
    if (std::abs(m) > n) {
        special::set_error("sph_harm", SF_ERROR_ARG, "m should not be greater than n");
        return std::numeric_limits<T>::quiet_NaN();
    }

    return special::sph_harm(m, n, theta, phi);
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

} // namespace
