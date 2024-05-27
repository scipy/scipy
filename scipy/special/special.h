#pragma once

#include "Python.h"

#include "special/legendre.h"
#include "special/specfun.h"
#include "special/sph_harm.h"

namespace {

template <typename T, typename OutputVec1, typename OutputVec2>
void lpn(T z, OutputVec1 p, OutputVec2 p_jac) {
    special::legendre_all(z, p);
    special::legendre_all_jac(z, p, p_jac);
}

template <typename T, typename OutputMat1, typename OutputMat2>
void lpmn(T x, bool m_signbit, OutputMat1 p, OutputMat2 p_jac) {
    special::assoc_legendre_all(x, m_signbit, p);
    special::assoc_legendre_all_jac(x, m_signbit, p, p_jac);
}

template <typename T>
std::complex<T> sph_harm(long m, long n, T theta, T phi) {
    if (std::abs(m) > n) {
        special::set_error("sph_harm", SF_ERROR_ARG, "m should not be greater than n");
        return NAN;
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
