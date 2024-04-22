#pragma once

#include "Python.h"

#include "special/legendre.h"
#include "special/specfun.h"
#include "special/sph_harm.h"

namespace {

template <typename T>
T lpn(long n, T z) {
    return special::legendre_p(n, z);
}

template <typename T>
void lpn(long n, T z, T &p, T &p_jac) {
    special::legendre_p(n, z, p, p_jac);
}

template <typename T>
void lpn(long n, T z, T &p, T &p_jac, T &p_hess) {
    special::legendre_p(n, z, p, p_jac, p_hess);
}

template <typename T, typename OutputVec>
void lpn_all(T z, OutputVec p) {
    special::legendre_p_all(z, p);
}

template <typename T, typename OutputVec1, typename OutputVec2>
void lpn_all(T z, OutputVec1 p, OutputVec2 p_jac) {
    special::legendre_p_all(z, p, p_jac);
}

template <typename T, typename OutputVec1, typename OutputVec2, typename OutputVec3>
void lpn_all(T z, OutputVec1 p, OutputVec2 p_jac, OutputVec3 p_hess) {
    special::legendre_p_all(z, p, p_jac, p_hess);
}

template <typename T, typename OutputMat1>
void lpmn_all(T z, OutputMat1 p) {
    special::assoc_legendre_p_all(z, p);
}

template <typename T, typename OutputMat1, typename OutputMat2>
void lpmn_all(T z, OutputMat1 p, OutputMat2 p_jac) {
    special::assoc_legendre_p_all(z, p, p_jac);
}

template <typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void lpmn_all(T z, OutputMat1 p, OutputMat2 p_jac, OutputMat3 p_hess) {
    special::assoc_legendre_p_all(z, p, p_jac, p_hess);
}

template <typename T>
std::complex<T> sph_harm(long m, long n, T theta, T phi) {
    if (n < 0) {
        special::set_error("sph_harm", SF_ERROR_ARG, "n should not be negative");
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::abs(m) > n) {
        special::set_error("sph_harm", SF_ERROR_ARG, "m should not be greater than n");
        return std::numeric_limits<T>::quiet_NaN();
    }

    return special::sph_harm(n, m, theta, phi);
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
