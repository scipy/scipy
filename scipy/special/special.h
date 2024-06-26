#pragma once

#include "Python.h"

#include "special/bessel.h"
#include "special/legendre.h"
#include "special/specfun.h"
#include "special/sph_harm.h"

namespace {

template <typename T>
T legendre_p(long long int n, T z) {
    return special::legendre_p(n, z);
}

template <typename T>
void legendre_p(long long int n, T z, T &res, T &res_jac) {
    special::legendre_p(n, z, std::tie(res, res_jac));
}

template <typename T>
void legendre_p(long long int n, T z, T &res, T &res_jac, T &res_hess) {
    special::legendre_p(n, z, std::tie(res, res_jac, res_hess));
}

template <typename T, typename OutputVec1>
void legendre_p_all(T z, OutputVec1 res) {
    special::legendre_p_all(z, std::tie(res));
}

template <typename T, typename OutputVec1, typename OutputVec2>
void legendre_p_all(T z, OutputVec1 res, OutputVec2 res_jac) {
    special::legendre_p_all(z, std::tie(res, res_jac));
}

template <typename T, typename OutputVec1, typename OutputVec2, typename OutputVec3>
void legendre_p_all(T z, OutputVec1 res, OutputVec2 res_jac, OutputVec3 res_hess) {
    special::legendre_p_all(z, std::tie(res, res_jac, res_hess));
}

using special::assoc_legendre_norm;
using special::assoc_legendre_unnorm;

template <typename NormPolicy, typename T>
T assoc_legendre_p(NormPolicy norm, long long int n, long long int m, T z, long long int branch_cut) {
    return special::assoc_legendre_p(norm, n, m, z, branch_cut);
}

template <typename NormPolicy, typename T>
void assoc_legendre_p(NormPolicy norm, long long int n, long long int m, T z, long long int branch_cut, T &res,
                      T &res_jac) {
    special::assoc_legendre_p(norm, n, m, z, branch_cut, std::tie(res, res_jac));
}

template <typename NormPolicy, typename T>
void assoc_legendre_p(NormPolicy norm, long long int n, long long int m, T z, long long int branch_cut, T &res,
                      T &res_jac, T &res_hess) {
    special::assoc_legendre_p(norm, n, m, z, branch_cut, std::tie(res, res_jac, res_hess));
}

template <typename NormPolicy, typename T, typename OutputMat1>
void assoc_legendre_p_all(NormPolicy norm, T z, long long int branch_cut, OutputMat1 res) {
    special::assoc_legendre_p_all(norm, z, branch_cut, std::tie(res));
}

template <typename NormPolicy, typename T, typename OutputMat1, typename OutputMat2>
void assoc_legendre_p_all(NormPolicy norm, T z, long long int branch_cut, OutputMat1 res, OutputMat2 res_jac) {
    special::assoc_legendre_p_all(norm, z, branch_cut, std::tie(res, res_jac));
}

template <typename NormPolicy, typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void assoc_legendre_p_all(NormPolicy norm, T z, long long int branch_cut, OutputMat1 res, OutputMat2 res_jac,
                          OutputMat3 res_hess) {
    special::assoc_legendre_p_all(norm, z, branch_cut, std::tie(res, res_jac, res_hess));
}

template <typename T>
T sph_legendre_p(long long int n, long long int m, T theta) {
    return special::sph_legendre_p(n, m, theta);
}

template <typename T>
void sph_legendre_p(long long int n, long long int m, T theta, T &res, T &res_jac) {
    special::sph_legendre_p(n, m, theta, std::tie(res, res_jac));
}

template <typename T>
void sph_legendre_p(long long int n, long long int m, T theta, T &res, T &res_jac, T &res_hess) {
    special::sph_legendre_p(n, m, theta, std::tie(res, res_jac, res_hess));
}

template <typename T, typename OutputMat1>
void sph_legendre_p_all(T theta, OutputMat1 res) {
    special::sph_legendre_p_all(theta, std::tie(res));
}

template <typename T, typename OutputMat1, typename OutputMat2>
void sph_legendre_p_all(T theta, OutputMat1 res, OutputMat2 res_jac) {
    special::sph_legendre_p_all(theta, std::tie(res, res_jac));
}

template <typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void sph_legendre_p_all(T theta, OutputMat1 res, OutputMat2 res_jac, OutputMat3 res_hess) {
    special::sph_legendre_p_all(theta, std::tie(res, res_jac, res_hess));
}

template <typename T>
std::complex<T> sph_harm(long long int m, long long int n, T theta, T phi) {
    if (n < 0) {
        special::set_error("sph_harm", SF_ERROR_ARG, "n should not be negative");
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::abs(m) > n) {
        special::set_error("sph_harm", SF_ERROR_ARG, "m should not be greater than n");
        return std::numeric_limits<T>::quiet_NaN();
    }

    return special::sph_harm_y(n, m, phi, theta);
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

template <typename T>
std::complex<T> sph_harm_y(long long int n, long long int m, T theta, T phi) {
    return special::sph_harm_y(n, m, theta, phi);
}

template <typename T>
void sph_harm_y(long long int n, long long int m, T theta, T phi, std::complex<T> &res, std::complex<T> (&res_jac)[2]) {
    special::sph_harm_y(n, m, theta, phi, std::tie(res, res_jac));
}

template <typename T>
void sph_harm_y(long long int n, long long int m, T theta, T phi, std::complex<T> &res, std::complex<T> (&res_jac)[2],
                std::complex<T> (&res_hess)[2][2]) {
    special::sph_harm_y(n, m, theta, phi, std::tie(res, res_jac, res_hess));
}

template <typename T, typename OutputMat1>
void sph_harm_y_all(T theta, T phi, OutputMat1 res) {
    special::sph_harm_y_all(theta, phi, std::tie(res));
}

template <typename T, typename OutputMat1, typename OutputMat2>
void sph_harm_y_all(T theta, T phi, OutputMat1 res, OutputMat2 res_jac) {
    special::sph_harm_y_all(theta, phi, std::tie(res, res_jac));
}

template <typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void sph_harm_y_all(T theta, T phi, OutputMat1 res, OutputMat2 res_jac, OutputMat3 res_hess) {
    special::sph_harm_y_all(theta, phi, std::tie(res, res_jac, res_hess));
}

} // namespace
