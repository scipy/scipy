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
    special::grad<T, 1> tmp;
    special::legendre_p(n, z, tmp);

    std::tie(res, res_jac) = tmp.refs_as_tuple();
}

template <typename T>
void legendre_p(long long int n, T z, T &res, T &res_jac, T &res_hess) {
    special::grad<T, 2> tmp;
    special::legendre_p(n, z, tmp);

    std::tie(res, res_jac, res_hess) = tmp.refs_as_tuple();
}

template <typename T, typename OutputVec1>
void legendre_p_all(T z, OutputVec1 res) {
    special::grad<OutputVec1, 0> tmp{res};
    special::legendre_p_all(z, tmp);
}

template <typename T, typename OutputVec1, typename OutputVec2>
void legendre_p_all(T z, OutputVec1 res, OutputVec2 res_jac) {
    special::grad<OutputVec1, 1> tmp{res, res_jac};
    special::legendre_p_all(z, tmp);
}

template <typename T, typename OutputVec1, typename OutputVec2, typename OutputVec3>
void legendre_p_all(T z, OutputVec1 res, OutputVec2 res_jac, OutputVec3 res_hess) {
    special::grad<OutputVec1, 2> tmp{res, res_jac, res_hess};
    special::legendre_p_all(z, tmp);
}

using special::assoc_legendre_norm;
using special::assoc_legendre_unnorm;

template <typename NormPolicy, typename T>
T assoc_legendre_p(NormPolicy norm, long long int n, long long int m, T z) {
    return special::assoc_legendre_p(norm, n, m, z);
}

template <typename NormPolicy, typename T>
void assoc_legendre_p(NormPolicy norm, long long int n, long long int m, T z, T &res, T &res_jac) {
    special::grad<T, 1> tmp;
    special::assoc_legendre_p(norm, n, m, z, tmp);

    std::tie(res, res_jac) = tmp.refs_as_tuple();
}

template <typename NormPolicy, typename T>
void assoc_legendre_p(NormPolicy norm, long long int n, long long int m, T z, T &res, T &res_jac, T &res_hess) {
    special::grad<T, 2> tmp;
    special::assoc_legendre_p(norm, n, m, z, tmp);

    std::tie(res, res_jac, res_hess) = tmp.refs_as_tuple();
}

template <typename NormPolicy, typename T, typename OutputMat1>
void assoc_legendre_p_all(NormPolicy norm, T z, OutputMat1 res) {
    special::grad<OutputMat1, 0> tmp(res);
    special::assoc_legendre_p_all(norm, z, tmp);
}

template <typename NormPolicy, typename T, typename OutputMat1, typename OutputMat2>
void assoc_legendre_p_all(NormPolicy norm, T z, OutputMat1 res, OutputMat2 res_jac) {
    special::grad<OutputMat1, 1> tmp(res, res_jac);
    special::assoc_legendre_p_all(norm, z, tmp);
}

template <typename NormPolicy, typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void assoc_legendre_p_all(NormPolicy norm, T z, OutputMat1 res, OutputMat2 res_jac, OutputMat3 res_hess) {
    special::grad<OutputMat1, 2> tmp(res, res_jac, res_hess);
    special::assoc_legendre_p_all(norm, z, tmp);
}

template <typename NormPolicy, typename T>
std::complex<T>
multi_assoc_legendre_p(NormPolicy norm, long long int n, long long int m, long long int type, std::complex<T> z) {
    return special::multi_assoc_legendre_p(norm, n, m, type, z);
}

template <typename NormPolicy, typename T>
void multi_assoc_legendre_p(
    NormPolicy norm, long long int n, long long int m, long long int type, std::complex<T> z, std::complex<T> &res,
    std::complex<T> &res_jac
) {
    special::grad<std::complex<T>, 1> tmp;
    special::multi_assoc_legendre_p(norm, n, m, type, z, tmp);

    std::tie(res, res_jac) = tmp.refs_as_tuple();
}

template <typename NormPolicy, typename T>
void multi_assoc_legendre_p(
    NormPolicy norm, long long int n, long long int m, long long int type, std::complex<T> z, std::complex<T> &res,
    std::complex<T> &res_jac, std::complex<T> &res_hess
) {
    special::grad<std::complex<T>, 2> tmp;
    special::multi_assoc_legendre_p(norm, n, m, type, z, tmp);

    std::tie(res, res_jac, res_hess) = tmp.refs_as_tuple();
}

template <typename NormPolicy, typename T, typename OutputMat1>
void multi_assoc_legendre_p_all(NormPolicy norm, int type, T z, OutputMat1 res) {
    special::grad<OutputMat1, 0> tmp(res);
    special::multi_assoc_legendre_p_all(norm, type, z, tmp);
}

template <typename NormPolicy, typename T, typename OutputMat1, typename OutputMat2>
void multi_assoc_legendre_p_all(NormPolicy norm, int type, T z, OutputMat1 res, OutputMat2 res_jac) {
    special::grad<OutputMat1, 1> tmp(res, res_jac);
    special::multi_assoc_legendre_p_all(norm, type, z, tmp);
}

template <typename NormPolicy, typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void multi_assoc_legendre_p_all(
    NormPolicy norm, int type, T z, OutputMat1 res, OutputMat2 res_jac, OutputMat3 res_hess
) {
    special::grad<OutputMat1, 2> tmp(res, res_jac, res_hess);
    special::multi_assoc_legendre_p_all(norm, type, z, tmp);
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
