#pragma once

#include "Python.h"

#include "special/legendre.h"
#include "special/specfun.h"
#include "special/sph_harm.h"

namespace {

template <typename T>
T lpn(long long n, T z) {
    return special::legendre_p(n, z);
}

template <typename T>
void lpn(long long n, T z, T &res, T &res_jac) {
    special::legendre_p(n, z, res, res_jac);
}

template <typename T>
void lpn(long long n, T z, T &res, T &res_jac, T &res_hess) {
    special::legendre_p(n, z, res, res_jac, res_hess);
}

template <typename T, typename OutputVec>
void lpn_all(T z, OutputVec res) {
    special::legendre_p_all(z, res);
}

template <typename T, typename OutputVec1, typename OutputVec2>
void lpn_all(T z, OutputVec1 res, OutputVec2 res_jac) {
    special::legendre_p_all(z, res, res_jac);
}

template <typename T, typename OutputVec1, typename OutputVec2, typename OutputVec3>
void lpn_all(T z, OutputVec1 res, OutputVec2 res_jac, OutputVec3 res_hess) {
    special::legendre_p_all(z, res, res_jac, res_hess);
}

template <typename T>
T lpmn(long m, long n, T z) {
    int type;
    if (std::abs(z) <= 1) {
        type = 2;
    } else {
        type = 3;
    }

    return special::assoc_legendre_p(n, m, type, z);
}

template <typename T>
void lpmn(long m, long n, T z, T &res, T &res_jac) {
    int type;
    if (std::abs(z) <= 1) {
        type = 2;
    } else {
        type = 3;
    }

    special::assoc_legendre_p(n, m, type, z, res, res_jac);
}

template <typename T>
void lpmn(long m, long n, T z, T &res, T &res_jac, T &res_hess) {
    int type;
    if (std::abs(z) <= 1) {
        type = 2;
    } else {
        type = 3;
    }

    special::assoc_legendre_p(n, m, type, z, res, res_jac, res_hess);
}

template <typename T, typename OutputMat1>
void lpmn_all(T z, OutputMat1 res) {
    int type;
    if (std::abs(z) <= 1) {
        type = 2;
    } else {
        type = 3;
    }

    special::assoc_legendre_p_all(type, z, res);
}

template <typename T, typename OutputMat1, typename OutputMat2>
void lpmn_all(T z, OutputMat1 res, OutputMat2 res_jac) {
    int type;
    if (std::abs(z) <= 1) {
        type = 2;
    } else {
        type = 3;
    }

    special::assoc_legendre_p_all(type, z, res, res_jac);
}

template <typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void lpmn_all(T z, OutputMat1 res, OutputMat2 res_jac, OutputMat3 res_hess) {
    int type;
    if (std::abs(z) <= 1) {
        type = 2;
    } else {
        type = 3;
    }

    special::assoc_legendre_p_all(type, z, res, res_jac, res_hess);
}

template <typename T>
std::complex<T> clpmn(long m, long n, long type, std::complex<T> z) {
    return special::assoc_legendre_p(n, m, type, z);
}

template <typename T>
void clpmn(long m, long n, long type, T z, T &res, T &res_jac) {
    special::assoc_legendre_p(n, m, type, z, res, res_jac);
}

template <typename T>
void clpmn(long m, long n, long type, T z, T &res, T &res_jac, T &res_hess) {
    special::assoc_legendre_p(n, m, type, z, res, res_jac, res_hess);
}

template <typename T, typename OutputMat1>
void clpmn_all(long type, std::complex<T> z, OutputMat1 res) {
    special::assoc_legendre_p_all(type, z, res);
}

template <typename T, typename OutputMat1, typename OutputMat2>
void clpmn_all(long type, std::complex<T> z, OutputMat1 res, OutputMat2 res_jac) {
    special::assoc_legendre_p_all(type, z, res, res_jac);
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
