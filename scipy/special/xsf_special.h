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

template <typename T>
T legendre_p(long long int n, T z) {
    return xsf::legendre_p(n, z);
}

template <typename T>
void legendre_p(long long int n, T z, T &res, T &res_jac) {
    dual_assign_grad(xsf::legendre_p(n, xsf::make_dual<1>(z)), std::tie(res, res_jac));
}

template <typename T>
void legendre_p(long long int n, T z, T &res, T &res_jac, T &res_hess) {
    dual_assign_grad(xsf::legendre_p(n, xsf::make_dual<2>(z)), std::tie(res, res_jac, res_hess));
}

template <typename T, typename OutputVec1>
void legendre_p_all(T z, OutputVec1 res) {
    xsf::legendre_p_all(xsf::make_dual<0>(z), std::tie(res));
}

template <typename T, typename OutputVec1, typename OutputVec2>
void legendre_p_all(T z, OutputVec1 res, OutputVec2 res_jac) {
    xsf::legendre_p_all(xsf::make_dual<1>(z), std::tie(res, res_jac));
}

template <typename T, typename OutputVec1, typename OutputVec2, typename OutputVec3>
void legendre_p_all(T z, OutputVec1 res, OutputVec2 res_jac, OutputVec3 res_hess) {
    xsf::legendre_p_all(xsf::make_dual<2>(z), std::tie(res, res_jac, res_hess));
}

using xsf::assoc_legendre_norm;
using xsf::assoc_legendre_unnorm;

template <typename NormPolicy, typename T>
T assoc_legendre_p(NormPolicy norm, long long int n, long long int m, T z, long long int branch_cut) {
    return xsf::assoc_legendre_p(norm, n, m, z, branch_cut);
}

template <typename NormPolicy, typename T>
void assoc_legendre_p(NormPolicy norm, long long int n, long long int m, T z, long long int branch_cut, T &res,
                      T &res_jac) {
    xsf::dual_assign_grad(xsf::assoc_legendre_p(norm, n, m, xsf::make_dual<1>(z), branch_cut), std::tie(res, res_jac));
}

template <typename NormPolicy, typename T>
void assoc_legendre_p(NormPolicy norm, long long int n, long long int m, T z, long long int branch_cut, T &res,
                      T &res_jac, T &res_hess) {
    xsf::dual_assign_grad(xsf::assoc_legendre_p(norm, n, m, xsf::make_dual<2>(z), branch_cut),
                          std::tie(res, res_jac, res_hess));
}

template <typename NormPolicy, typename T, typename OutputMat1>
void assoc_legendre_p_all(NormPolicy norm, T z, long long int branch_cut, OutputMat1 res) {
    xsf::assoc_legendre_p_all(norm, xsf::make_dual<0>(z), branch_cut, std::tie(res));
}

template <typename NormPolicy, typename T, typename OutputMat1, typename OutputMat2>
void assoc_legendre_p_all(NormPolicy norm, T z, long long int branch_cut, OutputMat1 res, OutputMat2 res_jac) {
    xsf::assoc_legendre_p_all(norm, xsf::make_dual<1>(z), branch_cut, std::tie(res, res_jac));
}

template <typename NormPolicy, typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void assoc_legendre_p_all(NormPolicy norm, T z, long long int branch_cut, OutputMat1 res, OutputMat2 res_jac,
                          OutputMat3 res_hess) {
    xsf::assoc_legendre_p_all(norm, xsf::make_dual<2>(z), branch_cut, std::tie(res, res_jac, res_hess));
}

template <typename T>
T sph_legendre_p(long long int n, long long int m, T theta) {
    return xsf::sph_legendre_p(n, m, xsf::make_dual<0>(theta)).front();
}

template <typename T>
void sph_legendre_p(long long int n, long long int m, T theta, T &res, T &res_jac) {
    xsf::dual_assign_grad(xsf::sph_legendre_p(n, m, xsf::make_dual<1>(theta)), std::tie(res, res_jac));
}

template <typename T>
void sph_legendre_p(long long int n, long long int m, T theta, T &res, T &res_jac, T &res_hess) {
    xsf::dual_assign_grad(xsf::sph_legendre_p(n, m, xsf::make_dual<2>(theta)), std::tie(res, res_jac, res_hess));
}

template <typename T, typename OutputMat1>
void sph_legendre_p_all(T theta, OutputMat1 res) {
    xsf::sph_legendre_p_all(xsf::make_dual<0>(theta), std::tie(res));
}

template <typename T, typename OutputMat1, typename OutputMat2>
void sph_legendre_p_all(T theta, OutputMat1 res, OutputMat2 res_jac) {
    xsf::sph_legendre_p_all(xsf::make_dual<1>(theta), std::tie(res, res_jac));
}

template <typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void sph_legendre_p_all(T theta, OutputMat1 res, OutputMat2 res_jac, OutputMat3 res_hess) {
    xsf::sph_legendre_p_all(xsf::make_dual<2>(theta), std::tie(res, res_jac, res_hess));
}

template <typename T>
std::complex<T> sph_harm_y(long long int n, long long int m, T theta, T phi) {
    xsf::dual<xsf::dual<T, 0>, 0> theta_dual;
    theta_dual[0][0] = theta;

    xsf::dual<xsf::dual<T, 0>, 0> phi_dual;
    phi_dual[0][0] = phi;

    std::complex<T> res;
    dual_assign_grad(xsf::sph_harm_y(n, m, theta_dual, phi_dual), std::tie(res));

    return res;
}

template <typename T>
void sph_harm_y(long long int n, long long int m, T theta, T phi, std::complex<T> &res, std::complex<T> (&res_jac)[2]) {
    xsf::dual<xsf::dual<T, 1>, 1> theta_dual;
    theta_dual[0][0] = theta;
    theta_dual[0][1] = 0;
    theta_dual[1][0] = 1;
    theta_dual[1][1] = 0;

    xsf::dual<xsf::dual<T, 1>, 1> phi_dual;
    phi_dual[0][0] = phi;
    phi_dual[0][1] = 1;
    phi_dual[1][0] = 0;
    phi_dual[1][1] = 0;

    dual_assign_grad(xsf::sph_harm_y(n, m, theta_dual, phi_dual), std::tie(res, res_jac));
}

template <typename T>
void sph_harm_y(long long int n, long long int m, T theta, T phi, std::complex<T> &res, std::complex<T> (&res_jac)[2],
                std::complex<T> (&res_hess)[2][2]) {
    xsf::dual<xsf::dual<T, 2>, 2> theta_dual;
    theta_dual[0][0] = theta;
    theta_dual[0][1] = 0;
    theta_dual[0][2] = 0;
    theta_dual[1][0] = 1;
    theta_dual[1][1] = 0;
    theta_dual[1][2] = 0;
    theta_dual[2][0] = 0;
    theta_dual[2][1] = 0;
    theta_dual[2][2] = 0;

    xsf::dual<xsf::dual<T, 2>, 2> phi_dual;
    phi_dual[0][0] = phi;
    phi_dual[0][1] = 1;
    phi_dual[0][2] = 0;
    phi_dual[1][0] = 0;
    phi_dual[1][1] = 0;
    phi_dual[1][2] = 0;
    phi_dual[2][0] = 0;
    phi_dual[2][1] = 0;
    phi_dual[2][2] = 0;

    dual_assign_grad(xsf::sph_harm_y(n, m, theta_dual, phi_dual), std::tie(res, res_jac, res_hess));
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

    return sph_harm_y(n, m, phi, theta);
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

template <typename T, typename OutputMat1>
void sph_harm_y_all(T theta, T phi, OutputMat1 res) {
    xsf::dual<xsf::dual<T, 0>, 0> theta_dual;
    theta_dual[0][0] = theta;

    xsf::dual<xsf::dual<T, 0>, 0> phi_dual;
    phi_dual[0][0] = phi;

    xsf::sph_harm_y_all(theta_dual, phi_dual, std::tie(res));
}

template <typename T, typename OutputMat1, typename OutputMat2>
void sph_harm_y_all(T theta, T phi, OutputMat1 res, OutputMat2 res_jac) {
    xsf::dual<xsf::dual<T, 1>, 1> theta_dual;
    theta_dual[0][0] = theta;
    theta_dual[0][1] = 0;
    theta_dual[1][0] = 1;
    theta_dual[1][1] = 0;

    xsf::dual<xsf::dual<T, 1>, 1> phi_dual;
    phi_dual[0][0] = phi;
    phi_dual[0][1] = 1;
    phi_dual[1][0] = 0;
    phi_dual[1][1] = 0;

    xsf::sph_harm_y_all(theta_dual, phi_dual, std::tie(res, res_jac));
}

template <typename T, typename OutputMat1, typename OutputMat2, typename OutputMat3>
void sph_harm_y_all(T theta, T phi, OutputMat1 res, OutputMat2 res_jac, OutputMat3 res_hess) {
    xsf::dual<xsf::dual<T, 2>, 2> theta_dual;
    theta_dual[0][0] = theta;
    theta_dual[0][1] = 0;
    theta_dual[0][2] = 0;
    theta_dual[1][0] = 1;
    theta_dual[1][1] = 0;
    theta_dual[1][2] = 0;
    theta_dual[2][0] = 0;
    theta_dual[2][1] = 0;
    theta_dual[2][2] = 0;

    xsf::dual<xsf::dual<T, 2>, 2> phi_dual;
    phi_dual[0][0] = phi;
    phi_dual[0][1] = 1;
    phi_dual[0][2] = 0;
    phi_dual[1][0] = 0;
    phi_dual[1][1] = 0;
    phi_dual[1][2] = 0;
    phi_dual[2][0] = 0;
    phi_dual[2][1] = 0;
    phi_dual[2][2] = 0;

    xsf::sph_harm_y_all(theta_dual, phi_dual, std::tie(res, res_jac, res_hess));
}

} // namespace
