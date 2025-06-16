#pragma once

#include "Python.h"

#include <xsf/bessel.h>
#include <xsf/sph_harm.h>

// This header exists to add behaviors to special functions from the xsf library,
// either because they involve some Python-specific features or because there are
// legacy behaviors that need to be preserved

namespace {

template <typename T>
std::complex<T> sph_harm(long long int m, long long int n, T theta, T phi) {
    PyGILState_STATE gstate = PyGILState_Ensure();
    PyErr_WarnEx(
        PyExc_DeprecationWarning,
        "`scipy.special.sph_harm` is deprecated as of SciPy 1.15.0 and will be "
        "removed in SciPy 1.17.0. Please use `scipy.special.sph_harm_y` instead.",
        1
    );
    PyGILState_Release(gstate);
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

} // namespace
