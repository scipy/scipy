#pragma once

#include <Python.h>
#include <cmath>

#include <xsf/cephes/pdtr.h>

namespace special {

inline void _legacy_cast_check(double x, double y) noexcept
{
    if (static_cast<int>(x) != x || static_cast<int>(y) != y) {
        PyGILState_STATE save = PyGILState_Ensure();
        // Purposefully ignore the raised PyError --- assume the ufunc will collect it.
        PyErr_WarnEx(PyExc_RuntimeWarning, "floating point number truncated to an integer", 1);
        PyGILState_Release(save);
    }
}

inline double pdtri_unsafe(double k, double y) noexcept
{
    if (std::isnan(k)) {
        return k;
    }
    _legacy_cast_check(k, 0);
    return xsf::cephes::pdtri(static_cast<int>(k), y);
}

} // namespace special
