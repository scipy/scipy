#pragma once

#include "cephes/cbrt.h"

namespace xsf {

XSF_HOST_DEVICE inline double cbrt(double x) { return cephes::detail::cbrt(x); }

XSF_HOST_DEVICE inline float cbrt(float x) { return cbrt(static_cast<double>(x)); }

} // namespace xsf
