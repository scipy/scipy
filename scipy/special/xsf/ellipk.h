#pragma once

#include "cephes/ellpk.h"
#include "config.h"

namespace xsf {

XSF_HOST_DEVICE inline double ellipk(double m) { return cephes::ellpk(1.0 - m); }

} // namespace xsf
