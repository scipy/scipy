#pragma once

#include "cephes/ellpk.h"
#include "config.h"

namespace special {

SPECFUN_HOST_DEVICE inline double ellipk(double m) { return cephes::ellpk(1.0 - m); }
} // namespace special
