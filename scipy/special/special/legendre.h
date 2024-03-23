#pragma once

#include "config.h"
#include "specfun/specfun.h"

namespace special {

inline void lpn(double x, std::span<double> pn, std::span<double> pd) {
    specfun::lpn(pn.size() - 1, x, pn.data(), pd.data());
}

} // namespace special