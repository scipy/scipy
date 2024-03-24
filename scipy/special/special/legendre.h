#pragma once

#include "config.h"
#include "mdspan.h"
#include "specfun/specfun.h"

namespace special {

template <typename T>
void lpn(T x, std::mdspan<T, std::dextents<int, 1>> pn, std::mdspan<T, std::dextents<int, 1>> pd);

template <>
void lpn(double x, std::mdspan<double, std::dextents<int, 1>> pn, std::mdspan<double, std::dextents<int, 1>> pd) {
    specfun::lpn(pn.size() - 1, x, pn.data_handle(), pd.data_handle());
}

template <>
void lpn(std::complex<double> z, std::mdspan<std::complex<double>, std::dextents<int, 1>> pn,
         std::mdspan<std::complex<double>, std::dextents<int, 1>> pd) {
    specfun::clpn(pn.size() - 1, z, pn.data_handle(), pd.data_handle());
}

} // namespace special
