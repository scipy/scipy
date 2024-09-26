#pragma once

#include "config.h"

namespace xsf {
namespace numbers {

    template <typename T>
    std::complex<T> i_v;

    template <>
    std::complex<float> i_v<float> = std::literals::complex_literals::operator""if(1.0L);

    template <>
    std::complex<double> i_v<double> = std::literals::complex_literals::operator""i(1.0L);

} // namespace numbers
} // namespace xsf
