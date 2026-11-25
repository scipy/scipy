#pragma once

#include <vector>

#include <xsf/mdspan.h>

namespace special {

// Helper to wrap a 1D std::vector in a contiguous mdspan
template <typename T>
auto as_mdspan(std::vector<T> &vec) {
    return xsf::cxx::mdspan<T, xsf::cxx::dextents<ptrdiff_t, 1>>(
        vec.data(), vec.size());
}

template <typename T>
auto as_mdspan(const std::vector<T> &vec) {
    return xsf::cxx::mdspan<const T, xsf::cxx::dextents<ptrdiff_t, 1>>(
        vec.data(), vec.size());
}

} // namespace special
