#pragma once

#include <xsf/third_party/kokkos/mdspan.hpp>

namespace special {

// Helper to wrap a 1D std::vector in a contiguous mdspan
template <typename T>
auto as_mdspan(std::vector<T> &vec) {
    return std::mdspan<T, std::dextents<ptrdiff_t, 1>>(vec.data(), vec.size());
}

template <typename T>
auto as_mdspan(const std::vector<T> &vec) {
    return std::mdspan<const T, std::dextents<ptrdiff_t, 1>>(vec.data(), vec.size());
}

} // namespace special
