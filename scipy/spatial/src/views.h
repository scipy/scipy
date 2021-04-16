#pragma once

#include <vector>
#include <array>
#include <cstdint>

struct ArrayDescriptor {
    ArrayDescriptor(int64_t ndim):
        ndim(ndim), shape(ndim, 1), strides(ndim, 0) {
    }

    int64_t ndim;
    int64_t element_size;
    std::vector<int64_t> shape, strides;
};

template <typename T>
struct StridedView2D {
    std::array<int64_t, 2> shape;
    std::array<int64_t, 2> strides;
    T* data;

    T& operator()(int64_t i, int64_t j) {
        return data[i * strides[0] + j * strides[1]];
    }
};
