#pragma once

#include <vector>
#include <array>
#include <cstdint>

struct ArrayDescriptor {
    ArrayDescriptor(intptr_t ndim):
        ndim(ndim), shape(ndim, 1), strides(ndim, 0) {
    }

    intptr_t ndim;
    intptr_t element_size;
    std::vector<intptr_t> shape, strides;
};

template <typename T>
struct StridedView2D {
    std::array<intptr_t, 2> shape;
    std::array<intptr_t, 2> strides;
    T* data;

    T& operator()(intptr_t i, intptr_t j) {
        return data[i * strides[0] + j * strides[1]];
    }
};
