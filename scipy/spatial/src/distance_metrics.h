#pragma once

#include <cmath>
#include "views.h"

struct MinkowskiDistance {
    double p_;

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y) const {
        T p = p_;
        for (intptr_t i = 0; i < x.shape[0]; ++i) {
            T p_dist = 0;
            for (intptr_t j = 0; j < x.shape[1]; ++j) {
                auto diff = std::abs(x(i, j) - y(i, j));
                p_dist += std::pow(diff, p);
            }
            out(i, 0) = std::pow(p_dist, 1 / p);
        }
    }

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y, StridedView2D<const T> w) const {
        T p = p_;
        for (intptr_t i = 0; i < x.shape[0]; ++i) {
            T p_dist = 0;
            for (intptr_t j = 0; j < x.shape[1]; ++j) {
                auto diff = std::abs(x(i, j) - y(i, j));
                p_dist += w(i, j) * std::pow(diff, p);
            }
            out(i, 0) = std::pow(p_dist, 1 / p);
        }
    }
};

struct EuclideanDistance {
    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y) const {
        for (intptr_t i = 0; i < x.shape[0]; ++i) {
            T sqdist = 0;
            for (intptr_t j = 0; j < x.shape[1]; ++j) {
                auto diff = std::abs(x(i, j) - y(i, j));
                sqdist += diff * diff;
            }
            out(i, 0) = std::sqrt(sqdist);
        }
    }

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y, StridedView2D<const T> w) const {
        for (intptr_t i = 0; i < x.shape[0]; ++i) {
            T sqdist = 0;
            for (intptr_t j = 0; j < x.shape[1]; ++j) {
                auto diff = std::abs(x(i, j) - y(i, j));
                sqdist += w(i, j) * (diff * diff);
            }
            out(i, 0) = std::sqrt(sqdist);
        }
    }
};

struct ChebyshevDistance {
    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y) const {
        for (intptr_t i = 0; i < x.shape[0]; ++i) {
            T dist = 0;
            for (intptr_t j = 0; j < x.shape[1]; ++j) {
                auto diff = std::abs(x(i, j) - y(i, j));
                if (diff > dist) {
                    dist = diff;
                }
            }
            out(i, 0) = dist;
        }
    }

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y, StridedView2D<const T> w) const {
        for (intptr_t i = 0; i < x.shape[0]; ++i) {
            T dist = 0;
            for (intptr_t j = 0; j < x.shape[1]; ++j) {
                auto diff = std::abs(x(i, j) - y(i, j));
                if (w(i, j) > 0 && diff > dist) {
                    dist = diff;
                }
            }
            out(i, 0) = dist;
        }
    }
};

struct CityBlockDistance {
    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y) const {
        if (x.strides[1] == 1 && y.strides[1] == 1) {
            for (intptr_t i = 0; i < x.shape[0]; ++i) {
                const T* x_data = x.data + i * x.strides[0];
                const T* y_data = y.data + i * y.strides[0];

                T dist = 0;
                for (intptr_t j = 0; j < x.shape[1]; ++j) {
                    dist += std::abs(x_data[j] - y_data[j]);
                }
                out(i, 0) = dist;
            }
        } else {
            for (intptr_t i = 0; i < x.shape[0]; ++i) {
                T dist = 0;
                for (intptr_t j = 0; j < x.shape[1]; ++j) {
                    dist += std::abs(x(i, j) - y(i, j));
                }
                out(i, 0) = dist;
            }
        }
    }

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y, StridedView2D<const T> w) const {
        for (intptr_t i = 0; i < x.shape[0]; ++i) {
            T dist = 0;
            for (intptr_t j = 0; j < x.shape[1]; ++j) {
                auto diff = std::abs(x(i, j) - y(i, j));
                dist += w(i, j) * diff;
            }
            out(i, 0) = dist;
        }
    }
};
