#pragma once

#include <cmath>
#include "views.h"

#ifdef __GNUC__
#define ALWAYS_INLINE inline __attribute__((always_inline))
#define INLINE_LAMBDA __attribute__((always_inline))
#elif defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#define INLINE_LAMBDA
#else
#define ALWAYS_INLINE inline
#define INLINE_LAMBDA
#endif

struct Identity {
    template <typename T>
    ALWAYS_INLINE T operator() (T && val) const {
        return std::forward<T>(val);
    }
};

struct Plus {
    template <typename T>
    ALWAYS_INLINE T operator()(T a, T b) const {
        return a + b;
    }
};

// Helper to force a fixed bound loop to be completely unrolled
template <int unroll>
struct ForceUnroll{
    template <typename Func>
    ALWAYS_INLINE void operator()(const Func& f) const {
        ForceUnroll<unroll - 1>{}(f);
        f(unroll - 1);
    }
};

template <>
struct ForceUnroll<1> {
    template <typename Func>
    ALWAYS_INLINE void operator()(const Func& f) const {
        f(0);
    }
};

template <int ilp_factor=4,
          typename T,
          typename TransformFunc,
          typename ProjectFunc = Identity,
          typename ReduceFunc = Plus>
void transform_reduce_2d_(
    StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y,
    const TransformFunc& map,
    const ProjectFunc& project = Identity{},
    const ReduceFunc& reduce = Plus{}) {
    // Result type of calling map
    using AccumulateType = typename std::decay<decltype(
        map(std::declval<T>(), std::declval<T>()))>::type;
    intptr_t xs = x.strides[1], ys = y.strides[1];

    intptr_t i = 0;
    if (xs == 1 && ys == 1) {
        for (; i + (ilp_factor - 1) < x.shape[0]; i += ilp_factor) {
            const T* x_rows[ilp_factor];
            const T* y_rows[ilp_factor];
            ForceUnroll<ilp_factor>{}([&](int k) {
                x_rows[k] = &x(i + k, 0);
                y_rows[k] = &y(i + k, 0);
            });

            AccumulateType dist[ilp_factor] = {};
            for (intptr_t j = 0; j < x.shape[1]; ++j) {
                ForceUnroll<ilp_factor>{}([&](int k) {
                    auto val = map(x_rows[k][j], y_rows[k][j]);
                    dist[k] = reduce(dist[k], val);
                });
            }

            ForceUnroll<ilp_factor>{}([&](int k) {
                out(i + k, 0) = project(dist[k]);
            });
        }
    } else {
        for (; i + (ilp_factor - 1) < x.shape[0]; i += ilp_factor) {
            const T* x_rows[ilp_factor];
            const T* y_rows[ilp_factor];
            ForceUnroll<ilp_factor>{}([&](int k) {
                x_rows[k] = &x(i + k, 0);
                y_rows[k] = &y(i + k, 0);
            });

            AccumulateType dist[ilp_factor] = {};
            for (intptr_t j = 0; j < x.shape[1]; ++j) {
                auto x_offset = j * xs;
                auto y_offset = j * ys;
                ForceUnroll<ilp_factor>{}([&](int k) {
                    auto val = map(x_rows[k][x_offset], y_rows[k][y_offset]);
                    dist[k] = reduce(dist[k], val);
                });
            }

            ForceUnroll<ilp_factor>{}([&](int k) {
                out(i + k, 0) = project(dist[k]);
            });
        }
    }
    for (; i < x.shape[0]; ++i) {
        const T* x_row = &x(i, 0);
        const T* y_row = &y(i, 0);
        AccumulateType dist = {};
        for (intptr_t j = 0; j < x.shape[1]; ++j) {
            auto val = map(x_row[j * xs], y_row[j * ys]);
            dist = reduce(dist, val);
        }
        out(i, 0) = project(dist);
    }
}

template <int ilp_factor=2, typename T,
          typename TransformFunc,
          typename ProjectFunc = Identity,
          typename ReduceFunc = Plus>
void transform_reduce_2d_(
    StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y,
    StridedView2D<const T> w, const TransformFunc& map,
    const ProjectFunc& project = Identity{},
    const ReduceFunc& reduce = Plus{}) {
    intptr_t i = 0;
    intptr_t xs = x.strides[1], ys = y.strides[1], ws = w.strides[1];
    // Result type of calling map
    using AccumulateType = typename std::decay<decltype(
        map(std::declval<T>(), std::declval<T>(), std::declval<T>()))>::type;

    for (; i + (ilp_factor - 1) < x.shape[0]; i += ilp_factor) {
        const T* x_rows[ilp_factor];
        const T* y_rows[ilp_factor];
        const T* w_rows[ilp_factor];
        ForceUnroll<ilp_factor>{}([&](int k) {
            x_rows[k] = &x(i + k, 0);
            y_rows[k] = &y(i + k, 0);
            w_rows[k] = &w(i + k, 0);
        });

        AccumulateType dist[ilp_factor] = {};
        for (intptr_t j = 0; j < x.shape[1]; ++j) {
            ForceUnroll<ilp_factor>{}([&](int k) {
                auto val = map(x_rows[k][j * xs], y_rows[k][j * ys], w_rows[k][j * ws]);
                dist[k] = reduce(dist[k], val);
            });
        }

        ForceUnroll<ilp_factor>{}([&](int k) {
            out(i + k, 0) = project(dist[k]);
        });
    }
    for (; i < x.shape[0]; ++i) {
        const T* x_row = &x(i, 0);
        const T* y_row = &y(i, 0);
        const T* w_row = &w(i, 0);
        AccumulateType dist = {};
        for (intptr_t j = 0; j < x.shape[1]; ++j) {
            auto val = map(x_row[j * xs], y_row[j * ys], w_row[j * ws]);
            dist = reduce(dist, val);
        }
        out(i, 0) = project(dist);
    }
}

struct MinkowskiDistance {
    double p_;

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y) const {
        const T p = static_cast<T>(p_);
        const T invp = static_cast<T>(1.0 / p_);
        transform_reduce_2d_(out, x, y, [p](T x, T y) INLINE_LAMBDA {
            auto diff = std::abs(x - y);
            return std::pow(diff, p);
        },
        [invp](T x) { return std::pow(x, invp); });
    }

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y, StridedView2D<const T> w) const {
        const T p = static_cast<T>(p_);
        const T invp = static_cast<T>(1.0 / p_);
        transform_reduce_2d_(out, x, y, w, [p](T x, T y, T w) INLINE_LAMBDA {
            auto diff = std::abs(x - y);
            return w * std::pow(diff, p);
        },
        [invp](T x) { return std::pow(x, invp); });
    }
};

struct EuclideanDistance {
    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y) const {
        transform_reduce_2d_(out, x, y, [](T x, T y) INLINE_LAMBDA {
            auto diff = std::abs(x - y);
            return diff * diff;
        },
        [](T x) { return std::sqrt(x); });
    }

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y, StridedView2D<const T> w) const {
        transform_reduce_2d_(out, x, y, w, [](T x, T y, T w) INLINE_LAMBDA {
            auto diff = std::abs(x - y);
            return w * (diff * diff);
        },
        [](T x) { return std::sqrt(x); });
    }
};

struct ChebyshevDistance {
    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y) const {
        transform_reduce_2d_(out, x, y, [](T x, T y) INLINE_LAMBDA {
            return std::abs(x - y);
        },
        Identity{},
        [](T x, T y) { return std::max(x, y); });
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
        transform_reduce_2d_<2>(out, x, y, [](T x, T y) INLINE_LAMBDA {
            return std::abs(x - y);
        });
    }

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y, StridedView2D<const T> w) const {
        transform_reduce_2d_(out, x, y, w, [](T x, T y, T w) INLINE_LAMBDA {
            return w * std::abs(x - y);
        });
    }
};

struct SquareEuclideanDistance {
    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y) const {
        transform_reduce_2d_(out, x, y, [](T x, T y) INLINE_LAMBDA {
            auto diff = x - y;
            return diff * diff;
        });
    }

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y, StridedView2D<const T> w) const {
        transform_reduce_2d_(out, x, y, w, [](T x, T y, T w) INLINE_LAMBDA {
            auto diff = x - y;
            return w * diff * diff;
        });
    }
};

struct BraycurtisDistance {
    template <typename T>
    struct Acc {
        Acc(): diff(0), sum(0) {}
        T diff, sum;
    };

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y) const {
        // dist = abs(x - y).sum() / abs(x + y).sum()
        transform_reduce_2d_<2>(out, x, y, [](T x, T y) INLINE_LAMBDA {
            Acc<T> acc;
            acc.diff = std::abs(x - y);
            acc.sum = std::abs(x + y);
            return acc;
        },
        [](const Acc<T>& acc) INLINE_LAMBDA {
            return acc.diff / acc.sum;
        },
        [](const Acc<T>& a, const Acc<T>& b) INLINE_LAMBDA {
            Acc<T> acc;
            acc.diff = a.diff + b.diff;
            acc.sum = a.sum + b.sum;
            return acc;
        });
    }

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y, StridedView2D<const T> w) const {
        // dist = (w * abs(x - y)).sum() / (w * abs(x + y)).sum()
        transform_reduce_2d_(out, x, y, w, [](T x, T y, T w) INLINE_LAMBDA {
            Acc<T> acc;
            acc.diff = w * std::abs(x - y);
            acc.sum = w * std::abs(x + y);
            return acc;
        },
        [](const Acc<T>& acc) INLINE_LAMBDA {
            return acc.diff / acc.sum;
        },
        [](const Acc<T>& a, const Acc<T>& b) INLINE_LAMBDA {
            Acc<T> acc;
            acc.diff = a.diff + b.diff;
            acc.sum = a.sum + b.sum;
            return acc;
        });
    }
};

struct CanberraDistance {
    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y) const {
        // dist = (abs(x - y) / (abs(x) + abs(y))).sum()
        transform_reduce_2d_<2>(out, x, y, [](T x, T y) INLINE_LAMBDA {
            auto num = std::abs(x - y);
            auto denom = std::abs(x) + std::abs(y);
            // branchless replacement for (denom == 0) ? 0 : num / denom;
            return num / (denom + (denom == 0));
        });
    }

    template <typename T>
    void operator()(StridedView2D<T> out, StridedView2D<const T> x, StridedView2D<const T> y, StridedView2D<const T> w) const {
        // dist = (w * abs(x - y) / (abs(x) + abs(y))).sum()
        transform_reduce_2d_(out, x, y, w, [](T x, T y, T w) INLINE_LAMBDA {
            auto num = w * std::abs(x - y);
            auto denom = std::abs(x) + std::abs(y);
            // branchless replacement for (denom == 0) ? 0 : num / denom;
            return num / (denom + (denom == 0));
        });
    }
};
