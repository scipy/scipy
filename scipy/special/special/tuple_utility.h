#pragma once

#include <tuple>

#include "config.h"

namespace special {
namespace detail {

    template <typename T>
    void assign(T &dst, const T &src) {
        dst = src;
    }

    template <typename T, size_t N>
    void assign(T (&dst)[N], const T (&src)[N]) {
        for (size_t i = 0; i < N; ++i) {
            dst[i] = src[i];
        }
    }

    template <typename... T>
    struct tuple_wrapper {
        std::tuple<T...> underlying;

        tuple_wrapper(const T &...args) {
            std::apply([&args...](auto &...other_args) { (assign(other_args, args), ...); }, underlying);
        }
    };

    template <typename... T>
    struct tuple_assigner {
        std::tuple<T...> &assignee;

        tuple_assigner(std::tuple<T...> &assignee) : assignee(assignee) {}

        tuple_assigner &operator=(const tuple_wrapper<std::remove_reference_t<T>...> &other) {
            std::apply(
                [this](const auto &...other_args) {
                    std::apply([&other_args...](auto &...args) { (assign(args, other_args), ...); }, assignee);
                },
                other.underlying
            );

            return *this;
        }
    };

} // namespace detail

template <typename... T>
std::tuple<T &...> apply_tie(std::tuple<T...> &t) {
    return std::apply([](auto &...args) { return std::tie(args...); }, t);
}

template <typename... T>
detail::tuple_assigner<T...> tuple_assigner(std::tuple<T...> &t) {
    return t;
}

template <typename T, size_t K>
void dot(std::tuple<T (&)[K]> x, const std::tuple<T (&)[K]> y, std::tuple<T &> res) {
    const auto &[x0] = x;
    const auto &[y0] = y;
    auto &[res0] = res;

    res0 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
    }
}

template <typename T, size_t K>
void dot(std::tuple<T (&)[K], T (&)[K]> x, std::tuple<T (&)[K], T (&)[K]> y, std::tuple<T &, T &> res) {
    const auto &[x0, x1] = x;
    const auto &[y0, y1] = y;
    auto &[res0, res1] = res;

    res0 = 0;
    res1 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
        res1 += x0[k] * y1[k] + x1[k] * y0[k];
    }
}

template <typename T, size_t K>
void dot(
    std::tuple<T (&)[K], T (&)[K], T (&)[K]> x, std::tuple<T (&)[K], T (&)[K], T (&)[K]> y,
    std::tuple<T &, T &, T &> res
) {
    const auto &[x0, x1, x2] = x;
    const auto &[y0, y1, y2] = y;
    auto &[res0, res1, res2] = res;

    res0 = 0;
    res1 = 0;
    res2 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
        res1 += x0[k] * y1[k] + x1[k] * y0[k];
        res2 += x0[k] * y2[k] + T(2) * x1[k] * y1[k] + x2[k] * y0[k];
    }
}

} // namespace special
