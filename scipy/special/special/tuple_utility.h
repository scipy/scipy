#pragma once

#include <tuple>

#include "config.h"

namespace special {

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
struct initializer_tuple {
    std::tuple<T...> underlying;

    initializer_tuple(const T &...args) {
        std::apply([&args...](auto &...other_args) { (assign(other_args, args), ...); }, underlying);
    }
};

template <typename... T>
std::tuple<T &...> &tuple_assign(std::tuple<T &...> &t, const initializer_tuple<T...> &other) {
    std::apply(
        [&t](const auto &...other_args) {
            std::apply([&other_args...](auto &...args) { (assign(args, other_args), ...); }, t);
        },
        other.underlying
    );

    return t;
}

template <typename... T, typename Arg>
decltype(auto) tuple_access_each(std::tuple<T...> &t, Arg n) {
    return std::apply([&n](auto &...args) { return std::tie(args[n]...); }, t);
}

template <typename... T, typename... Args>
std::tuple<std::result_of_t<T(Args...)>...> tuple_call_each(std::tuple<T...> &t, Args &&...args) {
    return std::apply([&args...](auto &...elements) { return std::tie(elements(std::forward<Args>(args)...)...); }, t);
}

template <typename... T, typename Arg>
std::tuple<T...> &tuple_fill_each(std::tuple<T...> &t, Arg arg) {
    std::apply([arg](auto &...args) { (std::fill(std::begin(args), std::end(args), arg), ...); }, t);

    return t;
}

template <typename... T, typename Func>
void tuple_for_each(std::tuple<T...> &t, Func f) {
    std::apply([f](auto &...args) { (f(args), ...); }, t);
}

template <typename... T>
std::tuple<T &...> tuple_ref_each(std::tuple<T...> &t) {
    return std::apply([](auto &...args) { return std::tie(args...); }, t);
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
