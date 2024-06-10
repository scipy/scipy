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
const std::tuple<T &...> &tuple_assign(const std::tuple<T &...> &t, const initializer_tuple<T...> &other) {
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

} // namespace special
