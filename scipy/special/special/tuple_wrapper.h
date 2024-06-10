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
class tuple_wrapper {
  public:
    using tuple_type = std::tuple<T...>;

  protected:
    tuple_type m_underlying;

  public:
    tuple_wrapper() = default;

    tuple_wrapper(const T &...args) {
        std::apply([&args...](auto &...other_args) { (assign(other_args, args), ...); }, m_underlying);
    }

    template <typename... U>
    tuple_wrapper(const std::tuple<U...> &other) : m_underlying(other) {}

    tuple_wrapper(const tuple_wrapper &other) = default;

    tuple_wrapper(tuple_wrapper &&other) = default;

    tuple_type &underlying_tuple() { return m_underlying; }

    const tuple_type &underlying_tuple() const { return m_underlying; }

    tuple_wrapper &operator=(const tuple_wrapper &other) = default;

    tuple_wrapper &operator=(tuple_wrapper &&other) = default;
};

template <typename... T>
class tuple_wrapper<T &...> {
  public:
    using tuple_type = std::tuple<T &...>;

  protected:
    tuple_type m_underlying;

  public:
    tuple_wrapper() = default;

    tuple_wrapper(T &...args) : m_underlying(args...) {}

    tuple_wrapper(const tuple_type &other) : m_underlying(other) {}

    template <typename... U>
    tuple_wrapper(const std::tuple<U...> &other) : m_underlying(other) {}

    tuple_wrapper(const tuple_wrapper &other) = default;

    tuple_wrapper(tuple_wrapper &&other) = default;

    tuple_type &underlying_tuple() { return m_underlying; }

    const tuple_type &underlying_tuple() const { return m_underlying; }

    tuple_wrapper &operator=(const tuple_wrapper &other) = default;

    tuple_wrapper &operator=(tuple_wrapper &&other) = default;

    tuple_wrapper &operator=(const tuple_wrapper<T...> &other) {
        std::apply(
            [&other](auto &...args) {
                std::apply(
                    [&args...](const auto &...other_args) { (assign(args, other_args), ...); }, other.underlying_tuple()
                );
            },
            m_underlying
        );

        return *this;
    }
};

template <typename... T>
tuple_wrapper<T &...> apply_tie(tuple_wrapper<T...> &t) {
    return std::apply([](auto &...args) { return std::tie(args...); }, t.underlying_tuple());
}

template <typename... T>
std::tuple<T &...> apply_tie(std::tuple<T...> &t) {
    return std::apply([](auto &...args) { return std::tie(args...); }, t);
}

template <size_t I, typename... T>
std::tuple_element_t<I, tuple_wrapper<T...>> &get(tuple_wrapper<T...> &t) {
    return std::get<I>(t.underlying_tuple());
}

template <size_t I, typename... T>
const std::tuple_element_t<I, tuple_wrapper<T...>> &get(const tuple_wrapper<T...> &t) {
    return std::get<I>(t.underlying_tuple());
}

template <typename T, size_t K>
void dot(tuple_wrapper<T (&)[K]> x, const tuple_wrapper<T (&)[K]> y, tuple_wrapper<T &> res) {
    const auto &[x0] = x;
    const auto &[y0] = y;
    auto &[res0] = res;

    res0 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
    }
}

template <typename T, size_t K>
void dot(tuple_wrapper<T (&)[K], T (&)[K]> x, tuple_wrapper<T (&)[K], T (&)[K]> y, tuple_wrapper<T &, T &> res) {
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
    tuple_wrapper<T (&)[K], T (&)[K], T (&)[K]> x, tuple_wrapper<T (&)[K], T (&)[K], T (&)[K]> y,
    tuple_wrapper<T &, T &, T &> res
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

namespace std {

template <size_t I, typename... T>
struct tuple_element<I, special::tuple_wrapper<T...>> {
    using type = std::tuple_element_t<I, typename special::tuple_wrapper<T...>::tuple_type>;
};

template <typename... T>
struct tuple_size<special::tuple_wrapper<T...>> : std::integral_constant<size_t, sizeof...(T)> {};

} // namespace std
