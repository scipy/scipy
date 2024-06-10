#pragma once

#include <tuple>

#include "config.h"

namespace special {

template <typename T>
struct standardize {
    static_assert(!std::is_array_v<T>, "asf");
    using type = T;
};

template <typename T, size_t N>
struct standardize<T[N]> {
    using type = std::array<T, N>;
};

template <typename T, size_t N>
struct standardize<const T[N]> {
    using type = std::array<const T, N>;
};

template <typename T, size_t N>
struct standardize<T (&)[N]> {
    using type = std::array<T, N> &;
};

template <typename T, size_t N>
struct standardize<const T (&)[N]> {
    using type = const std::array<T, N> &;
};

template <typename T>
using standardize_t = typename standardize<T>::type;

template <typename... T>
class grad {
  public:
    using tuple_type = std::tuple<standardize_t<T>...>;

  protected:
    tuple_type m_underlying;

  public:
    grad() = default;

    grad(const standardize_t<T> &...args) : m_underlying(args...) {}

    template <typename... U>
    grad(const std::tuple<U...> &other) : m_underlying(other) {}

    grad(const grad &other) = default;

    grad(grad &&other) = default;

    tuple_type &underlying_tuple() { return m_underlying; }

    const tuple_type &underlying_tuple() const { return m_underlying; }

    grad &operator=(const grad &other) = default;

    grad &operator=(grad &&other) = default;

    grad<T &...> refs() {
        return std::apply([](auto &...args) { return std::tie(args...); }, underlying_tuple());
    }

    grad<const T &...> refs() const { return crefs(); }

    grad<const T &...> crefs() const {
        return std::apply([](const auto &...args) { return std::tie(args...); }, underlying_tuple());
    }
};

template <typename... T>
class grad<T &...> {
  public:
    using tuple_type = std::tuple<standardize_t<T &>...>;

  protected:
    tuple_type m_underlying;

  public:
    grad() = default;

    grad(const standardize_t<T &> &...args) : m_underlying(args...) {}

    template <typename... U>
    grad(const std::tuple<U...> &other) : m_underlying(other) {}

    grad(const grad &other) = default;

    grad(grad &&other) = default;

    tuple_type &underlying_tuple() { return m_underlying; }

    const tuple_type &underlying_tuple() const { return m_underlying; }

    grad &operator=(const grad &other) = default;

    grad &operator=(grad &&other) = default;

    grad<T &...> refs() {
        return std::apply([](auto &...args) { return std::tie(args...); }, underlying_tuple());
    }

    grad<const T &...> refs() const { return crefs(); }

    grad<const T &...> crefs() const {
        return std::apply([](const auto &...args) { return std::tie(args...); }, underlying_tuple());
    }

    grad &operator=(const grad<T...> &other) {
        underlying_tuple() = other.underlying_tuple();

        return *this;
    }
};

template <size_t I, typename... T>
std::tuple_element_t<I, grad<T...>> &get(grad<T...> &t) {
    return std::get<I>(t.underlying_tuple());
}

template <size_t I, typename... T>
const std::tuple_element_t<I, grad<T...>> &get(const grad<T...> &t) {
    return std::get<I>(t.underlying_tuple());
}

template <typename T, size_t K>
void dot(grad<T (&)[K]> x, const grad<T (&)[K]> y, grad<T &> res) {
    const auto &[x0] = x;
    const auto &[y0] = y;
    auto &[res0] = res;

    res0 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
    }
}

template <typename T, size_t K>
void dot(grad<T (&)[K], T (&)[K]> x, grad<T (&)[K], T (&)[K]> y, grad<T &, T &> res) {
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
void dot(grad<T (&)[K], T (&)[K], T (&)[K]> x, grad<T (&)[K], T (&)[K], T (&)[K]> y, grad<T &, T &, T &> res) {
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
struct tuple_element<I, special::grad<T...>> {
    using type = std::tuple_element_t<I, typename special::grad<T...>::tuple_type>;
};

template <typename... T>
struct tuple_size<special::grad<T...>> : std::integral_constant<size_t, sizeof...(T)> {};

} // namespace std
