#pragma once

#include <tuple>

#include "config.h"

namespace special {

template <typename T, size_t N>
class grad;

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

} // namespace special

namespace std {

template <size_t I, typename T, size_t N>
struct tuple_element<I, special::grad<T, N>> {
    //    using type = T;
    using type = special::standardize_t<T>;
};

template <typename T, size_t N>
struct tuple_size<special::grad<T, N>> : std::integral_constant<size_t, N + 1> {};

} // namespace std

namespace special {
namespace detail {

    template <typename T>
    void assign(T &dst, const T &src) {
        dst = src;
    }

    template <typename T, size_t K>
    void assign(T (&dst)[K], const T (&src)[K]) {
        for (size_t k = 0; k < K; ++k) {
            dst[k] = src[k];
        }
    }

    template <typename T, size_t K>
    void assign(std::array<T, K> &dst, const std::array<T, K> &src) {
        for (size_t k = 0; k < K; ++k) {
            dst[k] = src[k];
        }
    }

    template <typename Indices, typename T>
    class base_grad;

    template <size_t... I, typename T>
    class base_grad<std::index_sequence<I...>, T> {
      public:
        using tuple_type = std::tuple<std::tuple_element_t<I, grad<T, sizeof...(I) - 1>>...>;

      protected:
        tuple_type m_underlying;

      public:
        base_grad() = default;

        base_grad(const std::tuple_element_t<I, tuple_type> &...args) : m_underlying(args...) {}

        template <typename... U>
        base_grad(const std::tuple<U...> &other) : m_underlying(other) {}

        base_grad(const base_grad &other) = default;

        base_grad(base_grad &&other) = default;

        tuple_type &underlying_tuple() { return m_underlying; }

        const tuple_type &underlying_tuple() const { return m_underlying; }

        base_grad &operator=(const base_grad &other) = default;

        base_grad &operator=(base_grad &&other) = default;
    };

} // namespace detail

template <typename T, size_t N>
class grad : public detail::base_grad<std::make_index_sequence<N + 1>, T> {
  public:
    using base = detail::base_grad<std::make_index_sequence<N + 1>, T>;
    using tuple_type = typename base::tuple_type;

  public:
    grad(const base &other) : base(other) {}

  public:
    using base::base;

    grad<T &, N> refs() {
        return std::apply([](auto &...args) { return std::tie(args...); }, base::underlying_tuple());
    }

    grad<const T &, N> refs() const { return crefs(); }

    grad<const T &, N> crefs() const {
        return std::apply([](const auto &...args) { return std::tie(args...); }, base::underlying_tuple());
    }
};

template <typename T, size_t N>
class grad<T &, N> : public detail::base_grad<std::make_index_sequence<N + 1>, T &> {
  public:
    using base = detail::base_grad<std::make_index_sequence<N + 1>, T &>;
    using tuple_type = typename base::tuple_type;

  public:
    grad(const base &other) : base(other) {}

  public:
    using base::base;

    grad<T &, N> refs() {
        return std::apply([](auto &...args) { return std::tie(args...); }, base::underlying_tuple());
    }

    grad<const T &, N> refs() const { return crefs(); }

    grad<const T &, N> crefs() const {
        return std::apply([](const auto &...args) { return std::tie(args...); }, base::underlying_tuple());
    }

    grad &operator=(const grad<T, N> &other) {
        base::underlying_tuple() = other.underlying_tuple();

        return *this;
    }
};

template <typename T0, typename... T>
grad(T0 &, T &...) -> grad<T0 &, sizeof...(T)>;

template <size_t I, typename T, size_t N>
std::tuple_element_t<I, grad<T, N>> &get(grad<T, N> &t) {
    return std::get<I>(t.underlying_tuple());
}

template <size_t I, typename T, size_t N>
const std::tuple_element_t<I, grad<T, N>> &get(const grad<T, N> &t) {
    return std::get<I>(t.underlying_tuple());
}

template <typename T, size_t K>
void dot(grad<T (&)[K], 0> x, const grad<T (&)[K], 0> y, grad<T &, 0> res) {
    const auto &[x0] = x;
    const auto &[y0] = y;
    auto &[res0] = res;

    res0 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
    }
}

template <typename T, size_t K>
void dot(grad<T (&)[K], 1> x, grad<T (&)[K], 1> y, grad<T &, 1> res) {
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
void dot(grad<T (&)[K], 2> x, grad<T (&)[K], 2> y, grad<T &, 2> res) {
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
