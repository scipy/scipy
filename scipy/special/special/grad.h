#pragma once

#include <tuple>

#include "config.h"

namespace special {

template <typename T, size_t N>
class grad;

} // namespace special

namespace std {

template <size_t I, typename T, size_t N>
struct tuple_element<I, special::grad<T, N>> {
    using type = T;
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

        base_grad(const std::tuple_element_t<I, tuple_type> &...args) {
            (assign(std::get<I>(m_underlying), args), ...);
        }

        template <typename... U>
        base_grad(const std::tuple<U...> &other) : m_underlying(other) {}

        base_grad(const base_grad &other) = default;

        base_grad(base_grad &&other) = default;

        base_grad &operator=(const base_grad &other) {
            (assign(std::get<I>(m_underlying), std::get<I>(other.m_underlying)), ...);

            return *this;
        }

        base_grad &operator=(base_grad &&other) = default; // { return *this; }

        tuple_type &underlying_tuple() { return m_underlying; }

        const tuple_type &underlying_tuple() const { return m_underlying; }
    };

    template <size_t... I, typename T>
    class base_grad<std::index_sequence<I...>, T &> {
      public:
        using indices = std::index_sequence<I...>;
        using tuple_type = std::tuple<std::tuple_element_t<I, grad<T &, sizeof...(I) - 1>>...>;

      protected:
        tuple_type m_underlying;

      public:
        base_grad() = default;

        base_grad(const std::tuple_element_t<I, tuple_type> &...args) : m_underlying(args...) {}

        template <typename... U>
        base_grad(const std::tuple<U...> &other) : m_underlying(other) {}

        base_grad(const base_grad &other) = default;

        base_grad(const base_grad<indices, T> &other){};

        base_grad(base_grad &&other) = default;

        base_grad &operator=(base_grad &&other) = default; // { return *this; }
                                                           // { return *this; }

        tuple_type &underlying_tuple() { return m_underlying; }

        const tuple_type &underlying_tuple() const { return m_underlying; }

        base_grad &operator=(base_grad<indices, T> &other) { return *this; }

        base_grad &operator=(const base_grad &other) {
            (assign(std::get<I>(m_underlying), std::get<I>(other.m_underlying)), ...);

            return *this;
        }

        //        base_grad &operator=(const base_grad<indices, T> &other) {
        //            (assign(std::get<I>(m_underlying), std::get<I>(other.m_underlying)), ...);

        //          return *this;
        //    }
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
        std::apply(
            [&other](auto &...args) {
                std::apply(
                    [&args...](const auto &...other_args) { (detail::assign(args, other_args), ...); },
                    other.underlying_tuple()
                );
            },
            base::underlying_tuple()
        );

        return *this;
    }
};

template <typename T0, typename... T>
grad(T0 &, T &...) -> grad<T0 &, sizeof...(T)>;

template <size_t I, typename T, size_t N>
T &get(grad<T, N> &t) {
    return std::get<I>(t.underlying_tuple());
}

template <size_t I, typename T, size_t N>
const T &get(const grad<T, N> &t) {
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
