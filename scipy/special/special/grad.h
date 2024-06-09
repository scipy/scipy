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

    template <typename T, size_t I>
    struct grad_leaf {
        T value;

      public:
        grad_leaf() = default;

        grad_leaf(const T &other) : value(other) {}
    };

    template <typename T, size_t I, size_t K>
    struct grad_leaf<T[K], I> {
        T value[K];

      public:
        grad_leaf() = default;

        grad_leaf(const T (&other)[K]) {
            for (size_t k = 0; k < K; ++k) {
                value[k] = other[k];
            }
        }
    };

    template <typename Indices, typename T>
    class grad_sequence;

    template <size_t... I, typename T>
    class grad_sequence<std::index_sequence<I...>, T> : public grad_leaf<T, I>... {
        using indices = std::index_sequence<I...>;
        static constexpr size_t N = sizeof...(I) - 1;

        //        std::tuple<T...> m_tup;

      public:
        grad_sequence() = default;

        grad_sequence(const std::tuple_element_t<I, grad<T, N>> &...args) : grad_leaf<T, I>(args)... {}

        template <typename... U>
        grad_sequence(const std::tuple<U...> &other) : grad_leaf<T, I>(std::get<I>(other))... {}

        grad_sequence(const grad_sequence &other) = default;

        grad_sequence(grad_sequence &&other) = default;

        grad_sequence &operator=(const grad_sequence &other) = default; // { return *this; }

        grad_sequence &operator=(grad_sequence &&other) = default; // { return *this; }

        //        T &value() { return static_cast<grad_leaf<T, 0> *>(this)->value; }

        //      const T &value() const { return static_cast<const grad_leaf<T, 0> *>(this)->value; }

        // underlying_tuple
        std::tuple<std::tuple_element_t<I, grad<T, N>>...> as_tuple() const {
            return std::forward_as_tuple(static_cast<const grad_leaf<T, I> *>(this)->value...);
        }

        grad_sequence<indices, T &> refs() { return {static_cast<grad_leaf<T, I> *>(this)->value...}; }

        grad_sequence<indices, const T &> refs() const { return crefs(); }

        grad_sequence<indices, const T &> crefs() const {
            return {static_cast<const grad_leaf<T, I> *>(this)->value...};
        }
    };

} // namespace detail

template <typename T, size_t N>
class grad : public detail::grad_sequence<std::make_index_sequence<N + 1>, T> {
    using base = detail::grad_sequence<std::make_index_sequence<N + 1>, T>;

  public:
    grad(const base &other) : base(other) {}

  public:
    using base::base;

    grad<T &, N> refs() { return base::refs(); }

    grad<const T &, N> refs() const { return base::refs(); }

    grad<const T &, N> crefs() const { return base::refs(); }

    decltype(auto) refs_as_tuple() {
        const auto &tmp = base::refs();
        return tmp.as_tuple();
    }

    decltype(auto) refs_as_tuple() const { return crefs_as_tuple(); }

    decltype(auto) crefs_as_tuple() const {
        const auto &tmp = base::crefs();
        return tmp.as_tuple();
    }
};

/*
template <typename T, size_t N>
class grad<T &, N> : public detail::grad_tuple<T &, std::make_index_sequence<N + 1>> {
  public:
    using detail::grad_tuple<T &, std::make_index_sequence<N + 1>>::grad_tuple;

    grad &operator=(const grad<T, N> &other) { return *this; }

};
*/

template <size_t I, typename T, size_t N>
T &get(grad<T, N> &t) {
    return static_cast<detail::grad_leaf<T, I> *>(&t)->value;
}

template <size_t I, typename T, size_t N>
const T &get(const grad<T, N> &t) {
    return static_cast<const detail::grad_leaf<T, I> *>(&t)->value;
}

template <typename T, size_t K>
void dot(const grad<T[K], 0> &x, const grad<T[K], 0> &y, grad<T, 0> &res) {
    const auto &[x0] = x;
    const auto &[y0] = y;
    auto &[res0] = res;

    res0 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
    }
}

template <typename T, size_t K>
void dot(const grad<T[K], 1> &x, const grad<T[K], 1> &y, grad<T, 1> &res) {
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
void dot(const grad<T[K], 2> &x, const grad<T[K], 2> &y, grad<T, 2> &res) {
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
