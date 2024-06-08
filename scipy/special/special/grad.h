#pragma once

#include <tuple>

#include "config.h"

namespace special {
namespace detail {

    template <typename T, size_t I>
    struct grad_tuple_leaf {
        using value_type = T;

        T value;

      public:
        grad_tuple_leaf() = default;

        grad_tuple_leaf(T other) : value(other) {}
    };

    template <typename T, size_t I, size_t K>
    struct grad_tuple_leaf<T[K], I> {
        using value_type = T[K];

        T value[K];

      public:
        grad_tuple_leaf() = default;

        grad_tuple_leaf(T (&other)[K]) {}

        grad_tuple_leaf &operator=(const T (&other)[K]) {
            for (size_t k = 0; k < K; ++k) {
                value[k] = other[k];
            }

            return *this;
        }
    };

    template <typename T, typename Indices>
    class grad_tuple;

    template <typename T, size_t... I>
    class grad_tuple<T, std::index_sequence<I...>> : public grad_tuple_leaf<T, I>... {
      public:
        grad_tuple() = default;

        template <typename... Args>
        explicit grad_tuple(Args &&...args) : grad_tuple_leaf<T, I>(std::forward<Args>(args))... {}

        std::tuple<typename grad_tuple_leaf<T, I>::value_type &...> refs() {
            return std::tie(static_cast<grad_tuple_leaf<T, I> *>(this)->value...);
        }

        std::tuple<const typename grad_tuple_leaf<T, I>::value_type &...> refs() const { return crefs(); }

        std::tuple<const typename grad_tuple_leaf<T, I>::value_type &...> crefs() const {
            return std::tie(static_cast<const grad_tuple_leaf<T, I> *>(this)->value...);
        }

        template <typename... Args>
        grad_tuple &operator=(const std::tuple<Args...> &other) {
            ((static_cast<grad_tuple_leaf<T, I> *>(this)->value = std::get<I>(other)), ...);

            return *this;
        }

        grad_tuple &operator=(std::initializer_list<T> other) {
            ((*static_cast<grad_tuple_leaf<T, I> *>(this) = *(other.begin() + I)), ...);

            return *this;
        }
    };

} // namespace detail

template <typename T, size_t N>
class grad_tuple : public detail::grad_tuple<T, std::make_index_sequence<N + 1>> {
  public:
    using detail::grad_tuple<T, std::make_index_sequence<N + 1>>::grad_tuple;

    using detail::grad_tuple<T, std::make_index_sequence<N + 1>>::operator=;
};

template <size_t I, typename T, size_t N>
T &get(grad_tuple<T, N> &t) {
    return static_cast<detail::grad_tuple_leaf<T, I> *>(&t)->value;
}

template <size_t I, typename T, size_t N>
const T &get(const grad_tuple<T, N> &t) {
    return static_cast<const detail::grad_tuple_leaf<T, I> *>(&t)->value;
}

template <typename T, size_t K>
void dot(const grad_tuple<T[K], 0> &x, const grad_tuple<T[K], 0> &y, grad_tuple<T, 0> &res) {
    const auto &[x0] = x;
    const auto &[y0] = y;
    auto &[res0] = res;

    res0 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
    }
}

template <typename T, size_t K>
void dot(const grad_tuple<T[K], 1> &x, const grad_tuple<T[K], 1> &y, grad_tuple<T, 1> &res) {
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
void dot(const grad_tuple<T[K], 2> &x, const grad_tuple<T[K], 2> &y, grad_tuple<T, 2> &res) {
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

template <size_t I, typename T, size_t N>
struct tuple_element<I, special::grad_tuple<T, N>> {
    using type = T;
};

template <typename T, size_t N>
struct tuple_size<special::grad_tuple<T, N>> : std::integral_constant<size_t, N + 1> {};

} // namespace std
