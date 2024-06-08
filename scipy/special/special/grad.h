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
    };

} // namespace detail

template <typename T, size_t N>
class grad_tuple : public detail::grad_tuple<T, std::make_index_sequence<N + 1>> {
  public:
    using detail::grad_tuple<T, std::make_index_sequence<N + 1>>::grad_tuple;
};

template <size_t I, typename T, size_t N>
T &get(grad_tuple<T, N> &t) {
    return static_cast<detail::grad_tuple_leaf<T, I> *>(&t)->value;
}

template <size_t I, typename T, size_t N>
const T &get(const grad_tuple<T, N> &t) {
    return static_cast<const detail::grad_tuple_leaf<T, I> *>(&t)->value;
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
