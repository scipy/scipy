#pragma once

#include <tuple>

#include "config.h"

namespace special {

template <size_t I, typename Res, size_t NArgs = 1>
struct grad_tuple_element;

template <size_t I, typename Res>
struct grad_tuple_element<I, Res, 1> {
    using type = Res;
};

template <size_t I, typename Res, size_t NArgs>
using grad_tuple_element_t = typename grad_tuple_element<I, Res, NArgs>::type;

namespace detail {

    template <typename Res, typename I, size_t NArgs>
    struct grad_tuple;

    template <typename Res, size_t... I, size_t NArgs>
    struct grad_tuple<Res, std::index_sequence<I...>, NArgs> {
        using type = std::tuple<grad_tuple_element_t<I, Res, NArgs>...>;
    };

} // namespace detail

template <typename Res, size_t N, size_t NArgs = 1>
using grad_tuple_t = typename detail::grad_tuple<Res, std::make_index_sequence<N + 1>, NArgs>::type;

namespace tuples {
    namespace detail {

        template <typename T, typename U>
        void assign(T &dst, const U &src) {
            dst = src;
        }

        template <typename T, size_t N, typename U>
        void assign(T (&dst)[N], const U (&src)[N]) {
            for (size_t i = 0; i < N; ++i) {
                assign(dst[i], src[i]);
            }
        }

        template <typename T, typename U>
        void fill(T &dst, const U &src) {
            dst = src;
        }

        template <typename T, size_t N, typename U>
        void fill(T (&dst)[N], const U &src) {
            for (size_t i = 0; i < N; ++i) {
                fill(dst[i], src);
            }
        }

    } // namespace detail

    template <typename... T>
    struct initializer_tuple {
        std::tuple<T...> underlying;

        initializer_tuple(const T &...args) {
            std::apply([&args...](auto &...other_args) { (detail::assign(other_args, args), ...); }, underlying);
        }
    };

    template <typename... T, typename Index>
    decltype(auto) access(std::tuple<T...> &t, Index &&index) {
        return std::apply([&index](auto &...args) { return std::tie(args[std::forward<Index>(index)]...); }, t);
    }

    template <typename... T>
    const std::tuple<T &...> &assign(const std::tuple<T &...> &t, const initializer_tuple<T...> &other) {
        std::apply(
            [&t](const auto &...other_args) {
                std::apply([&other_args...](auto &...args) { (detail::assign(args, other_args), ...); }, t);
            },
            other.underlying
        );

        return t;
    }

    template <typename... T, typename Arg>
    void fill(std::tuple<T...> t, Arg arg) {
        std::apply([arg](auto &...args) { (detail::fill(args, arg), ...); }, t);
    }

    template <typename... T, typename Func>
    void for_each(std::tuple<T...> t, Func f) {
        std::apply([f](auto &...args) { (f(args), ...); }, t);
    }

    template <typename... T>
    std::tuple<T &...> ref(std::tuple<T...> &t) {
        return std::apply([](auto &...args) { return std::tie(args...); }, t);
    }

    template <typename... T, typename... Args>
    std::tuple<std::result_of_t<T(Args...)>...> call(std::tuple<T...> &t, Args &&...args) {
        return std::apply(
            [&args...](auto &...elements) { return std::tie(elements(std::forward<Args>(args)...)...); }, t
        );
    }

    template <typename... T, typename... Args>
    std::tuple<std::result_of_t<T(Args...)>...> call(const std::tuple<T...> &t, Args &&...args) {
        return std::apply(
            [&args...](auto &...elements) { return std::tie(elements(std::forward<Args>(args)...)...); }, t
        );
    }

    template <typename... T, typename... U>
    std::tuple<std::tuple<T, U>...> zip(std::tuple<T...> t, std::tuple<U...> u) {
        return std::apply(
            [&t](auto &...other_args) {
                return std::apply(
                    [&other_args...](auto &...args) {
                        return std::make_tuple(std::forward_as_tuple(args, other_args)...);
                    },
                    t
                );
            },
            u
        );
    }

} // namespace tuples
} // namespace special
