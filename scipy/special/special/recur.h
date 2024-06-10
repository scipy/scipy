#pragma once

#include "config.h"
#include "tuple_wrapper.h"

namespace special {

template <typename T, size_t K>
void forward_recur_shift_left(T (&res)[K]) {
    for (size_t k = 1; k < K; ++k) {
        res[k - 1] = res[k];
    }
}

template <typename... T, size_t K>
void forward_recur_shift_left(tuple_wrapper<T (&)[K]...> &res) {
    std::apply([](auto &...args) { (forward_recur_shift_left(args), ...); }, res.underlying_tuple());
}

template <typename T, size_t K>
void forward_recur_rotate_left(T (&res)[K]) {
    T tmp = res[0];
    forward_recur_shift_left(res);
    res[K - 1] = tmp;
}

template <typename... T, size_t K>
void forward_recur_rotate_left(tuple_wrapper<T (&)[K]...> &res) {
    std::apply([](auto &...args) { (forward_recur_rotate_left(args), ...); }, res.underlying_tuple());
}

/**
 * Compute a forward recurrence that depends on K previous values.
 *
 * @param first begin iterator
 * @param last end iterator
 * @param r recurrence
 * @param res values, initialised to the leading K values
 * @param f a function to be called as f(it, res)
 */
template <typename InputIt, typename Recurrence, typename... T, ssize_t K, typename Func>
void forward_recur(InputIt first, InputIt last, Recurrence r, tuple_wrapper<T (&)[K]...> res, Func f) {
    InputIt it = first;
    while (it - first != K && it != last) {
        forward_recur_rotate_left(res);

        f(it, res);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            tuple_wrapper<T[K]...> coef;
            r(it, apply_tie(coef));

            tuple_wrapper<T...> tmp;
            dot(apply_tie(coef), res, apply_tie(tmp));

            std::apply([](auto &...args) { (forward_recur_shift_left(args), ...); }, res.underlying_tuple());
            std::apply([](auto &...args) { return std::tie(args[K - 1]...); }, res.underlying_tuple()) =
                tmp.underlying_tuple();

            f(it, res);
            ++it;
        }
    }
}

} // namespace special
